"""
download.py — fetch planned units from the open S3 buckets, verified.

Not a step. Called by `index.py get` AFTER plan.py has printed the
preview and the --max-gb check has passed. Nothing here decides WHAT to
fetch; it only fetches what the plan lists, and proves each file arrived
intact.

    python scripts/index.py get LEC --product metadata --match final-nuclei

Decisions made here, and why:

1. ASK BEFORE EVERY DOWNLOAD, UNLESS TOLD NOT TO.
   --confirm-gb defaults to 0: any real transfer, however small, prints
   the plan and asks y/N. --yes skips the question (for scripts). With
   no terminal to ask and no --yes, it refuses rather than assuming.

2. THE FILE MUST STILL BE WHAT THE CATALOG SAID.
   Immediately before transferring, HEAD the object. If its size or ETag
   differs from the catalog row, the file changed on S3 after the last
   refresh: stop, say so, and fetch nothing. Downloading a different
   file than the one previewed would make the preview a lie.

3. WRITE TO .part, VERIFY, THEN RENAME.
   Bytes stream into <file>.part. Only after every check passes is it
   renamed to its real name (an atomic rename on the same disk). A
   dropped connection, a wrong size or a wrong checksum deletes the
   .part file, so a truncated file can never sit under the real name —
   the failure the h5ad_path size check once caught by luck.

4. VERIFICATION IS AS STRONG AS S3 ALLOWS, AND SAYS WHICH.
   These buckets publish no checksums, so the ETag is the only one:
   - single-part upload: the ETag IS the file's MD5 -> 'md5 match';
   - multipart upload ('<hash>-N'): the ETag is the MD5 of the parts'
     MD5s. HEAD with PartNumber=1 returns part 1's exact size, so the
     ETag is recomputed exactly -> 'multipart etag match';
   - if part size cannot be learned: 'size only', recorded as such.
   The wording in the ledger never claims more than was checked.

5. A SHA-256 OF EVERY FILE GOES IN A LOCAL LEDGER.
   data_cache/_ledger.csv records time, bucket, key, size, ETag, S3
   last-modified, SHA-256 and the verification level. It answers "which
   bytes did I analyse" later, and re-checking a file needs no network.
   The path is written relative to data_cache/ (where the ledger lives), so
   the ledger survives moving or renaming the project folder and holds no
   user name. Older rows with absolute paths are still read.

6. FILES ONLY, FOR NOW.
   Collapsed units (a Zarr store, 1.1 million plot images) are refused
   with a message; fetching folders of thousands of objects needs its
   own safeguards and comes later.

7. NOTHING IS OPENED OR PRINTED.
   Bytes go from the network to disk through two hash functions. No
   data file content is ever displayed.

8. THERE MUST BE ROOM.
   Free disk space is checked against the total (plus 5%) before the
   first byte is fetched.
"""

import csv
import datetime
import hashlib
import math
import pathlib
import shutil
import sys
import time

import pathrules
import sources

CHUNK = 8 * 1024 * 1024
LEDGER_FIELDS = ["downloaded_at", "bucket", "key", "size", "etag",
                 "s3_last_modified", "sha256", "verified", "path"]


class DownloadError(RuntimeError):
    pass


def human(n):
    for unit, div in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= div:
            return f"{n / div:,.2f} {unit}"
    return f"{n} B"


# ---- checksums ------------------------------------------------------------

class ETagHasher:
    """Recomputes an S3 ETag while bytes stream past.

    part_size=None -> single-part: plain MD5 of the whole file.
    part_size=P    -> multipart: MD5 of the concatenated per-part MD5
                      digests, suffixed with '-<number of parts>'.
    """

    def __init__(self, part_size=None):
        self.part_size = part_size
        self.whole = hashlib.md5()
        self.part = hashlib.md5()
        self.in_part = 0
        self.digests = []

    def update(self, data):
        if self.part_size is None:
            self.whole.update(data)
            return
        view = memoryview(data)
        while view:
            take = min(len(view), self.part_size - self.in_part)
            self.part.update(view[:take])
            self.in_part += take
            view = view[take:]
            if self.in_part == self.part_size:
                self.digests.append(self.part.digest())
                self.part, self.in_part = hashlib.md5(), 0

    def hexdigest(self):
        if self.part_size is None:
            return self.whole.hexdigest()
        digests = self.digests + ([self.part.digest()] if self.in_part else [])
        return f"{hashlib.md5(b''.join(digests)).hexdigest()}-{len(digests)}"


def multipart_count(etag):
    """'abc-95' -> 95; 'abc' -> None."""
    if "-" in etag:
        tail = etag.rsplit("-", 1)[1]
        if tail.isdigit():
            return int(tail)
    return None


def part_size_for(client, bucket, key, etag, size):
    """Exact part size of a multipart object, or None if it cannot be learned."""
    n = multipart_count(etag)
    if n is None:
        return None
    try:
        head = client.head_object(Bucket=bucket, Key=key, PartNumber=1)
        p = int(head["ContentLength"])
    except Exception:
        return None
    # A real part size reproduces the part count; anything else is suspect.
    return p if p > 0 and math.ceil(size / p) == n else None


# ---- one file -------------------------------------------------------------

def fetch_one(client, row, dest, progress=None):
    """Download, verify and place one file. Returns a ledger dict."""
    bucket, key = row["bucket"], row["unit"]
    want_size, want_etag = int(row["size"]), row["etag"]

    head = client.head_object(Bucket=bucket, Key=key)
    now_size = int(head["ContentLength"])
    now_etag = head.get("ETag", "").strip('"')
    if now_size != want_size or (want_etag and now_etag != want_etag):
        raise DownloadError(
            f"{key} changed on S3 since the catalog was built "
            f"(catalog: {want_size:,} bytes, ETag {want_etag}; now: {now_size:,} "
            f"bytes, ETag {now_etag}). Nothing was downloaded. Refresh the "
            f"catalog, re-check the plan, then try again.")

    part_size = part_size_for(client, bucket, key, now_etag, now_size)
    multipart = multipart_count(now_etag) is not None
    hasher = ETagHasher(part_size) if (part_size or not multipart) else None
    sha = hashlib.sha256()

    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    got = 0
    try:
        body = client.get_object(Bucket=bucket, Key=key)["Body"]
        with open(part, "wb") as fh:
            while True:
                data = body.read(CHUNK)
                if not data:
                    break
                fh.write(data)
                sha.update(data)
                if hasher:
                    hasher.update(data)
                got += len(data)
                if progress:
                    progress(got, want_size)
        if got != want_size:
            raise DownloadError(
                f"{key}: received {got:,} bytes, expected {want_size:,} "
                f"(connection dropped?). The partial file was deleted.")
        if hasher:
            calc = hasher.hexdigest()
            if calc != now_etag:
                raise DownloadError(
                    f"{key}: checksum mismatch (computed {calc}, S3 ETag "
                    f"{now_etag}). The file was deleted; do not use it.")
            verified = "multipart etag match" if multipart else "md5 match"
        else:
            verified = "size only"
        part.replace(dest)
    except BaseException:
        part.unlink(missing_ok=True)
        raise

    return {"downloaded_at": datetime.datetime.now(datetime.timezone.utc)
                                     .isoformat(timespec="seconds"),
            "bucket": bucket, "key": key, "size": got, "etag": now_etag,
            "s3_last_modified": str(head.get("LastModified", "")),
            "sha256": sha.hexdigest(), "verified": verified, "path": str(dest)}


def ledger_path(dest, data_dir):
    """How a downloaded file is written in the ledger: relative to data_cache/.

    The ledger lives in data_cache/, so a relative path keeps working when the
    project folder is moved or renamed, and carries no user name or home folder.
    """
    try:
        return pathlib.Path(dest).resolve().relative_to(pathlib.Path(data_dir).resolve()) \
            .as_posix()
    except ValueError:
        return str(dest)


def resolve_ledger_path(p, data_dir):
    """The file a ledger row points to, for rows written either way.

    Rows written before paths were relative hold an absolute path. If that path
    no longer exists (the folder was moved), the part after the last
    'data_cache' folder is looked up under the current data_cache/.
    """
    p = pathlib.Path(p)
    if not p.is_absolute():
        return pathlib.Path(data_dir) / p
    if not p.exists() and "data_cache" in p.parts:
        i = len(p.parts) - 1 - p.parts[::-1].index("data_cache")
        moved = pathlib.Path(data_dir).joinpath(*p.parts[i + 1:])
        if moved.exists():
            return moved
    return p


def read_ledger(data_dir):
    """Ledger rows, each with 'file': the resolved local path. [] if no ledger."""
    led = pathlib.Path(data_dir) / "_ledger.csv"
    if not led.exists():
        return []
    with open(led, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["file"] = resolve_ledger_path(r["path"], data_dir)
    return rows


def append_ledger(path, entry):
    new = not path.exists()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LEDGER_FIELDS)
        if new:
            w.writeheader()
        w.writerow(entry)


# ---- a whole plan ---------------------------------------------------------

def run(plan, data_dir, client=None, yes=False, force=False, confirm_gb=0.0,
        ask=input, isatty=None, out=print, disk_usage=shutil.disk_usage, show_progress=None):
    """Fetch every unit in a plan built by plan.build(). Returns exit code.

    The percentage meter is for a person watching: by default it shows only when
    stderr is a terminal and the caller has not redirected `out` (tests always do).
    """
    if show_progress is None:
        show_progress = out is print and sys.stderr.isatty()
    units = plan["units"]
    bad_bucket = sorted({u["bucket"] for u in units} - set(pathrules.BUCKETS))
    if bad_bucket:
        out(f"\n  REFUSED: {', '.join(bad_bucket)} is not an open bucket this tool "
            f"may download from.")
        return 4
    folders = [u for u in units if u["kind"] != "file"]
    if folders:
        out(f"\n  REFUSED: {len(folders)} collapsed folder unit(s) in the selection "
            f"(e.g. {folders[0]['unit']}). Downloading folders of many objects "
            f"is not supported yet; narrow the selection to files.")
        return 4

    todo = [u for u in units if force or u["state"] != "present"]
    if not todo:
        out("\n  Nothing to download: every selected file is already present at "
            "the catalog's size. --force re-downloads.")
        return 0
    need = sum(int(u["size"]) for u in todo)

    free = disk_usage(data_dir if data_dir.exists() else data_dir.parent).free
    if need * 1.05 > free:
        out(f"\n  REFUSED: {human(need)} needed (+5% margin) but only {human(free)} "
            f"free on the disk holding {data_dir}.")
        return 5

    if need > confirm_gb * 1e9 and not yes:
        if not (isatty if isatty is not None else sys.stdin.isatty()):
            out(f"\n  {human(need)} needs your confirmation, but there is no terminal "
                f"to ask. Re-run with --yes to confirm.")
            return 6
        reply = ask(f"\n  Download {len(todo)} file(s), {human(need)}? [y/N] ")
        if reply.strip().lower() not in ("y", "yes"):
            out("\n  Aborted. Nothing was downloaded.")
            return 1

    client = client or sources.anonymous_client()
    ledger = data_dir / "_ledger.csv"
    out(f"\n{'=' * 78}\nDOWNLOAD\n{'=' * 78}")
    done, t0 = 0, time.time()
    for i, u in enumerate(todo, 1):
        out(f"  [{i}/{len(todo)}] {human(int(u['size'])):>10s}  s3://{u['bucket']}/{u['unit']}")
        last = [0.0]

        def progress(got, total, last=last):
            if not show_progress:
                return          # progress is for a person watching, not for logs
            if total and (got / total - last[0] >= 0.1 or got == total):
                last[0] = got / total
                print(f"      {100 * got / total:5.1f}%  {human(got)}", file=sys.stderr,
                      flush=True)
        try:
            entry = fetch_one(client, u, u["dest"], progress)
        except DownloadError as exc:
            out(f"\n  FAILED: {exc}")
            out(f"  {i - 1} of {len(todo)} file(s) completed before this one.")
            return 7
        entry["path"] = ledger_path(u["dest"], data_dir)
        append_ledger(ledger, entry)
        done += entry["size"]
        out(f"      ok — {entry['verified']}; sha256 {entry['sha256'][:16]}…")
        out(f"      -> {u['dest']}")

    secs = max(time.time() - t0, 1e-6)
    out(f"\n  downloaded: {len(todo)} file(s), {human(done)} in {secs:.0f} s "
        f"({human(done / secs)}/s)")
    out(f"  ledger:     {ledger}")
    return 0
