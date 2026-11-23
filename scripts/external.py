"""
external.py — sources that are not S3 buckets: open web files and controlled-access studies.

Not a step. Run through index.py:

    python scripts/index.py sources                 # list them (no network)
    python scripts/index.py sources --check         # HEAD each open file: size, date
    python scripts/index.py get --external web-donor-metadata --dry-run

The list itself is catalog/external_sources.csv, written by hand from the
SEA-AD data page and reviewed like code. The latest check of each open
file is catalog/external_status.csv.

Decisions made here, and why:

1. A HAND-WRITTEN LIST, BECAUSE A WEBSITE CANNOT BE LISTED.
   An S3 bucket can be listed; a web page cannot, reliably. So each
   source is one reviewed row: what it is, its exact URL, open or
   controlled, and the page it was listed on. Adding a source is a
   one-line change anyone can review.

2. A HEAD REQUEST IS THE WEB'S 'LISTING'.
   `sources --check` asks each open URL for its headers only (size, last
   modified, ETag, content type). No content is transferred, so, like
   listing a bucket, it needs no approval. It sends a generic client name
   and nothing about you.

3. CONTROLLED ACCESS IS LINKED, NEVER FETCHED.
   Rows marked 'controlled' (the two Synapse studies) carry a URL to
   request access. `get` refuses them by design, whatever flags are given.
   The tool has no Synapse client and stores no credentials.

4. OPEN FILES GO THROUGH THE SAME GUARDS AS S3 FILES.
   Preview with the size from a fresh HEAD, --max-gb, y/N or --yes,
   a .part file renamed only after checks pass, a SHA-256 in the same
   ledger. The size must match the HEAD answer. Verification wording is
   honest: a web server's ETag is not promised to be an MD5, so the
   ledger says 'md5 matches ETag' only when it does, else 'size only'.

5. FILES LAND IN data_cache/web/<host>/<file name>.
   data_cache/ is already gitignored, so nothing downloaded can be
   committed by accident.
"""

import csv
import datetime
import hashlib
import pathlib
import re
import sys
import urllib.parse
import urllib.request

import index
import plan

ROOT = pathlib.Path(__file__).resolve().parent.parent
SOURCES_FILE = ROOT / "catalog" / "external_sources.csv"
STATUS_FILE = ROOT / "catalog" / "external_status.csv"
FIELDS = ["id", "title", "access", "url", "format", "modality", "regions", "product",
          "provider", "listed_on", "notes"]
STATUS_FIELDS = ["id", "size", "etag", "last_modified", "content_type", "checked"]
ACCESS = ("open", "controlled")
USER_AGENT = "sea-ad-catalog/1.0 (structure-only data indexer)"
CHUNK = 1024 * 1024


class SourceError(ValueError):
    pass


# ---- the list -------------------------------------------------------------------

def load(path=None):
    path = pathlib.Path(path or SOURCES_FILE)
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != FIELDS:
            raise SourceError(f"{path.name}: columns must be {', '.join(FIELDS)}")
        rows = list(reader)
    seen = set()
    for n, r in enumerate(rows, start=2):
        where = f"{path.name} line {n} ({r['id'] or 'no id'})"
        if not r["id"] or r["id"] in seen:
            raise SourceError(f"{where}: id missing or repeated")
        seen.add(r["id"])
        if r["access"] not in ACCESS:
            raise SourceError(f"{where}: access must be open or controlled")
        if not r["url"].startswith(("https://", "http://")):
            raise SourceError(f"{where}: url must start with https:// or http://")
        if r["access"] == "open" and not r["url"].startswith("https://"):
            raise SourceError(f"{where}: open files must be fetched over https")
        if r["access"] == "controlled" and r["id"].startswith("syn") and \
                not re.fullmatch(r"syn\d+", r["id"]):
            raise SourceError(f"{where}: a Synapse id looks like syn12345678")
    return rows


def load_status(path=None):
    path = pathlib.Path(path or STATUS_FILE)
    if not path.exists():
        return {}
    with open(path, newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}


def destination(row, data_dir):
    u = urllib.parse.urlparse(row["url"])
    return data_dir / "web" / u.netloc / pathlib.PurePosixPath(u.path).name


# ---- HEAD -------------------------------------------------------------------------

def head(url, opener=None):
    """Headers of a URL without its content: {size, etag, last_modified, content_type}."""
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    with (opener or urllib.request.urlopen)(req, timeout=30) as resp:
        h = resp.headers
        size = h.get("Content-Length")
        return {"size": int(size) if size and size.isdigit() else None,
                "etag": (h.get("ETag") or "").strip('"'),
                "last_modified": h.get("Last-Modified") or "",
                "content_type": h.get("Content-Type") or ""}


def check(sources_path=None, status_path=None, opener=None, out=print):
    rows = [r for r in load(sources_path) if r["access"] == "open"]
    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    status, problems = load_status(status_path), 0
    out(f"\n{'=' * 78}\nCHECK (HEAD only: headers, no content)\n{'=' * 78}")
    for r in rows:
        try:
            h = head(r["url"], opener)
        except Exception as exc:                 # report every failure, keep going
            problems += 1
            out(f"  FAILED  {r['id']}: {type(exc).__name__}: {exc}")
            continue
        old = status.get(r["id"])
        changed = old and (str(old["size"]) != str(h["size"]) or old["etag"] != h["etag"])
        status[r["id"]] = {"id": r["id"], "size": h["size"] if h["size"] is not None else "",
                           "etag": h["etag"], "last_modified": h["last_modified"],
                           "content_type": h["content_type"], "checked": now}
        size = index.human(h["size"]) if h["size"] is not None else "size unknown"
        out(f"  ok      {r['id']:24s} {size:>11s}  {h['last_modified']}"
            + ("   CHANGED since last check" if changed else ""))
    path = pathlib.Path(status_path or STATUS_FILE)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=STATUS_FIELDS)
        w.writeheader()
        for k in sorted(status):
            w.writerow(status[k])
    out(f"\n  wrote {path}")
    return 1 if problems else 0


# ---- listing ------------------------------------------------------------------------

def report(sources_path=None, status_path=None, out=print):
    rows, status = load(sources_path), load_status(status_path)
    out(f"\n{'=' * 78}\nOTHER SOURCES (not in the S3 buckets)\n{'=' * 78}")
    for access in ACCESS:
        block = [r for r in rows if r["access"] == access]
        if not block:
            continue
        label = ("OPEN: direct downloads (get --external <id>)" if access == "open" else
                 "CONTROLLED ACCESS: linked only, never downloaded by this tool")
        out(f"\n  {label}")
        for r in block:
            st = status.get(r["id"], {})
            size = ("—" if access == "controlled" else
                    index.human(int(st["size"])) if st.get("size") else "not checked")
            out(f"    {r['id']:24s} {size:>11s}  {r['title']}"
                + (f"  [{r['regions']}]" if r["regions"] else ""))
            if access == "controlled":
                out(f"    {'':24s} {'':>11s}  request access: {r['url']}")
    checked = sorted({s["checked"] for s in status.values()})
    out(f"\n  sizes from: {checked[-1] if checked else 'never checked (sources --check)'}")
    return rows


# ---- download -----------------------------------------------------------------------

def run_get(pattern, data_dir=None, max_gb=50.0, dry_run=False, yes=False, force=False,
            sources_path=None, opener=None, ask=input, isatty=None, out=print):
    import fnmatch
    import download
    data_dir = data_dir or plan.DATA_DIR
    rows = load(sources_path)
    chosen = [r for r in rows if fnmatch.fnmatch(r["id"], pattern)]
    if not chosen:
        raise SystemExit(f"\nNo source matches {pattern!r}. See: python scripts/index.py sources")
    controlled = [r for r in chosen if r["access"] == "controlled"]
    if controlled:
        c = controlled[0]
        raise SystemExit(f"\nREFUSED: {c['id']} ({c['title']}) is controlled-access data. This "
                         f"tool links to it but never downloads it. Request access at:\n  "
                         f"{c['url']}")

    out(f"\n{'=' * 78}\nDOWNLOAD PLAN — open web files (sizes from a fresh HEAD)\n{'=' * 78}")
    todo, need = [], 0
    for r in chosen:
        h = head(r["url"], opener)
        dest = destination(r, data_dir)
        have = dest.stat().st_size if dest.exists() else None
        state = ("present" if h["size"] is not None and have == h["size"] else
                 "partial" if have is not None else "missing")
        out(f"  {index.human(h['size'] or 0):>11s}  {state:8s} {r['id']}  ->  {dest}")
        if force or state != "present":
            todo.append((r, h, dest))
            need += h["size"] or 0
    out(f"\n  TO DOWNLOAD: {len(todo)} file(s), {index.human(need)}")
    if need > max_gb * 1e9:
        out(f"\n  REFUSED: exceeds --max-gb {max_gb:g}")
        return 2
    if dry_run:
        out("\n  --dry-run: nothing was transferred.")
        return 0
    if not todo:
        out("\n  Nothing to download.")
        return 0
    if not yes:
        if not (isatty if isatty is not None else sys.stdin.isatty()):
            out(f"\n  {index.human(need)} needs your confirmation, but there is no terminal "
                f"to ask. Re-run with --yes to confirm.")
            return 6
        if ask(f"\n  Download {len(todo)} file(s), {index.human(need)}? [y/N] ").strip() \
                .lower() not in ("y", "yes"):
            out("\n  Aborted. Nothing was downloaded.")
            return 1
    ledger = data_dir / "_ledger.csv"
    for r, h, dest in todo:
        try:
            entry = fetch(r, h, dest, opener)
        except download.DownloadError as exc:
            out(f"\n  FAILED: {exc}")
            return 7
        entry["path"] = download.ledger_path(dest, data_dir)
        download.append_ledger(ledger, entry)
        out(f"  ok — {entry['verified']}; sha256 {entry['sha256'][:16]}…  -> {dest}")
    return 0


def fetch(row, h, dest, opener=None):
    """Stream to .part, check size (and MD5 if the ETag is one), rename, ledger entry."""
    import download
    req = urllib.request.Request(row["url"], headers={"User-Agent": USER_AGENT})
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    sha, md5, got = hashlib.sha256(), hashlib.md5(), 0
    try:
        with (opener or urllib.request.urlopen)(req, timeout=60) as resp, open(part, "wb") as fh:
            while True:
                data = resp.read(CHUNK)
                if not data:
                    break
                fh.write(data)
                sha.update(data)
                md5.update(data)
                got += len(data)
        if h["size"] is not None and got != h["size"]:
            raise download.DownloadError(
                f"{row['id']}: received {got:,} bytes, the server announced {h['size']:,}. "
                f"The partial file was deleted.")
        etag = h["etag"]
        if re.fullmatch(r"[0-9a-f]{32}", etag or ""):
            verified = "md5 matches ETag" if md5.hexdigest() == etag else "size only"
        else:
            verified = "size only"
        part.replace(dest)
    except BaseException:
        part.unlink(missing_ok=True)
        raise
    u = urllib.parse.urlparse(row["url"])
    return {"downloaded_at": datetime.datetime.now(datetime.timezone.utc)
                                     .isoformat(timespec="seconds"),
            "bucket": "web", "key": f"{u.netloc}{u.path}", "size": got, "etag": h["etag"],
            "s3_last_modified": h["last_modified"], "sha256": sha.hexdigest(),
            "verified": verified, "path": str(dest)}
