"""
peek.py — read an image's header straight from S3, without downloading the image.

Not a step. Run through index.py:

    python scripts/index.py peek --match "1193169898/DAPI_Max.tif" --dry-run
    python scripts/index.py peek --match "1193169898/DAPI_Max.tif"

A spatial-transcriptomics image is 1–7 GB; its header is a few kilobytes.
peek asks S3 for just those bytes ("range requests"), parses them with the
same tested code inspect uses for files on disk (images.py), and writes
the same report to results/inspect/. Nothing is saved to data_cache/.

Decisions made here, and why:

1. A PARTIAL READ IS STILL A DOWNLOAD, SO IT GETS THE SAME GUARDS.
   Your rule: partial reads only on files you approved. peek prints a
   preview (which files, their size, the most it will read from each) and
   asks y/N; --yes confirms for scripts; with no terminal and no --yes it
   refuses. --dry-run shows the preview and stops.

2. A HARD BYTE LIMIT PER FILE.
   Reads stop at --max-kb (default 4,096 KB = 4 MB) per file. A header
   that needs more (a huge OME-XML description, a file with thousands of
   pages) is reported as incomplete instead of quietly fetching more.

3. THE FILE MUST STILL BE WHAT THE CATALOG SAID, FOR EVERY BYTE.
   HEAD first: size and ETag must match the catalog. Each range request
   then carries If-Match: <ETag>, so if the object is replaced halfway
   through, S3 refuses rather than mixing bytes from two versions.

4. ONLY FORMATS WITH A HEADER READER.
   TIFF-family files (.tif, .tiff, .svs, .ome.tif). Anything else is
   refused with a message: peeking at a CSV's first bytes would be
   reading its contents.

5. EVERY PEEK IS LOGGED.
   data_cache/_peek_log.csv records time, bucket, key, ETag, bytes read and
   number of requests, so "which bytes did this report come from" has an
   answer even though no file was kept.

6. AT MOST 20 FILES AT A TIME.
   Peeking is for understanding a file type, not for sweeping a bucket.
"""

import csv
import datetime
import sys

import images
import index
import inspect_data
import plan

MAX_FILES = 20
BLOCK = 64 * 1024
PEEKABLE = ("tiff", "svs")
LOG_FIELDS = ["peeked_at", "bucket", "key", "etag", "size", "bytes_read", "requests"]


class PeekLimit(images.HeaderError):
    pass


class RangeReader:
    """read(offset, n) answered by S3 range requests, in 64 KB blocks, capped."""

    def __init__(self, client, bucket, key, size, etag, max_bytes):
        self.client, self.bucket, self.key = client, bucket, key
        self.size, self.etag, self.max_bytes = size, etag, max_bytes
        self.blocks, self.bytes_read, self.requests = {}, 0, 0

    def _block(self, i):
        if i not in self.blocks:
            start = i * BLOCK
            end = min(start + BLOCK, self.size) - 1
            if self.bytes_read + (end - start + 1) > self.max_bytes:
                raise PeekLimit(f"header needs more than the {self.max_bytes // 1024:,} KB "
                                f"limit (--max-kb)")
            kw = {"IfMatch": f'"{self.etag}"'} if self.etag else {}
            body = self.client.get_object(Bucket=self.bucket, Key=self.key,
                                          Range=f"bytes={start}-{end}", **kw)["Body"].read()
            self.requests += 1
            self.bytes_read += len(body)
            self.blocks[i] = body
        return self.blocks[i]

    def read(self, offset, n):
        if offset < 0 or offset >= self.size:
            raise images.HeaderError(f"offset {offset} outside the file")
        n = min(n, self.size - offset)
        out = bytearray()
        pos = offset
        while len(out) < n:
            i = pos // BLOCK
            block = self._block(i)
            start = pos - i * BLOCK
            take = block[start:start + (n - len(out))]
            if not take:
                break
            out += take
            pos += len(take)
        return bytes(out)


def select(match, bucket=None, catalog_dir=None, out=print):
    rows, _ = plan.select(None, match, None, None, bucket, False, False,
                          catalog_dir or index.CATALOG_DIR)
    files = [r for r in rows if r["kind"] == "file"]
    bad = [r for r in files if r["format"] not in PEEKABLE]
    if bad:
        raise SystemExit(f"\nREFUSED: peek reads image headers only (formats: "
                         f"{', '.join(PEEKABLE)}). Not peekable: {bad[0]['unit']} "
                         f"({bad[0]['format']}). Narrow --match.")
    if not files:
        raise SystemExit("\nNothing to peek at: no file matches.")
    if len(files) > MAX_FILES:
        raise SystemExit(f"\nREFUSED: {len(files)} files match; peek takes at most "
                         f"{MAX_FILES}. Narrow --match.")
    return files


def run(match, bucket=None, max_kb=4096, dry_run=False, yes=False, client=None,
        catalog_dir=None, data_dir=None, out_dir=None, ask=input, isatty=None, out=print):
    files = select(match, bucket, catalog_dir, out)
    cap = max_kb * 1024
    out(f"\n{'=' * 78}\nPEEK — read image headers from S3; nothing is saved\n{'=' * 78}")
    for r in files:
        out(f"  {index.human(r['size']):>11s}  s3://{r['bucket']}/{r['unit']}")
    out(f"\n  {len(files)} file(s); at most {max_kb:,} KB read from each "
        f"(≤ {index.human(cap * len(files))} in total)")
    if dry_run:
        out("\n  --dry-run: nothing was read.")
        return 0
    if not yes:
        if not (isatty if isatty is not None else sys.stdin.isatty()):
            out("\n  A partial read needs your confirmation, but there is no terminal to "
                "ask. Re-run with --yes to confirm.")
            return 6
        if ask(f"\n  Read the headers of {len(files)} file(s)? [y/N] ").strip().lower() \
                not in ("y", "yes"):
            out("\n  Aborted. Nothing was read.")
            return 1

    import sources
    client = client or sources.anonymous_client()
    data_dir = data_dir or plan.DATA_DIR
    code = 0
    for r in files:
        bucket_, key = r["bucket"], r["unit"]
        head = client.head_object(Bucket=bucket_, Key=key)
        size, etag = int(head["ContentLength"]), head.get("ETag", "").strip('"')
        if size != r["size"] or (r["etag"] and etag != r["etag"]):
            out(f"\n  SKIPPED {key}: changed on S3 since the catalog was built "
                f"(catalog {r['size']:,} bytes / {r['etag']}; now {size:,} / {etag}). "
                f"Refresh the catalog first.")
            code = 7
            continue
        reader = RangeReader(client, bucket_, key, size, etag, cap)
        report = {"file": f"s3://{bucket_}/{key}", "bucket": bucket_, "key": key,
                  "size": size, "detected_format": "tiff",
                  "extension": key.rsplit(".", 1)[-1].lower(), "warnings": []}
        try:
            report["image"] = images.tiff_header(reader)
            report["status"] = (f"inspected by peek (header only: {reader.bytes_read:,} bytes "
                                f"in {reader.requests} range request(s); no pixels; nothing "
                                f"saved)")
        except images.HeaderError as exc:
            report["status"] = "header not read completely"
            report["warnings"].append(f"{exc}")
            code = code or 7
        report["etag"] = etag
        md = inspect_data.to_markdown(report)
        inspect_data.write_report(report, md, bucket_, key,
                                  out_dir or inspect_data.OUT_DIR)
        log(data_dir, bucket_, key, etag, size, reader)
        out("\n" + md)
        for w in report.get("written", []):
            out(f"  wrote {w}")
    return code


def log(data_dir, bucket, key, etag, size, reader):
    path = data_dir / "_peek_log.csv"
    new = not path.exists()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LOG_FIELDS)
        if new:
            w.writeheader()
        w.writerow({"peeked_at": datetime.datetime.now(datetime.timezone.utc)
                    .isoformat(timespec="seconds"), "bucket": bucket, "key": key,
                    "etag": etag, "size": size, "bytes_read": reader.bytes_read,
                    "requests": reader.requests})
