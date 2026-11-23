"""
sources.py — read-only, anonymous listing of the open SEA-AD S3 buckets.

Not a step. The catalog (step 4) uses S3Source.walk() to turn a bucket
into units it can classify. This module only LISTS. It has no download
or read method at all; those come later, behind the size preview.

    python scripts/sources.py ls   sea-ad-single-cell-profiling MTG/
    python scripts/sources.py walk sea-ad-single-cell-profiling MTG/ATACseq/
    python scripts/sources.py walk sea-ad-quantitative-neuropathology \\
        middle-temporal-gyrus/H19.33.004/ --max-keys-per-unit 2000

Decisions made here, and why:

1. ANONYMOUS ACCESS, ALWAYS.
   The client is built with UNSIGNED signatures, the boto3 equivalent of
   `aws s3 ls --no-sign-request`. Credentials on this machine (if any)
   are never sent, so a listing can never act as, or be billed to, your
   AWS account, and it is guaranteed to see only what the public sees.

2. LISTING ONLY.
   The class exposes list operations and nothing else. The test suite
   uses botocore's Stubber, which fails on any S3 call it was not told
   to expect, so an accidental GetObject would fail the tests.

3. THE WALK GOES ONE FOLDER LEVEL AT A TIME.
   S3 has no folders; Delimiter='/' makes it group keys by the next '/'
   so they can be treated as folders. Walking level by level means a
   folder can be skipped (web-viewer assets) or collapsed (a Zarr store,
   1.1 million plot images) as soon as its name is seen, using
   pathrules.folder_rule(). A flat listing would have to page through
   every key first. In the single-cell bucket that is ~1.9 million keys,
   ~1.13 million of them plot images.

4. COLLAPSED UNITS ARE SIZED, AND A CAP MAKES THE SIZE A LOWER BOUND.
   S3 cannot report a folder's size; the only way is to list every key
   in it and add them up. The plots folder took ~4 minutes to count.
   --max-keys-per-unit stops counting early. The unit is then marked
   complete=False and its size is a minimum ("at least"), never shown as
   exact. Default is no cap: correct and slow beats fast and wrong.

5. SKIPPED FOLDERS ARE REPORTED, NOT LISTED.
   A skipped folder yields one SkippedFolder record with the reason.
   Its contents are never listed, so its size is unknown, and the record
   says so rather than reporting zero.

6. FOLDER-MARKER KEYS ARE DROPPED.
   Some buckets contain zero-byte keys equal to a folder name
   ('MTG/RNAseq/'). They are console artefacts, not files, and listing a
   folder returns its own marker. They are counted, not yielded.

7. ETag AND LastModified ARE KEPT ON EVERY FILE.
   These buckets have no manifest and no release numbers. The ETag and
   the modification time are the only way to tell later whether a file
   changed after it was catalogued (plan section 3). The surrounding
   quotes S3 puts on ETags are stripped.
"""

import argparse
import dataclasses
import sys

import pathrules

DEFAULT_REGION = "us-west-2"


@dataclasses.dataclass(frozen=True)
class FileObject:
    key: str
    size: int
    etag: str
    last_modified: str     # ISO 8601, UTC


@dataclasses.dataclass(frozen=True)
class CollapsedUnit:
    prefix: str            # ends with '/'
    n_objects: int
    total_size: int
    last_modified: str     # newest object seen
    complete: bool         # False -> n_objects and total_size are minimums
    format: str
    product: str


@dataclasses.dataclass(frozen=True)
class SkippedFolder:
    prefix: str
    reason: str            # contents were not listed; size unknown


def anonymous_client(region=DEFAULT_REGION):
    """An S3 client that never signs requests (no credentials sent)."""
    try:
        import boto3
        from botocore import UNSIGNED
        from botocore.config import Config
    except ImportError as exc:
        raise SystemExit(
            f"\n{exc}. boto3/botocore are in requirements-lock.txt but not "
            f"in this Python ({sys.prefix}).\n"
            f"Activate the project environment first:  conda activate seaad")

    return boto3.client(
        "s3", region_name=region,
        config=Config(signature_version=UNSIGNED,
                      retries={"max_attempts": 10, "mode": "standard"}))


def _iso(dt):
    return dt.isoformat() if hasattr(dt, "isoformat") else str(dt)


class S3Source:
    """Read-only lister for one bucket."""

    def __init__(self, bucket, client=None, region=DEFAULT_REGION):
        self.bucket = bucket
        self.client = client or anonymous_client(region)
        self.folder_markers = 0

    # ---- one level ------------------------------------------------------

    def list_level(self, prefix=""):
        """Immediate sub-folders and files under prefix (one level).

        Returns (folders, files): folder prefixes ending '/', and
        FileObject records. Pages through every result.
        """
        folders, files = [], []
        token = None
        while True:
            kw = {"Bucket": self.bucket, "Prefix": prefix, "Delimiter": "/"}
            if token:
                kw["ContinuationToken"] = token
            page = self.client.list_objects_v2(**kw)
            folders.extend(p["Prefix"] for p in page.get("CommonPrefixes", []))
            for o in page.get("Contents", []):
                if o["Key"].endswith("/"):
                    self.folder_markers += 1
                    continue
                files.append(FileObject(o["Key"], o["Size"],
                                        o.get("ETag", "").strip('"'),
                                        _iso(o["LastModified"])))
            if not page.get("IsTruncated"):
                return folders, files
            token = page["NextContinuationToken"]

    # ---- a whole folder, flat (for collapsed units) ---------------------

    def summarise(self, prefix, rule, max_keys=None):
        """Count and size every key under prefix, up to max_keys."""
        n = size = 0
        newest = ""
        token = None
        while True:
            kw = {"Bucket": self.bucket, "Prefix": prefix}
            if token:
                kw["ContinuationToken"] = token
            page = self.client.list_objects_v2(**kw)
            for o in page.get("Contents", []):
                if o["Key"].endswith("/"):
                    continue
                n += 1
                size += o["Size"]
                newest = max(newest, _iso(o["LastModified"]))
                if max_keys and n >= max_keys:
                    return CollapsedUnit(prefix, n, size, newest, False,
                                         rule.format, rule.product)
            if not page.get("IsTruncated"):
                return CollapsedUnit(prefix, n, size, newest, True,
                                     rule.format, rule.product)
            token = page["NextContinuationToken"]

    # ---- recursive walk -------------------------------------------------

    def walk(self, prefix="", max_keys_per_unit=None, progress=None):
        """Yield FileObject, CollapsedUnit and SkippedFolder records.

        Every key under prefix is accounted for exactly once: as a file,
        inside one collapsed unit, or inside one skipped folder.
        """
        if prefix:
            rule = pathrules.folder_rule(self.bucket, prefix)
            if rule.action == "skip":
                yield SkippedFolder(prefix, rule.reason)
                return
            if rule.action == "collapse":
                yield self.summarise(prefix, rule, max_keys_per_unit)
                return
        folders, files = self.list_level(prefix)
        if progress:
            progress(prefix, len(folders), len(files))

        # A Hive-partitioned table shows up as a folder whose children
        # are 'key=value' folders; that needs the children's names. The
        # check comes BEFORE yielding this level's files, because the
        # flat summary below counts them too.
        names = [f[len(prefix):].rstrip("/") for f in folders]
        if prefix and any("=" in n for n in names):
            rule = pathrules.folder_rule(self.bucket, prefix, names)
            yield self.summarise(prefix, rule, max_keys_per_unit)
            return
        yield from files
        for f in folders:
            yield from self.walk(f, max_keys_per_unit, progress)


# ========================= cli ==========================================

def _human(n):
    for unit, div in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= div:
            return f"{n / div:,.2f} {unit}"
    return f"{n} B"


def cmd_ls(args):
    src = S3Source(args.bucket)
    folders, files = src.list_level(args.prefix)
    for f in folders:
        rule = pathrules.folder_rule(args.bucket, f)
        note = "" if rule.action == "descend" else f"   [{rule.action}: {rule.reason or rule.product}]"
        print(f"  {'folder':>11s}  {f}{note}")
    for o in files:
        print(f"  {_human(o.size):>11s}  {o.key}")
    print(f"\n  {len(folders)} folders, {len(files)} files "
          f"({_human(sum(o.size for o in files))}) directly under "
          f"s3://{args.bucket}/{args.prefix}")
    return 0


def cmd_walk(args):
    src = S3Source(args.bucket)
    n_files = n_units = n_skip = 0
    size = 0
    lower_bound = False
    for rec in src.walk(args.prefix, args.max_keys_per_unit):
        if isinstance(rec, FileObject):
            n_files += 1
            size += rec.size
            if args.show_files:
                print(f"  {_human(rec.size):>11s}  {rec.key}")
        elif isinstance(rec, CollapsedUnit):
            n_units += 1
            size += rec.total_size
            lower_bound |= not rec.complete
            approx = "" if rec.complete else "at least "
            print(f"  {'UNIT':>11s}  {rec.prefix}  [{rec.product}: "
                  f"{approx}{rec.n_objects:,} objects, {approx}{_human(rec.total_size)}]")
        else:
            n_skip += 1
            print(f"  {'SKIPPED':>11s}  {rec.prefix}  [{rec.reason}; not listed]")
    approx = "at least " if lower_bound else ""
    print(f"\n  s3://{args.bucket}/{args.prefix}")
    print(f"  files: {n_files:,}   collapsed units: {n_units:,}   "
          f"skipped folders: {n_skip:,}   folder markers dropped: "
          f"{src.folder_markers:,}")
    print(f"  total listed size: {approx}{_human(size)}"
          + ("   (a capped unit makes this a minimum)" if lower_bound else ""))
    return 0


def build_parser():
    ap = argparse.ArgumentParser(
        prog="sources.py",
        description="Read-only, anonymous listing of open SEA-AD S3 buckets. "
                    "Nothing is downloaded.")
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("ls", help="one level: sub-folders and files")
    p.add_argument("bucket")
    p.add_argument("prefix", nargs="?", default="")
    p.set_defaults(func=cmd_ls)

    p = sub.add_parser("walk", help="recursive walk with skip/collapse rules")
    p.add_argument("bucket")
    p.add_argument("prefix", nargs="?", default="")
    p.add_argument("--max-keys-per-unit", type=int, default=None, metavar="N",
                   help="stop counting a collapsed unit after N keys; its "
                        "size is then reported as a minimum (default: no cap)")
    p.add_argument("--show-files", action="store_true",
                   help="print every file, not just units and totals")
    p.set_defaults(func=cmd_walk)
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "walk" and not args.prefix:
        raise SystemExit(
            "\nRefusing to walk a whole bucket from the command line; the "
            "single-cell bucket alone is ~1.9 million keys.\n"
            "Give a prefix (e.g. MTG/ATACseq/). Whole-bucket walks belong "
            "to the catalog refresh in step 4.")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
