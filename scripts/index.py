"""
index.py — build and refresh the SEA-AD catalog from the open S3 buckets.

Not a pipeline step. Read-only against S3: it LISTS a bucket, classifies
what it finds, and writes three files under catalog/. Nothing is
downloaded and no data file is opened.

    python scripts/index.py refresh --bucket sea-ad-single-cell-profiling
    python scripts/index.py refresh --bucket sea-ad-single-cell-profiling --dry-run
    python scripts/index.py review  --bucket sea-ad-single-cell-profiling
    python scripts/index.py reclassify --bucket sea-ad-single-cell-profiling
    python scripts/index.py show MTG --modality rna      (see browse.py)
    python scripts/index.py search final-nuclei

Files written (decision 1):

    catalog/files/<bucket>.csv.gz   one row per UNIT (a file, a collapsed
                                    folder, or a skipped folder)
    catalog/catalog.csv             one row per (bucket, region, modality)
    catalog/refresh_log.csv         one row per refresh, appended

Decisions made here, and why:

1. TWO LEVELS: UNITS AND SUMMARY ROWS.
   The summary is what a person browses ("MTG, snRNA-seq: 44 units,
   1.2 TB"). The unit table is what a download preview needs ("which 3 of
   those 44, and how big is each"). Every summary number is a sum over
   unit rows, and a test checks the two always agree.

2. A UNIT COVERING SEVERAL REGIONS GETS ONE ROW, NOT ONE PER REGION.
   The neuropathology MTG+STG slides form one row with region 'MTG;STG',
   and the multiregion files one row listing all ten. Copying them into
   each region's row would make region totals add up to more than the
   bucket holds, which is the kind of quiet double count that makes
   later numbers wrong. Multiregion_2026 model outputs with no stated
   region form a 'multi-region (project)' row.

3. DEPRECATED UNITS ARE COUNTED SEPARATELY, NOT IN THE MAIN TOTALS.
   previous_objects/ data is kept in the unit table (flagged) and shown
   in the summary's deprecated_* columns. total_size is current data
   only, which is the number that answers "how much would I download".

4. SKIPPED FOLDERS ARE RECORDED WITH UNKNOWN SIZE.
   They appear in the unit table with kind 'skipped' and an empty size,
   and are not in any summary total. Empty means "not listed", not zero.

5. A CAPPED COUNT MARKS EVERY TOTAL IT TOUCHES.
   If --max-keys-per-unit stopped a count early, that unit's row and its
   summary row carry size_is_minimum=True. The default is no cap.

6. REFRESH IS WHOLE-BUCKET AND MANUAL, AND REPORTS WHAT CHANGED.
   A partial walk cannot say what was removed, so refresh always walks
   the whole bucket. It compares against the previous unit table and
   prints units added, removed and changed (ETag, size or object count
   differs). Nothing refreshes automatically: data changing silently
   under an analysis is what the ABC release pin exists to prevent, and
   these buckets have no pin, so the diff is the only warning you get.

7. OTHER BUCKETS' ROWS ARE LEFT ALONE.
   Refreshing one bucket rewrites only that bucket's summary rows.

8. ACCESS IS RECORDED PER ROW.
   Every bucket here is 'open'. Controlled-access sources (Synapse IDs
   you will supply) are added in a later step and can never be fetched.

9. RULE CHANGES DO NOT NEED A RE-LISTING.
   `reclassify` re-applies the current rules to the saved unit table
   with no S3 access: seconds instead of 36 minutes for the
   neuropathology bucket. It cannot see added or removed files, so it
   is not a substitute for `refresh`, and refresh_log.csv records which
   of the two produced each catalog ('mode' column).
"""

import argparse
import csv
import datetime
import gzip
import io
import pathlib
import sys
import time

import pathrules
import regions as regions_mod
import sources

ROOT = pathlib.Path(__file__).resolve().parent.parent
CATALOG_DIR = ROOT / "catalog"

UNIT_FIELDS = [
    "bucket", "unit", "kind", "n_objects", "size", "size_is_minimum",
    "etag", "last_modified", "regions", "region_from", "modality",
    "modality_from", "product", "format", "per_donor", "deprecated",
    "skip", "missing",
]
SUMMARY_FIELDS = [
    "bucket", "region", "region_name", "region_aliases", "modality",
    "access", "prefixes", "formats", "products", "n_units", "n_objects",
    "total_size", "size_is_minimum", "deprecated_units", "deprecated_size",
    "needs_review", "indexed_at",
]
LOG_FIELDS = [
    "bucket", "started", "finished", "seconds", "units", "current_size",
    "deprecated_size", "skipped_folders", "size_is_minimum", "added",
    "removed", "changed", "mode",
]


def human(n):
    for unit, div in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= div:
            return f"{n / div:,.2f} {unit}"
    return f"{n} B"


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


# ========================= units ========================================

def unit_rows(bucket, records, vocab):
    """Turn walker records into unit-table rows (dicts of strings)."""
    rows = []
    for rec in records:
        if isinstance(rec, sources.FileObject):
            c = pathrules.classify(bucket, rec.key, vocab)
            row = dict(kind="file", n_objects=1, size=rec.size,
                       size_is_minimum=False, etag=rec.etag,
                       last_modified=rec.last_modified)
        elif isinstance(rec, sources.CollapsedUnit):
            c = pathrules.classify_unit(bucket, rec.prefix, rec.format,
                                        rec.product, vocab)
            row = dict(kind="collapsed", n_objects=rec.n_objects,
                       size=rec.total_size,
                       size_is_minimum=not rec.complete, etag="",
                       last_modified=rec.last_modified)
        else:
            c = None
            row = dict(kind="skipped", n_objects="", size="",
                       size_is_minimum="", etag="", last_modified="")
        if c is None:
            row.update(bucket=bucket, unit=rec.prefix, regions="",
                       region_from="", modality="", modality_from="",
                       product="skipped", format="", per_donor=False,
                       deprecated=False, skip=rec.reason, missing="")
        elif c.skip:
            continue        # a skipped single key (e.g. checkpoint copy)
        else:
            row.update(bucket=bucket, unit=c.unit if row["kind"] == "file" else rec.prefix,
                       regions=";".join(c.regions), region_from=c.region_from,
                       modality=c.modality, modality_from=c.modality_from,
                       product=c.product, format=c.format,
                       per_donor=c.per_donor, deprecated=c.deprecated,
                       skip="", missing=";".join(c.missing))
        rows.append(row)
    return rows


def region_key(row):
    if row["regions"]:
        return row["regions"]
    if row["region_from"].startswith("project"):
        return "multi-region (project)"
    if "region" in row["missing"].split(";"):
        return "unclassified"
    return "(none)"


# ========================= summary ======================================

def summarise(bucket, rows, vocab, indexed_at):
    """One summary row per (region, modality) for this bucket."""
    groups = {}
    for r in rows:
        if r["kind"] == "skipped":
            continue
        key = (region_key(r), r["modality"] or "(none)")
        groups.setdefault(key, []).append(r)

    out = []
    for (region, modality), rs in sorted(groups.items()):
        current = [r for r in rs if not _true(r["deprecated"])]
        old = [r for r in rs if _true(r["deprecated"])]
        codes = [c for c in region.split(";") if c in vocab.regions]
        names = "; ".join(vocab.regions[c].name for c in codes)
        aliases = "; ".join(
            sorted({a for c in codes
                    for a in (c, vocab.regions[c].dhba_id,
                              *vocab.regions[c].aliases)}))
        out.append({
            "bucket": bucket, "region": region, "region_name": names,
            "region_aliases": aliases, "modality": modality,
            "access": "open",
            "prefixes": ";".join(sorted({r["unit"].split("/")[0] + "/"
                                         for r in rs if "/" in r["unit"]}))
                        or "(bucket root)",
            "formats": ";".join(sorted({r["format"] for r in current})),
            "products": ";".join(sorted({r["product"] for r in current})),
            "n_units": len(current),
            "n_objects": sum(int(r["n_objects"]) for r in current),
            "total_size": sum(int(r["size"]) for r in current),
            "size_is_minimum": any(_true(r["size_is_minimum"]) for r in rs),
            "deprecated_units": len(old),
            "deprecated_size": sum(int(r["size"]) for r in old),
            "needs_review": sum(1 for r in current if r["missing"]),
            "indexed_at": indexed_at,
        })
    return out


def _true(v):
    return v is True or str(v) == "True"


# ========================= files ========================================

def units_path(bucket, catalog_dir):
    return catalog_dir / "files" / f"{bucket}.csv.gz"


def read_units(path):
    if not path.exists():
        return None
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def write_units(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    # mtime=0 keeps the gzip bytes identical for identical content, so an
    # unchanged refresh produces no git diff.
    with open(tmp, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=UNIT_FIELDS)
                w.writeheader()
                for r in sorted(rows, key=lambda r: r["unit"]):
                    w.writerow({k: r[k] for k in UNIT_FIELDS})
    tmp.replace(path)


def write_summary(path, bucket, new_rows):
    """Replace this bucket's rows in catalog.csv, keep everyone else's."""
    kept = []
    if path.exists():
        with open(path, newline="") as fh:
            kept = [r for r in csv.DictReader(fh) if r["bucket"] != bucket]
    rows = kept + new_rows
    rows.sort(key=lambda r: (r["bucket"], r["region"], r["modality"]))
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def append_log(path, entry):
    rows = []
    if path.exists():
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh))
    # Older logs predate the 'mode' column; those rows were all refreshes.
    for r in rows:
        r.setdefault("mode", "refresh")
        if r["mode"] in ("", None):
            r["mode"] = "refresh"
    rows.append({**{f: "" for f in LOG_FIELDS}, **entry})
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LOG_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


# ========================= diff =========================================

def diff(old_rows, new_rows):
    """(added, removed, changed) unit names between two unit tables."""
    if old_rows is None:
        return None
    old = {r["unit"]: r for r in old_rows}
    new = {r["unit"]: r for r in new_rows}
    added = sorted(set(new) - set(old))
    removed = sorted(set(old) - set(new))
    changed = sorted(u for u in set(old) & set(new)
                     if any(str(old[u][f]) != str(new[u][f])
                            for f in ("etag", "size", "n_objects")))
    return added, removed, changed


# ========================= commands =====================================

def refresh(bucket, client=None, catalog_dir=CATALOG_DIR, max_keys_per_unit=None,
            dry_run=False, vocab=None, out=print, progress=None):
    """Walk one bucket, write the catalog files. Returns a result dict."""
    if bucket not in pathrules.BUCKETS:
        raise SystemExit(
            f"\nUnknown bucket {bucket!r}. Known: "
            f"{', '.join(sorted(pathrules.BUCKETS))}")
    vocab = vocab or regions_mod.default()
    src = sources.S3Source(bucket, client=client)
    started, t0 = now(), time.time()
    rows = unit_rows(bucket, src.walk("", max_keys_per_unit, progress), vocab)
    finished, seconds = now(), round(time.time() - t0)

    summary = summarise(bucket, rows, vocab, finished)
    upath = units_path(bucket, catalog_dir)
    d = diff(read_units(upath), rows)

    current = sum(s["total_size"] for s in summary)
    deprecated = sum(s["deprecated_size"] for s in summary)
    skipped = [r for r in rows if r["kind"] == "skipped"]
    minimum = any(s["size_is_minimum"] for s in summary)
    result = dict(rows=rows, summary=summary, diff=d, seconds=seconds,
                  current=current, deprecated=deprecated, skipped=skipped,
                  minimum=minimum, folder_markers=src.folder_markers)

    if not dry_run:
        write_units(upath, rows)
        write_summary(catalog_dir / "catalog.csv", bucket, summary)
        append_log(catalog_dir / "refresh_log.csv", {
            "bucket": bucket, "started": started, "finished": finished,
            "seconds": seconds, "units": len(rows) - len(skipped),
            "current_size": current, "deprecated_size": deprecated,
            "skipped_folders": len(skipped), "size_is_minimum": minimum,
            "added": "" if d is None else len(d[0]),
            "removed": "" if d is None else len(d[1]),
            "changed": "" if d is None else len(d[2]),
            "mode": "refresh",
        })
    return result


def reclassify(bucket, catalog_dir=CATALOG_DIR, dry_run=False, vocab=None):
    """Re-apply the CURRENT rules to the saved unit table. No S3 access.

    Decision 9: listing is slow (the neuropathology bucket took 36 min)
    and rules change often; the stored keys, sizes and ETags do not.
    This re-derives every label from the stored unit names, so a rule
    change costs seconds. It cannot discover new or removed files —
    that still needs refresh — and the log records which one ran.
    """
    if bucket not in pathrules.BUCKETS:
        raise SystemExit(f"\nUnknown bucket {bucket!r}.")
    vocab = vocab or regions_mod.default()
    upath = units_path(bucket, catalog_dir)
    old = read_units(upath)
    if old is None:
        raise SystemExit(f"\nNo catalog for {bucket} yet. Run refresh first.")

    rows = []
    for r in old:
        if r["kind"] == "skipped":
            rows.append(r)
            continue
        if r["kind"] == "file":
            c = pathrules.classify(bucket, r["unit"], vocab)
            if c.skip:
                continue            # a newly skippable key (e.g. Thumbs.db)
        else:
            c = pathrules.classify_unit(
                bucket, r["unit"], r["format"],
                pathrules.COLLAPSE_PRODUCT.get(r["format"], r["product"]), vocab)
        rows.append({**r,
                     "n_objects": int(r["n_objects"]), "size": int(r["size"]),
                     "size_is_minimum": _true(r["size_is_minimum"]),
                     "regions": ";".join(c.regions), "region_from": c.region_from,
                     "modality": c.modality, "modality_from": c.modality_from,
                     "product": c.product, "format": c.format,
                     "per_donor": c.per_donor, "deprecated": c.deprecated,
                     "skip": "", "missing": ";".join(c.missing)})

    stamp = now()
    summary = summarise(bucket, rows, vocab, stamp)
    before = {r["unit"]: r for r in old}
    changed = sum(1 for b in rows
                  if any(str(before[b["unit"]][f]) != str(b[f]) for f in
                         ("regions", "modality", "product", "format", "missing")))
    dropped = len(old) - len(rows)
    result = dict(rows=rows, summary=summary, relabelled=changed,
                  dropped=dropped, diff=None, seconds=0,
                  current=sum(s["total_size"] for s in summary),
                  deprecated=sum(s["deprecated_size"] for s in summary),
                  skipped=[r for r in rows if r["kind"] == "skipped"],
                  minimum=any(s["size_is_minimum"] for s in summary),
                  folder_markers=0)
    if not dry_run:
        write_units(upath, rows)
        write_summary(catalog_dir / "catalog.csv", bucket, summary)
        append_log(catalog_dir / "refresh_log.csv", {
            "bucket": bucket, "started": stamp, "finished": stamp,
            "seconds": 0, "units": len(rows) - len(result["skipped"]),
            "current_size": result["current"],
            "deprecated_size": result["deprecated"],
            "skipped_folders": len(result["skipped"]),
            "size_is_minimum": result["minimum"], "mode": "reclassify"})
    return result


def report(bucket, res, dry_run, catalog_dir, out=print):
    approx = "at least " if res["minimum"] else ""
    out(f"\n{'=' * 74}\nCATALOG — s3://{bucket}  ({res['seconds']} s)\n{'=' * 74}")
    out(f"  {'region':28s} {'modality':34s} {'units':>6s} {'current':>12s} {'deprecated':>12s}  review")
    for s in res["summary"]:
        region = s["region"] if len(s["region"]) <= 28 else s["region"][:25] + "..."
        out(f"  {region:28s} {s['modality']:34s} {s['n_units']:6d} "
            f"{('≥' if s['size_is_minimum'] else '') + human(s['total_size']):>12s} "
            f"{human(s['deprecated_size']) if s['deprecated_units'] else '-':>12s}  "
            f"{s['needs_review'] or ''}")
    n_units = sum(s["n_units"] for s in res["summary"])
    n_dep = sum(s["deprecated_units"] for s in res["summary"])
    out(f"\n  current:    {n_units:,} units, {approx}{human(res['current'])}")
    out(f"  deprecated: {n_dep:,} units, {human(res['deprecated'])}  (previous_objects/)")
    out(f"  skipped:    {len(res['skipped'])} folders, size not listed"
        + (": " + ", ".join(f"{r['unit']} ({r['skip']})" for r in res["skipped"])
           if res["skipped"] else ""))
    out(f"  folder-marker keys dropped: {res['folder_markers']}")
    review = sum(s["needs_review"] for s in res["summary"])
    out(f"  needing review: {review} current units  "
        f"(python scripts/index.py review --bucket {bucket})")

    d = res["diff"]
    if d is None:
        out("\n  First refresh of this bucket: nothing to compare against.")
    else:
        out(f"\n  Since the last refresh: {len(d[0])} added, {len(d[1])} removed, "
            f"{len(d[2])} changed")
        for label, items in (("added", d[0]), ("removed", d[1]), ("changed", d[2])):
            for u in items[:10]:
                out(f"    {label:8s} {u}")
            if len(items) > 10:
                out(f"    ... and {len(items) - 10} more {label}")

    if dry_run:
        out("\n  --dry-run: nothing was written.")
    else:
        upath = units_path(bucket, catalog_dir)
        out(f"\n  wrote {upath.relative_to(ROOT) if upath.is_relative_to(ROOT) else upath}"
            f"  ({human(upath.stat().st_size)})")
        out(f"  wrote {(catalog_dir / 'catalog.csv').name}, appended refresh_log.csv")


def cmd_refresh(args):
    def progress(prefix, n_folders, n_files):
        progress.n += 1
        if progress.n % 25 == 0:
            print(f"  ... {progress.n} folders listed (at {prefix or '/'})",
                  file=sys.stderr, flush=True)
    progress.n = 0
    res = refresh(args.bucket, max_keys_per_unit=args.max_keys_per_unit,
                  dry_run=args.dry_run, progress=progress)
    report(args.bucket, res, args.dry_run, CATALOG_DIR)
    return 0


def cmd_reclassify(args):
    res = reclassify(args.bucket, dry_run=args.dry_run)
    print(f"\n  re-applied current rules to the saved catalog (no S3 access)")
    print(f"  units relabelled: {res['relabelled']:,}   "
          f"now skipped (dropped): {res['dropped']:,}")
    report(args.bucket, res, args.dry_run, CATALOG_DIR)
    return 0


def cmd_review(args):
    rows = read_units(units_path(args.bucket, CATALOG_DIR))
    if rows is None:
        raise SystemExit(f"\nNo catalog for {args.bucket} yet. Run refresh first.")
    todo = [r for r in rows if r["missing"]
            and (args.all or not _true(r["deprecated"]))]
    print(f"\n{len(todo)} units with a field the rules could not determine"
          + ("" if args.all else " (current only; --all includes deprecated)"))
    for r in todo:
        size = human(int(r["size"])) if r["size"] else "?"
        print(f"  missing {r['missing']:24s} {size:>11s}  {r['unit']}")
    return 0


def cmd_show(args):
    import browse
    browse.show(args.region, modality=args.modality, product=args.product,
                bucket=args.bucket, include_deprecated=args.deprecated,
                files=args.files, limit=args.limit)
    return 0


def cmd_search(args):
    import browse
    browse.search(args.pattern, modality=args.modality, product=args.product,
                  bucket=args.bucket, include_deprecated=args.deprecated,
                  limit=args.limit)
    return 0


def cmd_get(args):
    import download
    import plan
    if args.external:
        import external
        if args.region or args.match:
            raise SystemExit("\n--external selects by source id only (see: index.py sources).")
        return external.run_get(args.external, max_gb=args.max_gb, dry_run=args.dry_run,
                                yes=args.yes, force=args.force)
    if args.descriptors:
        import descriptors
        if args.region:
            raise SystemExit("\n--descriptors selects folders with --match only.")
        return descriptors.run(args.match, bucket=args.bucket, dry_run=args.dry_run,
                               yes=args.yes, force=args.force, confirm_gb=args.confirm_gb)
    p, code = plan.preview(args.region, match=args.match, modality=args.modality,
                           product=args.product, bucket=args.bucket,
                           include_deprecated=args.deprecated,
                           include_shared=args.include_shared, max_gb=args.max_gb)
    if code:
        return code
    if args.dry_run:
        print("\n  --dry-run: nothing was transferred.")
        return 0
    return download.run(p, plan.DATA_DIR, yes=args.yes, force=args.force,
                        confirm_gb=args.confirm_gb)


def cmd_start(args):
    import start
    return start.run(args.region)


def cmd_sources(args):
    import external
    if args.check:
        return external.check()
    external.report()
    return 0


def cmd_report(args):
    import report
    report.run(args.region, all_regions=args.all)
    return 0


def cmd_peek(args):
    import peek
    return peek.run(args.match, bucket=args.bucket, max_kb=args.max_kb,
                    dry_run=args.dry_run, yes=args.yes)


def cmd_issues(args):
    import issues
    issues.report(include_resolved=args.all)
    return 0


def cmd_inspect(args):
    import inspect_data
    paths = list(args.paths)
    if args.downloaded:
        paths += inspect_data.downloaded_files()
    if not paths and args.downloaded:
        raise SystemExit("\nNothing to inspect: data_cache/_ledger.csv lists no downloaded file "
                         "that is still on disk.\nPreview and download one first; "
                         "python scripts/index.py start shows how.")
    if not paths:
        raise SystemExit("\nGive a file path, or --downloaded for every file in "
                         "data_cache/_ledger.csv.")
    failed = []
    for p in paths:
        try:
            report, md = inspect_data.inspect(p)
        except Exception as exc:          # one bad file must not stop the rest
            failed.append(p)
            print(f"\n  FAILED to inspect {p}: {type(exc).__name__}: {exc}")
            continue
        print("\n" + md)
        for w in report.get("written", []):
            print(f"  wrote {w}")
    if failed:
        print(f"\n  {len(failed)} of {len(paths)} file(s) could not be inspected (see FAILED "
              f"above); the others were.")
        return 1
    return 0


def _filter_flags(p, limit):
    p.add_argument("--modality", metavar="TEXT",
                   help="keep data types containing TEXT (e.g. rna, atac, spatial, merfish)")
    p.add_argument("--product", metavar="TEXT",
                   help="keep products containing TEXT (e.g. dataset, donor, metadata)")
    p.add_argument("--bucket", metavar="TEXT",
                   help="keep buckets containing TEXT (e.g. single, spatial, neuropath)")
    p.add_argument("--deprecated", action="store_true",
                   help="include superseded previous_objects/ data")
    p.add_argument("--limit", type=int, default=limit, metavar="N",
                   help=f"rows to list (default {limit})")


def build_parser():
    ap = argparse.ArgumentParser(
        prog="index.py",
        description="Find, preview, fetch and describe SEA-AD data. New here? Run: "
                    "python scripts/index.py start")
    sub = ap.add_subparsers(dest="command")

    p = sub.add_parser("start", help="START HERE: setup check, what is here, the workflow")
    p.add_argument("region", nargs="?", help="show the walk-through for this region, e.g. MTG")
    p.set_defaults(func=cmd_start)

    p = sub.add_parser("refresh", help="walk one whole bucket and rewrite its catalog rows")
    p.add_argument("--bucket", required=True, choices=sorted(pathrules.BUCKETS))
    p.add_argument("--max-keys-per-unit", type=int, default=None, metavar="N",
                   help="stop counting a collapsed unit after N keys; sizes then "
                        "become minimums (default: exact)")
    p.add_argument("--dry-run", action="store_true",
                   help="walk and report, but write nothing")
    p.set_defaults(func=cmd_refresh)

    p = sub.add_parser("reclassify",
                       help="re-apply current rules to the saved catalog (no S3 access)")
    p.add_argument("--bucket", required=True, choices=sorted(pathrules.BUCKETS))
    p.add_argument("--dry-run", action="store_true", help="report, but write nothing")
    p.set_defaults(func=cmd_reclassify)

    p = sub.add_parser("show", help="what the catalog holds: overview, or one region in detail")
    p.add_argument("region", nargs="?",
                   help="region code or name (MTG, PFC, hippocampus ...); omit for the overview")
    _filter_flags(p, 30)
    p.add_argument("--files", action="store_true", help="list the units, largest first")
    p.set_defaults(func=cmd_show)

    p = sub.add_parser("search", help="find units by path (glob, or a bare word)")
    p.add_argument("pattern")
    _filter_flags(p, 50)
    p.set_defaults(func=cmd_search)

    p = sub.add_parser("get", help="download plan for a selection (preview; see plan.py)")
    p.add_argument("region", nargs="?", help="region code or name (MTG, PFC ...)")
    p.add_argument("--match", metavar="GLOB", help="unit path glob, or a bare word")
    _filter_flags(p, 0)
    p.add_argument("--include-shared", action="store_true",
                   help="also include units shared with other regions")
    p.add_argument("--max-gb", type=float, default=50.0, metavar="N",
                   help="refuse plans larger than this (default: 50)")
    p.add_argument("--dry-run", action="store_true",
                   help="print the plan and stop; nothing is transferred")
    p.add_argument("--yes", "-y", action="store_true",
                   help="do not ask for confirmation (for scripts)")
    p.add_argument("--confirm-gb", type=float, default=0.0, metavar="N",
                   help="ask before downloads larger than this (default: 0 = always ask)")
    p.add_argument("--force", action="store_true",
                   help="re-download files already present")
    p.add_argument("--external", metavar="ID",
                   help="an open web file from catalog/external_sources.csv (glob allowed); "
                        "controlled-access sources are refused")
    p.add_argument("--descriptors", action="store_true",
                   help="for Zarr stores / Xenium bundles: fetch only the small files that "
                        "describe them (metadata), never the data")
    p.set_defaults(func=cmd_get)

    p = sub.add_parser("inspect", help="structure report for downloaded files (no contents)")
    p.add_argument("paths", nargs="*", help="file(s) under data_cache/")
    p.add_argument("--downloaded", action="store_true",
                   help="inspect every file recorded in data_cache/_ledger.csv")
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser("sources", help="sources outside the buckets: open web files, "
                                        "controlled-access studies (links only)")
    p.add_argument("--check", action="store_true",
                   help="HEAD each open file (headers only, no content): size and date")
    p.set_defaults(func=cmd_sources)

    p = sub.add_parser("report", help="one data-structure report per region (results/reports/)")
    p.add_argument("region", nargs="?", help="region name or code, e.g. LEC")
    p.add_argument("--all", action="store_true", help="every region, plus an index page")
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("peek", help="read image headers from S3 without downloading (TIFF/SVS)")
    p.add_argument("--match", required=True, help="glob or word selecting 1-20 image files")
    p.add_argument("--bucket", choices=sorted(pathrules.BUCKETS))
    p.add_argument("--max-kb", type=int, default=4096,
                   help="most bytes read per file, in KB (default 4096)")
    p.add_argument("--dry-run", action="store_true", help="show what would be read; read nothing")
    p.add_argument("--yes", "-y", action="store_true", help="confirm without asking")
    p.set_defaults(func=cmd_peek)

    p = sub.add_parser("issues", help="known problems in the upstream data (catalog/known_issues.csv)")
    p.add_argument("--all", action="store_true", help="include resolved issues")
    p.set_defaults(func=cmd_issues)

    p = sub.add_parser("review", help="units the rules could not fully classify")
    p.add_argument("--bucket", required=True, choices=sorted(pathrules.BUCKETS))
    p.add_argument("--all", action="store_true", help="include deprecated units")
    p.set_defaults(func=cmd_review)
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command is None:             # no command: the front door
        import start
        return start.run()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
