"""
descriptors.py — fetch only the small files that DESCRIBE a folder unit.

Not a step. Run through index.py:

    python scripts/index.py get --match "1417642733/combined_spatialdata_zarr/" --descriptors --dry-run

A SpatialData Zarr store is hundreds of objects and gigabytes; a Xenium
bundle is ~150 objects and ~15 GB. Their structure is written down in a few
small text files. This module lists the folder (listing needs no
approval), picks those files by name, and hands them to the normal
verified downloader (download.py): same preview, same y/N, same checksum
check, same ledger. The rest of the folder is never touched.

Decisions made here, and why:

1. WHICH FILES DESCRIBE WHICH FOLDER.
   - Zarr store: the consolidated metadata file (zmetadata or .zmetadata)
     holds every group's and array's description in one JSON file. If a
     store has none, every .zgroup / .zarray / .zattrs / zarr.json is
     taken instead (each a few hundred bytes).
   - Xenium bundle: experiment.xenium (run settings: pixel size, panel,
     software versions) and gene_panel.json (the probe panel design).
   Never taken: chunks, parquet parts, images, and summary tables such as
   metrics_summary.csv (its numbers are results, not structure).

2. THE SAME DOWNLOADER, SO THE SAME GUARDS.
   The chosen files become ordinary file rows and go through download.run:
   a HEAD re-check against the listing, ETag verification, .part files,
   the ledger, the y/N question. No second, weaker download path exists.

3. FILES LAND WHERE THEY LIVE ON S3.
   data_cache/<bucket>/<folder>/zmetadata, so the partial store keeps its
   real layout and the inspector can tell which folder it describes.

4. LIMITS.
   At most 10 folder units per call, and the descriptor files of a call
   must total under 50 MB (a real store's metadata is ~100 KB; more would
   mean the rules picked something they should not have).
"""

import plan
import index

DESCRIBABLE = ("zarr", "xenium bundle")
CONSOLIDATED = ("zmetadata", ".zmetadata")
ZARR_META = (".zgroup", ".zarray", ".zattrs", "zarr.json")
XENIUM_META = ("experiment.xenium", "gene_panel.json")
MAX_UNITS = 10
MAX_BYTES = 50_000_000


def pick(fmt, objects):
    """From [(relative key, size, etag)], the descriptor files of one unit."""
    if fmt == "zarr":
        top = [o for o in objects if o[0] in CONSOLIDATED]
        if top:
            return top[:1]
        return [o for o in objects if o[0].rsplit("/", 1)[-1] in ZARR_META]
    if fmt == "xenium bundle":
        return [o for o in objects if o[0] in XENIUM_META]
    return []


def list_unit(client, bucket, prefix):
    """Every object under a folder unit: (relative key, size, ETag). Listing only."""
    out, token = [], None
    while True:
        kw = {"Bucket": bucket, "Prefix": prefix}
        if token:
            kw["ContinuationToken"] = token
        page = client.list_objects_v2(**kw)
        for o in page.get("Contents", []):
            out.append((o["Key"][len(prefix):], int(o["Size"]), o.get("ETag", "").strip('"')))
        if not page.get("IsTruncated"):
            return out
        token = page["NextContinuationToken"]


def expand(client, units, data_dir):
    """Folder units -> file rows (shaped like plan rows) for their descriptors."""
    rows, empty = [], []
    for u in units:
        chosen = pick(u["format"], list_unit(client, u["bucket"], u["unit"]))
        if not chosen:
            empty.append(u["unit"])
        for rel, size, etag in chosen:
            key = u["unit"] + rel
            dest = plan.destination({"bucket": u["bucket"], "unit": key}, data_dir)
            have = dest.stat().st_size if dest.exists() else None
            rows.append({"bucket": u["bucket"], "unit": key, "kind": "file", "size": size,
                         "etag": etag, "format": "json" if rel.endswith(("json", "zmetadata",
                                                                         ".zmetadata", ".xenium"))
                         else "text", "product": "descriptor", "dest": dest,
                         "state": "present" if have == size else
                                  "partial" if have is not None else "missing",
                         "n_objects": 1, "deprecated": False, "missing": "",
                         "size_is_minimum": False, "parent": u["unit"]})
    return rows, empty


def run(match, bucket=None, dry_run=False, yes=False, force=False, client=None,
        catalog_dir=None, data_dir=None, out=print, **kw):
    import download
    import sources
    catalog_dir = catalog_dir or index.CATALOG_DIR
    data_dir = data_dir or plan.DATA_DIR
    rows, _ = plan.select(None, match, None, None, bucket, False, False, catalog_dir)
    units = [r for r in rows if r["kind"] == "collapsed" and r["format"] in DESCRIBABLE]
    other = [r for r in rows if r not in units]
    if other:
        raise SystemExit(f"\nREFUSED: --descriptors works on Zarr stores and Xenium bundles "
                         f"only; {other[0]['unit']} is {other[0]['kind']} / "
                         f"{other[0]['format'] or '?'}. Narrow --match.")
    if not units:
        raise SystemExit("\nNothing matches a Zarr store or Xenium bundle.")
    if len(units) > MAX_UNITS:
        raise SystemExit(f"\nREFUSED: {len(units)} folder units match; --descriptors takes at "
                         f"most {MAX_UNITS}. Narrow --match.")

    client = client or sources.anonymous_client()
    files, empty = expand(client, units, data_dir)
    total = sum(f["size"] for f in files)
    out(f"\n{'=' * 78}\nDESCRIPTOR FILES — the small files that describe each folder\n{'=' * 78}")
    for u in units:
        mine = [f for f in files if f["parent"] == u["unit"]]
        out(f"\n  s3://{u['bucket']}/{u['unit']}")
        out(f"    folder: {u['n_objects']:,} objects, {index.human(u['size'])} — only these "
            f"{len(mine)} file(s), {index.human(sum(f['size'] for f in mine))}, would be fetched:")
        for f in mine[:12]:
            out(f"    {index.human(f['size']):>11s}  {f['unit'][len(u['unit']):]}  [{f['state']}]")
        if len(mine) > 12:
            out(f"    ... {len(mine) - 12} more")
    for e in empty:
        out(f"\n  WARNING: no descriptor file found in {e}")
    out(f"\n  TOTAL: {len(files)} file(s), {index.human(total)} "
        f"(the folders hold {index.human(sum(u['size'] for u in units))})")
    if total > MAX_BYTES:
        out(f"\n  REFUSED: descriptor files total {index.human(total)}, over the "
            f"{index.human(MAX_BYTES)} sanity limit; the selection rules need a look.")
        return 2
    if dry_run:
        out("\n  --dry-run: nothing was transferred.")
        return 0
    return download.run({"units": files}, data_dir, client=client, yes=yes, force=force,
                        out=out, **kw)
