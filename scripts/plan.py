"""
plan.py — the download preview: exactly what would be fetched, before anything is.

Not a step. Builds a download plan from the saved catalog and prints it.
It contains NO download code; nothing in this module can transfer data.
Run through index.py:

    python scripts/index.py get MTG --modality rna --product metadata --dry-run
    python scripts/index.py get --match '*SEAAD_MTG_MERFISH*' --dry-run
    python scripts/index.py get LEC --product dataset --dry-run

Decisions made here, and why:

1. SELECTION REUSES `show`, SO WHAT YOU SAW IS WHAT YOU GET.
   A plan is built with the same region / --modality / --product /
   --bucket / --deprecated filters as `show`, plus --match for a path
   glob. Region selection includes units SHARED with other regions
   (e.g. the ten-region subclass objects) only if --include-shared is
   given, because those are large and are not "MTG data" in the narrow
   sense; the plan says how many were left out.

2. A BROAD SELECTION IS REFUSED.
   With no region and no --match, the selection would be the whole
   9 TB catalog. That is refused outright, like `step0b_fetch.py get
   <directory>` with no narrowing.

3. EVERY UNIT IS SHOWN WITH ITS SIZE, ITS DESTINATION AND ITS STATE.
   Destination is data_cache/<bucket>/<key>: the bucket's own layout,
   so a file's local path tells you exactly where it came from. State
   is 'missing', 'partial' (on disk but the wrong size — a truncated
   download), or 'present' (right size). Present units are excluded
   from the transfer total.

4. --max-gb IS CHECKED HERE, BEFORE A DOWNLOAD COULD EVER START.
   Default 50 GB, the same as step0b_fetch.py. Over the cap, the plan
   is printed and then refused, with the exact --max-gb needed.

5. WARNINGS ARE ABOUT INTERPRETATION, NOT JUST SIZE.
   The plan flags: collapsed units with very many objects (slow, and
   rarely what you want: the plots folder is 1.1 million images);
   duplicates, in two strengths (see duplicate_groups): byte-identical
   (same size and ETag — e.g. scCODA abundance objects repeated across
   analyses) and probable (same name and size but multipart ETags that
   cannot be compared — e.g. the two MERFISH h5ads); deprecated data;
   units whose labels are
   incomplete ('?'); and pickle files, which can run code when opened
   and which the inspect step will never open.

6. THE CATALOG'S AGE IS PRINTED.
   Sizes come from the last refresh, not from S3 right now. The plan
   says when that was. The download step (step 8) re-checks every file
   against S3 immediately before transferring it.

7. KNOWN UPSTREAM ISSUES ARE SHOWN, NOT ENFORCED.
   Any row of catalog/known_issues.csv that applies to a selected unit is
   printed with its severity. Even a 'problem' does not block the
   download: you may want the broken file precisely to look at it. The
   warning makes sure the choice is an informed one.
"""

import collections
import pathlib

import browse
import index
import issues as issues_mod
import regions as regions_mod

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data_cache"

MANY_OBJECTS = 1000


def destination(row, data_dir=DATA_DIR):
    return data_dir / row["bucket"] / row["unit"]


def local_state(row, data_dir=DATA_DIR):
    """('missing' | 'partial' | 'present', bytes on disk)."""
    p = destination(row, data_dir)
    if row["kind"] == "collapsed":
        if not p.is_dir():
            return "missing", 0
        have = sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
    else:
        if not p.is_file():
            return "missing", 0
        have = p.stat().st_size
    return ("present" if have == row["size"] else "partial"), have


def select(region=None, match=None, modality=None, product=None, bucket=None,
           include_deprecated=False, include_shared=False,
           catalog_dir=index.CATALOG_DIR, vocab=None):
    """Units to plan, plus what was left out and why."""
    if region is None and match is None:
        raise SystemExit(
            "\nRefusing to plan a download of the whole catalog (~9 TB).\n"
            "Give a region (e.g. MTG) and/or --match '<glob>', and narrow with\n"
            "--modality / --product. See what is there first:\n"
            "  python scripts/index.py show MTG")
    vocab = vocab or regions_mod.default()
    rows, hidden = browse.filter_units(browse.load_units(catalog_dir), modality,
                                       product, bucket, include_deprecated)
    left_out = collections.Counter()
    if region is not None:
        try:
            code = vocab.resolve(region)
        except regions_mod.RegionError as exc:
            raise SystemExit(f"\n{exc}")
        here = [r for r in rows if code in r["region_list"]]
        shared = [r for r in here if len(r["region_list"]) > 1]
        rows = here if include_shared else [r for r in here if len(r["region_list"]) == 1]
        if shared and not include_shared:
            left_out["shared with other regions (--include-shared)"] = len(shared)
        hidden = [r for r in hidden if code in r["region_list"]]
    if match is not None:
        import fnmatch
        pat = (match if any(ch in match for ch in "*?[") else f"*{match}*").lower()
        rows = [r for r in rows if fnmatch.fnmatch(r["unit"].lower(), pat)]
        hidden = [r for r in hidden if fnmatch.fnmatch(r["unit"].lower(), pat)]
    if hidden:
        left_out["deprecated (--deprecated)"] = len(hidden)
    return rows, left_out


def warnings_for(rows):
    out = []
    many = [r for r in rows if r["kind"] == "collapsed" and r["n_objects"] > MANY_OBJECTS]
    for r in many:
        out.append(f"{r['unit']} is {r['n_objects']:,} separate objects "
                   f"({r['product']}); downloading it is slow and rarely needed")
    out.extend(duplicate_warnings(rows))
    dep = [r for r in rows if r["deprecated"]]
    if dep:
        out.append(f"{len(dep)} unit(s) are deprecated (previous_objects/): "
                   f"superseded by a newer release")
    if any(r["missing"] for r in rows):
        out.append("units marked '?' have incomplete labels "
                   "(python scripts/index.py review --bucket <bucket>)")
    if any(r["format"] == "pickle" for r in rows):
        out.append("pickle files can execute code when opened; inspect will "
                   "report their size only and never open them")
    return out


def duplicate_groups(rows):
    """Two kinds of duplicate, kept apart because they mean different things.

    identical: same size AND same ETag. For a single-part upload the
    ETag is the file's MD5, so this is byte-identical content whatever
    the names are. (Same multipart ETag + size also means same content.)
    likely: same file NAME and size, but both ETags are multipart with
    DIFFERENT part counts, so the ETag cannot confirm or refute it.
    Same name, same size, different single-part MD5 = NOT a duplicate.
    """
    files = [r for r in rows if r["kind"] == "file" and r["etag"]]
    by_etag = collections.defaultdict(list)
    for r in files:
        by_etag[(r["size"], r["etag"])].append(r)
    identical = [rs for rs in by_etag.values() if len(rs) > 1]
    seen = {id(r) for rs in identical for r in rs}

    by_name = collections.defaultdict(list)
    for r in files:
        if id(r) not in seen and "-" in r["etag"]:
            by_name[(r["unit"].rsplit("/", 1)[-1], r["size"])].append(r)
    likely = [rs for rs in by_name.values()
              if len(rs) > 1 and len({r["etag"] for r in rs}) > 1]
    return identical, likely


def duplicate_warnings(rows, show=3):
    identical, likely = duplicate_groups(rows)
    out = []
    for label, groups, why in (
            ("identical files", identical,
             "same size and same ETag: byte-identical content"),
            ("probable duplicates", likely,
             "same name and size; multipart uploads, so the ETag cannot confirm")):
        if not groups:
            continue
        extra = sum(len(g) - 1 for g in groups)
        lines = [f"{len(groups)} group(s) of {label} ({why}); "
                 f"{extra} unit(s) you likely do not need" +
                 (f" — first {show} shown:" if len(groups) > show else ":")]
        for g in groups[:show]:
            lines.extend(f"          {r['unit']}" for r in g[:4])
            if len(g) > 4:
                lines.append(f"          ... and {len(g) - 4} more with the same content")
            lines.append("          --")
        out.append("\n".join(lines[:-1]))
    return out


def catalog_dates(rows, catalog_dir=index.CATALOG_DIR):
    import csv
    path = catalog_dir / "refresh_log.csv"
    if not path.exists():
        return {}
    last = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("mode", "refresh") in ("refresh", ""):
                last[r["bucket"]] = r["finished"]
    return {b: last.get(b, "unknown") for b in {r["bucket"] for r in rows}}


def build(rows, max_gb=50.0, data_dir=DATA_DIR):
    units = []
    for r in sorted(rows, key=lambda r: (r["bucket"], r["unit"])):
        state, have = local_state(r, data_dir)
        units.append({**r, "state": state, "have": have,
                      "dest": destination(r, data_dir)})
    need = sum(u["size"] for u in units if u["state"] != "present")
    return {"units": units,
            "selected": sum(u["size"] for u in units),
            "present": sum(u["size"] for u in units if u["state"] == "present"),
            "need": need,
            "over_cap": need > max_gb * 1e9,
            "minimum": any(u["size_is_minimum"] for u in units)}


def render(plan, left_out, warns, dates, max_gb, data_dir=DATA_DIR, out=print):
    units = plan["units"]
    h = index.human
    out(f"\n{'=' * 78}\nDOWNLOAD PLAN — preview only, nothing is transferred\n{'=' * 78}")
    if not units:
        out("  nothing matched." + (" Left out: " + ", ".join(
            f"{n} {why}" for why, n in left_out.items()) if left_out else ""))
        return 1
    out(f"  {'size':>11s}  {'state':8s} unit")
    for u in units:
        flag = "?" if u["missing"] else " "
        size = ("≥" if u["size_is_minimum"] else "") + h(u["size"])
        objs = f"  [{u['n_objects']:,} objects]" if u["kind"] == "collapsed" else ""
        out(f" {flag}{size:>11s}  {u['state']:8s} s3://{u['bucket']}/{u['unit']}{objs}")
    approx = "at least " if plan["minimum"] else ""
    n_present = sum(1 for u in units if u["state"] == "present")
    n_partial = sum(1 for u in units if u["state"] == "partial")
    out(f"\n  selected:      {len(units):,} units, {approx}{h(plan['selected'])}")
    out(f"  already here:  {n_present:,} units, {h(plan['present'])}")
    if n_partial:
        out(f"  partial:       {n_partial:,} units on disk at the wrong size "
            f"(truncated?) — will be re-fetched")
    out(f"  TO DOWNLOAD:   {len(units) - n_present:,} units, {approx}{h(plan['need'])}")
    out(f"  destination:   {data_dir}/<bucket>/<key>")
    for why, n in left_out.items():
        out(f"  left out:      {n:,} units {why}")
    for b, when in sorted(dates.items()):
        out(f"  sizes as of:   {when} (last refresh of {b})")
    for w in warns:
        out(f"  WARNING: {w}")
    if plan["over_cap"]:
        out(f"\n  REFUSED: {h(plan['need'])} exceeds --max-gb {max_gb:g}. Narrow the "
            f"selection, or raise the cap deliberately: --max-gb {plan['need'] / 1e9:.0f}")
        return 2
    return 0


def preview(region=None, match=None, modality=None, product=None, bucket=None,
            include_deprecated=False, include_shared=False, max_gb=50.0,
            catalog_dir=None, data_dir=None, out=print):
    # Resolved at call time, so tests can point both at temporary folders.
    catalog_dir = catalog_dir or index.CATALOG_DIR
    data_dir = data_dir or DATA_DIR
    rows, left_out = select(region, match, modality, product, bucket,
                            include_deprecated, include_shared, catalog_dir)
    plan = build(rows, max_gb, data_dir)
    known = issues_mod.matching(issues_mod.load(pathlib.Path(catalog_dir) / "known_issues.csv"),
                                rows)
    plan["issues"] = known
    warns = warnings_for(rows) + [
        f"known issue {i['id']} ({i['severity']}): {i['summary']} — {n:,} selected "
        f"unit(s); details: python scripts/index.py issues" for i, n in known]
    code = render(plan, left_out, warns,
                  catalog_dates(rows, catalog_dir), max_gb, data_dir, out)
    return plan, code
