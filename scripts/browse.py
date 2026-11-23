"""
browse.py — answer "what is there for region X / data type Y?" from the catalog.

Not a step. Reads only the saved catalog files under catalog/ — no S3
access, no data files — so every answer is instant and works offline.
Run through index.py:

    python scripts/index.py show                       # overview: every region
    python scripts/index.py show MTG                   # one region, all data types
    python scripts/index.py show MTG --modality rna    # narrow by data type
    python scripts/index.py show MTG --modality rna --files
    python scripts/index.py search 'final-nuclei'      # find units by name
    python scripts/index.py search '*MERFISH*'

Decisions made here, and why:

1. REGION-ONLY AND SHARED DATA ARE SHOWN SEPARATELY, NEVER ADDED TOGETHER.
   Asking for MEC returns MEC-only units, then, in a separate block, units
   that cover MEC *and* other regions (the HIP+MEC neuropathology slides,
   the ten-region subclass objects). Summing them into one MEC total
   would overstate how much MEC data exists and double count against HIP.

2. THE BREAKDOWN IS BY PRODUCT, NOT JUST BY DATA TYPE.
   "MTG snRNA-seq: 1.97 TB" is true but misleading: 1.38 TB of it is
   plot images. Each (bucket, modality) is split into products —
   dataset, per-donor object, metadata table, plots, ... — so the size
   of the thing you actually want is visible.

3. DEPRECATED DATA IS HIDDEN BY DEFAULT, BUT ALWAYS COUNTED.
   previous_objects/ files are superseded releases. They are left out
   unless --deprecated is given, and a line always says how much was
   left out, so nothing disappears silently.

4. REGION NAMES ARE RESOLVED LENIENTLY, DATA-TYPE FILTERS LOOSELY.
   The region goes through regions.resolve(): 'PFC', 'prefrontal' and
   'A9' all mean DFC, and an ambiguous word ('temporal') is refused with
   the candidates. --modality, --product and --bucket are case-
   insensitive substrings ('rna', 'spatial', 'neuropath'), because they
   only narrow a list you can see; nothing is fetched on the strength
   of them.

5. SIZES CARRY THEIR CAVEATS.
   A capped count prints with '≥'. Units the rules could not fully
   classify print with '?' so their labels are not mistaken for facts.
"""

import collections
import fnmatch
import pathlib

import index
import issues as issues_mod
import regions as regions_mod


# ---- loading ------------------------------------------------------------

def load_units(catalog_dir=index.CATALOG_DIR):
    """Every unit row from every bucket's catalog file, sizes as ints."""
    rows = []
    files_dir = catalog_dir / "files"
    if not files_dir.exists():
        raise SystemExit(
            f"\nNo catalog yet in {catalog_dir}. Build it with:\n"
            f"  python scripts/index.py refresh --bucket <bucket>")
    for path in sorted(files_dir.glob("*.csv.gz")):
        for r in index.read_units(path) or []:
            if r["kind"] == "skipped":
                continue
            r["size"] = int(r["size"])
            r["n_objects"] = int(r["n_objects"])
            r["deprecated"] = r["deprecated"] == "True"
            r["size_is_minimum"] = r["size_is_minimum"] == "True"
            r["region_list"] = [c for c in r["regions"].split(";") if c]
            rows.append(r)
    return rows


def _match(value, needle):
    return needle is None or needle.lower() in (value or "").lower()


def filter_units(rows, modality=None, product=None, bucket=None,
                 include_deprecated=False):
    kept, hidden = [], []
    for r in rows:
        if not (_match(r["modality"], modality) and _match(r["product"], product)
                and _match(r["bucket"], bucket)):
            continue
        (kept if include_deprecated or not r["deprecated"] else hidden).append(r)
    return kept, hidden


def _size(rows):
    total = sum(r["size"] for r in rows)
    approx = "≥" if any(r["size_is_minimum"] for r in rows) else ""
    return approx + index.human(total)


def _short_bucket(b):
    return b.replace("sea-ad-", "")


def region_label(codes):
    """'HIP+MEC', or 'HIP+…+V1C (10 regions)' — never silently truncated."""
    if len(codes) <= 3:
        return "+".join(codes)
    return f"{codes[0]}+…+{codes[-1]} ({len(codes)} regions)"


# ---- show ---------------------------------------------------------------

def breakdown(rows):
    """(bucket, modality, product) -> list of rows, largest first."""
    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r["bucket"], r["modality"] or "(none)", r["product"])].append(r)
    return sorted(groups.items(),
                  key=lambda kv: (kv[0][0], kv[0][1], -sum(r["size"] for r in kv[1])))


def show(region=None, modality=None, product=None, bucket=None,
         include_deprecated=False, files=False, limit=30,
         catalog_dir=index.CATALOG_DIR, vocab=None, out=print):
    """Print what the catalog holds. Returns a dict for tests."""
    vocab = vocab or regions_mod.default()
    rows, hidden = filter_units(load_units(catalog_dir), modality, product,
                                bucket, include_deprecated)
    if region is None:
        return _overview(rows, hidden, vocab, out, catalog_dir)

    try:
        code = vocab.resolve(region)
    except regions_mod.RegionError as exc:
        raise SystemExit(f"\n{exc}")
    only = [r for r in rows if r["region_list"] == [code]]
    shared = [r for r in rows if code in r["region_list"] and len(r["region_list"]) > 1]
    hidden_here = [r for r in hidden if code in r["region_list"]]

    title = f"{code} — {vocab.regions[code].name} ({vocab.regions[code].dhba_id})"
    filters = ", ".join(f"{k}~'{v}'" for k, v in
                        (("modality", modality), ("product", product), ("bucket", bucket)) if v)
    out(f"\n{'=' * 78}\n{title}" + (f"   [{filters}]" if filters else "") + f"\n{'=' * 78}")

    for label, block in ((f"{code} only", only),
                         (f"shared: {code} together with other regions", shared)):
        out(f"\n  {label}: {len(block):,} units, {_size(block)}")
        if not block:
            continue
        out(f"    {'bucket':28s} {'modality':34s} {'product':28s} {'units':>6s} {'size':>11s}")
        for (b, m, p), rs in breakdown(block):
            flag = " ?" if any(r["missing"] for r in rs) else ""
            extra = ""
            if label.startswith("shared"):
                combos = sorted({r["regions"] for r in rs})
                extra = "   with " + " | ".join(region_label(c.split(";")) for c in combos[:2]) + \
                        (" ..." if len(combos) > 2 else "")
            out(f"    {_short_bucket(b):28s} {m:34s} {p:28s} {len(rs):6d} {_size(rs):>11s}{flag}{extra}")

    if hidden_here:
        out(f"\n  hidden: {len(hidden_here):,} deprecated units, {_size(hidden_here)} "
            f"(previous_objects/; add --deprecated to include)")
    known = print_issues(only + shared, catalog_dir, out)
    if files:
        _list_files(only + shared, limit, out)
    else:
        out(f"\n  list the files:  python scripts/index.py show {region}"
            + (f" --modality {modality}" if modality else "") + " --files")
    return dict(code=code, only=only, shared=shared, hidden=hidden_here, issues=known)


def _overview(rows, hidden, vocab, out, catalog_dir=index.CATALOG_DIR):
    """Every region (and multi-region group) by modality."""
    by_key = collections.defaultdict(list)
    for r in rows:
        if r["region_list"]:
            key = region_label(r["region_list"])
        elif r["region_from"].startswith("project"):
            key = "(multi-region project)"
        elif "region" in r["missing"].split(";"):
            key = "(unclassified)"
        else:
            key = "(no region)"
        by_key[(key, r["modality"] or "(none)")].append(r)

    def order(k):
        region, _ = k
        return (0 if region in vocab.regions else 1 if "+" in region else 2, region, k[1])

    out(f"\n{'=' * 78}\nSEA-AD CATALOG — current data by region and data type\n{'=' * 78}")
    out(f"  {'region':30s} {'modality':34s} {'units':>6s} {'size':>11s}")
    last = None
    for key in sorted(by_key, key=order):
        region, modality = key
        rs = by_key[key]
        shown = region if region != last else ""
        last = region
        out(f"  {shown[:30]:30s} {modality:34s} {len(rs):6d} {_size(rs):>11s}")
    out(f"\n  total: {len(rows):,} units, {_size(rows)}"
        + (f"   (+ {len(hidden):,} deprecated units, {_size(hidden)}, hidden)" if hidden else ""))
    out("  rows like 'HIP+MEC' are units covering both regions; they are not")
    out("  repeated under HIP or MEC, so the column adds up to the total.")
    n_issues = len(issues_mod.load(pathlib.Path(catalog_dir) / "known_issues.csv"))
    if n_issues:
        out(f"  {n_issues} known upstream issue(s) recorded:  python scripts/index.py issues")
    out("\n  one region in detail:  python scripts/index.py show MTG")
    return dict(groups=by_key, hidden=hidden)


def _list_files(rows, limit, out):
    rows = sorted(rows, key=lambda r: -r["size"])
    out(f"\n  largest {min(limit, len(rows))} of {len(rows):,} units "
        f"(--limit N for more):")
    for r in rows[:limit]:
        size = ("≥" if r["size_is_minimum"] else "") + index.human(r["size"])
        flag = "?" if r["missing"] else " "
        kind = f"  [{r['n_objects']:,} objects]" if r["kind"] == "collapsed" else ""
        out(f"   {flag}{size:>11s}  s3://{r['bucket']}/{r['unit']}{kind}")


# ---- search -------------------------------------------------------------

def search(pattern, modality=None, product=None, bucket=None,
           include_deprecated=False, limit=50, catalog_dir=index.CATALOG_DIR,
           out=print):
    """Find units whose path matches a glob (or contains a word)."""
    pat = pattern if any(ch in pattern for ch in "*?[") else f"*{pattern}*"
    pat = pat.lower()
    rows, hidden = filter_units(load_units(catalog_dir), modality, product,
                                bucket, include_deprecated)
    hits = [r for r in rows if fnmatch.fnmatch(r["unit"].lower(), pat)]
    hidden_hits = [r for r in hidden if fnmatch.fnmatch(r["unit"].lower(), pat)]
    hits.sort(key=lambda r: (r["bucket"], r["unit"]))

    out(f"\nSEARCH {pattern!r}: {len(hits):,} units, {_size(hits)}"
        + (f"   (+ {len(hidden_hits):,} deprecated hidden; --deprecated to include)"
           if hidden_hits else ""))
    for r in hits[:limit]:
        size = ("≥" if r["size_is_minimum"] else "") + index.human(r["size"])
        flag = "?" if r["missing"] else " "
        where = region_label(r["region_list"]) if r["region_list"] else "-"
        out(f"  {flag}{size:>11s}  {where:24s} {r['modality'] or '-':34s} "
            f"{r['product']:26s} s3://{r['bucket']}/{r['unit']}")
    if len(hits) > limit:
        out(f"  ... {len(hits) - limit:,} more (--limit N)")
    if not hits and not hidden_hits:
        out("  no matches. A bare word is wrapped as *word*; matching ignores case.")
    known = print_issues(hits, catalog_dir, out)
    return dict(hits=hits, hidden=hidden_hits, issues=known)


def print_issues(rows, catalog_dir, out):
    """Print the known upstream issues that apply to these rows."""
    known = issues_mod.matching(
        issues_mod.load(pathlib.Path(catalog_dir) / "known_issues.csv"), rows)
    if known:
        out(f"\n  known upstream issues ({len(known)}; details:  python scripts/index.py issues)")
        for line in issues_mod.lines(known):
            out(f"    {line}")
    return known
