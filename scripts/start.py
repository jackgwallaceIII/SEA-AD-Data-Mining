"""
start.py — the front door: check the setup, show what is here, and walk through the workflow.

Not a step. Run through index.py (also what runs when no command is given):

    python scripts/index.py start
    python scripts/index.py start LEC        # the same walk-through, for one region
    python scripts/index.py                  # same as 'start'

Decisions made here, and why:

1. ONE COMMAND A NEWCOMER CAN RUN FIRST, SAFELY.
   It reads only local files (the catalog, the download ledger, the
   reports) and checks which Python packages are importable. It contacts
   nothing and downloads nothing, so it is safe to run before reading any
   documentation.

2. THE SETUP CHECK SAYS HOW TO FIX WHAT IS WRONG.
   Missing required packages are named with what each is needed for, and
   the exact commands to fix them (activate the conda environment, or
   install the pinned versions). The most common mistake, running from
   conda 'base', is detected and named. Exit code 3 when something
   required is missing, so a script can notice.

3. THE WORKFLOW IS FIVE VERBS, EACH WITH ITS REAL COMMAND.
   look -> read -> preview -> download -> describe. The commands shown are
   copied from what works, filled in for the region asked about (or LEC as
   the example), and written without comments so they paste cleanly into
   zsh.

4. THE SAFETY PROMISES ARE STATED UP FRONT.
   Nothing downloads without a preview and a yes; controlled-access data
   is never downloaded; reports show structure, never contents. A
   newcomer should know the limits before they type their first command.
"""

import importlib.util
import os
import pathlib
import sys

import index

ROOT = pathlib.Path(__file__).resolve().parent.parent

REQUIRED = (("boto3", "listing and downloading from the open S3 buckets"),
            ("pandas", "reading CSV tables"),
            ("numpy", "reading matrices"),
            ("h5py", "reading h5ad / HDF5 files"))
OPTIONAL = (("anndata", "writing the test suite's small h5ad files"),
            ("scipy", "the test suite's sparse matrices"),
            ("PIL", "one image test (Pillow)"),
            ("pyarrow", "Parquet row counts — not installed by design; adding it is "
                        "the maintainer's decision"))
MIN_PYTHON = (3, 10)


def environment(find_spec=importlib.util.find_spec, version=sys.version_info, env=None):
    env = os.environ if env is None else env
    missing = [(m, why) for m, why in REQUIRED if find_spec(m) is None]
    optional = [(m, why) for m, why in OPTIONAL if find_spec(m) is None]
    return {"python_ok": tuple(version[:2]) >= MIN_PYTHON,
            "python": ".".join(map(str, version[:3])),
            "conda_env": env.get("CONDA_DEFAULT_ENV", ""),
            "missing": missing, "optional_missing": optional}


def local_state(catalog_dir, data_dir, report_dir):
    import csv
    files = sorted((catalog_dir / "files").glob("*.csv.gz")) if catalog_dir.exists() else []
    st = {"buckets": [p.name[:-len(".csv.gz")] for p in files], "units": 0, "size": 0,
          "listed": {}, "downloads": 0, "download_size": 0, "reports": 0, "issues": 0}
    if files:
        import browse
        rows, _ = browse.filter_units(browse.load_units(catalog_dir), None, None, None, False)
        st["units"], st["size"] = len(rows), sum(r["size"] for r in rows)
        log = catalog_dir / "refresh_log.csv"
        if log.exists():
            with open(log, newline="") as fh:
                for r in csv.DictReader(fh):
                    if r.get("mode", "refresh") in ("refresh", ""):
                        st["listed"][r["bucket"]] = r["finished"][:10]
    import download
    latest = {str(r["file"]): int(r["size"] or 0) for r in download.read_ledger(data_dir)}
    present = {p: s for p, s in latest.items() if pathlib.Path(p).is_file()}
    st["downloads"], st["download_size"] = len(present), sum(present.values())
    if report_dir.exists():
        st["reports"] = len([p for p in report_dir.glob("*.md") if p.name != "README.md"])
    issues_file = catalog_dir / "known_issues.csv"
    if issues_file.exists():
        import issues
        st["issues"] = len(issues.load(issues_file))
    return st


def _example(code, catalog_dir):
    """A small, real file of the region to use in the example commands."""
    try:
        import browse
        rows, _ = browse.filter_units(browse.load_units(catalog_dir), None, None, None, False)
    except SystemExit:
        return None
    # A first download should be small and typical: the smallest RNA metadata
    # table or per-donor object (tens of MB), else any of those, else any file.
    mine = [r for r in rows if r["region_list"] == [code] and r["kind"] == "file"]
    kinds = ("cell metadata table", "per-donor object")
    pool = ([r for r in mine if r["product"] in kinds and r["modality"] == "snRNA-seq"]
            or [r for r in mine if r["product"] in kinds] or mine)
    return min(pool, key=lambda r: r["size"]) if pool else None


def run(region=None, catalog_dir=None, data_dir=None, report_dir=None, env_check=None,
        vocab=None, out=print):
    import regions as regions_mod
    catalog_dir = pathlib.Path(catalog_dir or index.CATALOG_DIR)
    data_dir = pathlib.Path(data_dir or ROOT / "data_cache")
    report_dir = pathlib.Path(report_dir or ROOT / "results" / "reports")
    vocab = vocab or regions_mod.default()
    code = "LEC"
    if region:
        try:
            code = vocab.resolve(region)
        except regions_mod.RegionError as exc:
            out(f"\n{exc}\n\nRegions: " + ", ".join(f"{c} ({r.name})"
                                                     for c, r in vocab.regions.items()))
            return 2
    env = env_check or environment()
    bar = "=" * 78
    out(f"\n{bar}\nSEA-AD DATA TOOL — start here\n{bar}")
    out("  Finds data from the Seattle Alzheimer's Disease Brain Cell Atlas across its\n"
        "  open S3 buckets and the SEA-AD website, previews downloads, fetches what you\n"
        "  approve, and writes structure reports (shapes, columns, types — never contents).")

    # 1. setup
    out(f"\n1. SETUP CHECK")
    ok = env["python_ok"] and not env["missing"]
    out(f"   {'ok ' if env['python_ok'] else 'NO '} Python {env['python']}"
        + ("" if env["python_ok"] else f" (needs {MIN_PYTHON[0]}.{MIN_PYTHON[1]} or newer)"))
    for m, why in REQUIRED:
        bad = any(m == x for x, _ in env["missing"])
        out(f"   {'NO ' if bad else 'ok '} {m:8s} {why}")
    for m, why in env["optional_missing"]:
        out(f"   --  {m:8s} optional, not installed: {why}")
    if env["conda_env"]:
        out(f"   conda environment: {env['conda_env']}")
    if not ok:
        out("\n   Something required is missing. In the project folder, run:")
        out("     conda activate seaad")
        out("   If that environment does not exist yet, create it and install the pinned "
            "versions:")
        out("     conda create -n seaad python=3.12 -y")
        out("     conda activate seaad")
        out("     pip install -r requirements-lock.txt")
        if env["conda_env"] == "base":
            out("   (You are in conda 'base', the usual cause: the packages live in 'seaad'.)")

    # 2. what is here
    st = local_state(catalog_dir, data_dir, report_dir)
    out(f"\n2. WHAT IS HERE (local files only; nothing was contacted)")
    if st["buckets"]:
        out(f"   catalog:   {st['units']:,} current units, {index.human(st['size'])}, in "
            f"{len(st['buckets'])} bucket(s); listed "
            + ", ".join(f"{b.replace('sea-ad-', '')} {d}" for b, d in sorted(st["listed"].items())))
    else:
        out("   catalog:   none yet. Build it by listing the buckets (no downloads):")
        out("     python scripts/index.py refresh --bucket sea-ad-single-cell-profiling")
    out(f"   downloads: {st['downloads']} file(s), {index.human(st['download_size'])} "
        f"in data_cache/")
    out(f"   reports:   {st['reports']} region page(s) in results/reports/; "
        f"{st['issues']} known upstream issue(s)")

    # 3. workflow
    ex = _example(code, catalog_dir) if st["buckets"] else None
    reg = vocab.regions[code]
    out(f"\n3. THE WORKFLOW — shown for {code} ({reg.name})"
        + ("" if region else "; add a region, e.g.  start MTG"))
    steps = [
        ("LOOK at what exists (catalog only, instant)",
         ["python scripts/index.py show", f"python scripts/index.py show {code}",
          "python scripts/index.py search final-nuclei", "python scripts/index.py sources"]),
        ("READ the region's structure report",
         [f"python scripts/index.py report {code}", f"open results/reports/{code}.md"]),
        ("PREVIEW a download (nothing is transferred)",
         [f"python scripts/index.py get --match \"{ex['unit']}\" --dry-run" if ex else
          f"python scripts/index.py get {code} --product \"cell metadata table\" --dry-run",
          "python scripts/index.py get --external web-donor-metadata --dry-run"]),
        ("DOWNLOAD it: the same command without --dry-run; it asks y/N first",
         [f"python scripts/index.py get --match \"{ex['unit']}\"" if ex else
          f"python scripts/index.py get {code} --product \"cell metadata table\""]),
        ("DESCRIBE what you downloaded (structure only)",
         ["python scripts/index.py inspect --downloaded", "python scripts/index.py report --all"]),
    ]
    for i, (title, cmds) in enumerate(steps, 1):
        out(f"\n   {i}. {title}")
        out("      conda activate seaad")
        for c in cmds:
            out(f"      {c}")
    if ex:
        out(f"\n   The example file is {index.human(ex['size'])} ({ex['product']}).")
    out("\n   Large folders and images without downloading them:")
    out("      python scripts/index.py get --descriptors --match \"<zarr or xenium folder>/\" --dry-run")
    out("      python scripts/index.py peek --match \"<image file>\" --dry-run")

    # 4. promises
    out(f"\n4. WHAT THIS TOOL WILL AND WILL NOT DO")
    for line in ("Nothing downloads without a size preview and your yes (y/N or --yes).",
                 "Controlled-access data (Synapse / AD Knowledge Portal) is linked, never "
                 "downloaded.",
                 "Reports show structure only: shapes, column names, types, counts. Donor IDs, "
                 "barcodes, free text and few-donor values are withheld.",
                 "Every download is checked against the server's checksum and recorded in "
                 "data_cache/_ledger.csv."):
        out(f"   - {line}")
    out(f"\n5. MORE: python scripts/index.py --help   (every command: "
        f"python scripts/index.py <command> --help)")
    return 0 if ok else 3
