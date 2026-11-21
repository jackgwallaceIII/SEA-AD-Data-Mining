"""
step0b_fetch.py — browse and fetch ABC Atlas data from the Allen S3 bucket.

This is the only script in the project that downloads data. Everything
downstream reads from the local abc_cache/ directory. abc_cache/ is
gitignored (~17 GB for this project), and it does not need to be shared:
every file is pinned by a versioned manifest with an md5 per file, so this
script rebuilds a byte-identical cache from a clean checkout. That is the
reason the manifest-pinned multiregion release was chosen over the 2024
flagship MTG atlas — see decision 2 in step0_discover.py.

    # Reproduce exactly what this project uses (the MTG pipeline):
    python scripts/step0b_fetch.py preset mtg

    # Browse what exists
    python scripts/step0b_fetch.py list
    python scripts/step0b_fetch.py list SEA-AD-Multiregion-10X
    python scripts/step0b_fetch.py regions
    python scripts/step0b_fetch.py search PFC
    python scripts/step0b_fetch.py search 'MTG*'

    # Fetch whatever you want
    python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X MTG-10X/raw
    python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X --kind metadata
    python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X 'HIP-10X/*'
    python scripts/step0b_fetch.py get Allen-CCF-2020 --kind image_volumes --dry-run

Run `--help` on any command for its flags.

Decisions made here, and why:

1. THE RELEASE IS PINNED, NOT RESOLVED AS "LATEST".
   AbcProjectCache.from_s3_cache() loads whichever manifest was last used
   in the cache directory, falling back to the newest on S3. Both are
   mutable state, so the data underneath an analysis can change without
   anything in the code changing. MANIFEST below is therefore hard-coded.
   --manifest overrides it for exploration, and says loudly that it has.
   To move the project to a newer release, change the constant
   deliberately and re-run every step; do not let it drift.

2. SELECTION IS EXPLICIT AND SIZE IS ALWAYS SHOWN FIRST.
   The whole catalog is 0.92 TB, and a single directory's expression
   matrices can be 144 GB — this project uses exactly one 15.19 GB file
   out of that. Every fetch therefore prints the file list with sizes and
   a total before transferring anything, refuses totals above --max-gb,
   and asks for confirmation above --confirm-gb. --dry-run prints the
   plan and exits. Nothing is ever fetched by directory sweep implicitly.

3. PER-FILE SIZES AND HASHES COME FROM THE MANIFEST, NOT THE API.
   The cache's helpers only report directory-level totals, so the
   manifest JSON is read directly. It also gives each file's version and
   md5, which is what makes "which data did I use" answerable later.

4. ONE KIND-AGNOSTIC ACCESSOR IS USED FOR EVERY FETCH.
   get_data_path() and get_metadata_path() are both deprecated in
   abc_atlas_access and forward to get_file_path(), which dispatches on
   directory + file_name alone. No file name maps to more than one
   extension anywhere in the manifest, so that pair is a unique key, and
   the same call works for matrices, metadata, image volumes and
   mapmycells alike.

5. HASHES ARE VERIFIED, NEVER SKIPPED.
   get_file_path checks each file's md5 against the manifest by default,
   and skip_hash_check is deliberately not exposed. A 15 GB download
   truncated by a dropped connection is the likeliest failure here, and a
   silently truncated h5ad produces wrong QC metrics rather than an
   error. Re-running is cheap: files already present at the right size
   are reported as cached and excluded from the transfer total.
"""

import argparse
import fnmatch
import json
import pathlib
import sys
import warnings

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

from seaad import REGION_ALIASES, REGION_INFO, region_description
from abc_atlas_access.abc_atlas_cache.cloud_cache import OutdatedManifestWarning

# Staying on an older release is the whole point of the pin (decision 1),
# so the library's multi-line warning about it is noise here. open_cache()
# reports the same fact in one line. Nothing else is silenced — in
# particular LargeDataSizeWarning is left on.
warnings.filterwarnings("ignore", category=OutdatedManifestWarning)

# Anchored to the project root via this file's own location, so the script
# runs from any working directory.
ROOT = pathlib.Path(__file__).resolve().parent.parent
CACHE_DIR = ROOT / "abc_cache"

# ---- The pin. See decision 1. Do not let this drift. --------------------
# step0_discover.py and step3_annotation.py carry the same value; this
# script is the authority on it.
MANIFEST = "releases/20261101/manifest.json"

KINDS = ("expression_matrices", "metadata", "image_volumes", "mapmycells")

# The directory cmd_regions reports on.
DIRECTORY_10X = "SEA-AD-Multiregion-10X"

# Named, reproducible file sets. A preset is the unit of "what this
# analysis used": it names every file explicitly so a collaborator gets
# the same bytes, rather than whatever a directory happens to contain.
PRESETS = {
    "mtg": {
        "description": "MTG pipeline (steps 1-3): one 15.19 GB matrix + 8 metadata tables",
        "files": [
            ("SEA-AD-Multiregion-10X", "MTG-10X/raw"),
            ("SEA-AD-Multiregion-10X", "cell_metadata"),
            ("SEA-AD-Multiregion-10X", "donor"),
            ("SEA-AD-Multiregion-10X", "disease"),
            ("SEA-AD-Multiregion-10X", "library"),
            ("SEA-AD-Multiregion-taxonomy", "cluster_annotation_term"),
            ("SEA-AD-Multiregion-taxonomy", "cluster_annotation_term_set"),
            ("SEA-AD-Multiregion-taxonomy", "cell_to_cluster_membership"),
            ("SEA-AD-Multiregion-taxonomy", "cluster_to_cluster_annotation_membership"),
        ],
    },
    "mtg-metadata": {
        "description": "MTG annotation tables only — enough to inspect, skips the 15.19 GB matrix",
        "files": [
            ("SEA-AD-Multiregion-10X", "cell_metadata"),
            ("SEA-AD-Multiregion-10X", "donor"),
            ("SEA-AD-Multiregion-10X", "disease"),
            ("SEA-AD-Multiregion-10X", "library"),
            ("SEA-AD-Multiregion-taxonomy", "cluster_annotation_term"),
            ("SEA-AD-Multiregion-taxonomy", "cluster_annotation_term_set"),
            ("SEA-AD-Multiregion-taxonomy", "cell_to_cluster_membership"),
            ("SEA-AD-Multiregion-taxonomy", "cluster_to_cluster_annotation_membership"),
        ],
    },
}


# ========================= helpers ======================================

def section(title):
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def human(n):
    """Bytes as a short human-readable string."""
    for unit, div in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= div:
            return f"{n / div:,.2f} {unit}"
    return f"{n} B"


def open_cache(manifest_name):
    """Connect to the S3-backed cache and pin the release."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = AbcProjectCache.from_s3_cache(CACHE_DIR)

    if manifest_name == "latest":
        manifest_name = cache.latest_manifest_file
        print(f"  --manifest latest resolved to {manifest_name}")

    try:
        cache.load_manifest(manifest_name)
    except Exception as exc:
        raise SystemExit(
            f"\nCould not load manifest {manifest_name!r}: {exc}\n"
            f"Available: {cache.list_manifest_file_names}"
        )

    if cache.current_manifest != manifest_name:
        raise SystemExit(
            f"\nPin failed: asked for {manifest_name}, "
            f"got {cache.current_manifest}."
        )

    if manifest_name != MANIFEST:
        print(f"\n  *** NOT the project pin. Using {manifest_name}; the")
        print(f"  *** project is pinned to {MANIFEST}. Results computed")
        print(f"  *** against a different release are not comparable.")
    elif cache.latest_manifest_file != MANIFEST:
        print(f"  (a newer release exists: {cache.latest_manifest_file} — "
              f"staying on the pin)")

    path = CACHE_DIR / manifest_name
    if not path.exists():
        raise SystemExit(f"Manifest not on disk after load: {path}")
    return cache, json.loads(path.read_text()), manifest_name


def walk(node, parts):
    """Yield (file_name, files_dict) for every leaf under a kind subtree.

    Nesting depth differs by kind — expression matrices are
    region/normalization while metadata is a flat name — so descend until
    a node carries a 'files' mapping and join the levels above it with
    '/', which is exactly the file_name the cache API expects.
    """
    if isinstance(node, dict) and isinstance(node.get("files"), dict):
        yield "/".join(parts), node["files"]
        return
    if isinstance(node, dict):
        for key, child in node.items():
            yield from walk(child, parts + [key])


def catalog(manifest, directory=None, kinds=KINDS):
    """Flatten the manifest into a list of file records."""
    out = []
    for d, by_kind in manifest["file_listing"].items():
        if directory and d != directory:
            continue
        for kind, tree in by_kind.items():
            if kind not in kinds:
                continue
            for name, files in walk(tree, []):
                for ext, entry in files.items():
                    out.append({
                        "directory": d, "kind": kind, "file_name": name,
                        "ext": ext, "size": entry["size"],
                        "version": entry["version"],
                        "md5": entry["file_hash"],
                        "relative_path": entry["relative_path"],
                        "url": entry["url"],
                    })
    return out


def local_state(rec):
    """Is this file already on disk at the manifest's size?"""
    p = CACHE_DIR / rec["relative_path"]
    if not p.exists():
        return "missing", 0
    have = p.stat().st_size
    if have == rec["size"]:
        return "cached", have
    return "partial", have


# ========================= commands =====================================

def cmd_list(args):
    _, manifest, mname = open_cache(args.manifest)

    if not args.directory:
        section(f"DIRECTORIES — {mname}")
        recs = catalog(manifest)
        by_dir = {}
        for r in recs:
            by_dir.setdefault(r["directory"], []).append(r)
        print(f"  {'directory':38s} {'files':>6s} {'size':>12s}  kinds")
        print(f"  {'-' * 38} {'-' * 6:>6s} {'-' * 12:>12s}  {'-' * 24}")
        for d in sorted(by_dir):
            rs = by_dir[d]
            kinds = sorted({r["kind"] for r in rs})
            short = ",".join(k.replace("expression_matrices", "expr")
                              .replace("image_volumes", "img")
                              .replace("mapmycells", "mmc") for k in kinds)
            print(f"  {d:38s} {len(rs):6d} {human(sum(r['size'] for r in rs)):>12s}  {short}")
        total = sum(r["size"] for r in recs)
        print(f"\n  {len(by_dir)} directories, {len(recs)} files, {human(total)} total")
        print(f"\n  Next: {pathlib.Path(sys.argv[0]).name} list <directory>")
        return 0

    recs = catalog(manifest, directory=args.directory,
                   kinds=(args.kind,) if args.kind else KINDS)
    if not recs:
        known = sorted(manifest["file_listing"])
        raise SystemExit(
            f"\nNo files for directory {args.directory!r}"
            + (f" and kind {args.kind!r}" if args.kind else "")
            + f".\nKnown directories:\n  " + "\n  ".join(known)
        )

    section(f"{args.directory} — {mname}")
    print(f"  {'kind':20s} {'file_name':44s} {'ext':7s} {'size':>11s}  state")
    print(f"  {'-' * 20} {'-' * 44} {'-' * 7} {'-' * 11:>11s}  {'-' * 7}")
    for r in sorted(recs, key=lambda r: (r["kind"], r["file_name"])):
        state, _ = local_state(r)
        print(f"  {r['kind']:20s} {r['file_name']:44s} {r['ext']:7s} "
              f"{human(r['size']):>11s}  {state}")
    print(f"\n  {len(recs)} files, {human(sum(r['size'] for r in recs))} total")
    return 0


def cmd_search(args):
    _, manifest, mname = open_cache(args.manifest)
    pat = args.pattern if any(c in args.pattern for c in "*?[") else f"*{args.pattern}*"
    pat = pat.lower()

    # Also match the region's anatomical name and any alias, so
    # `search PFC` and `search prefrontal` both find DFC-10X.
    bare = args.pattern.lower().strip("*")
    alias_hit = REGION_ALIASES.get(bare)
    anatomy = {code for code, (desc, _) in REGION_INFO.items()
               if bare and bare in desc.lower()}
    if alias_hit:
        anatomy.add(alias_hit)

    def matches(r):
        full = f"{r['directory']}/{r['file_name']}".lower()
        if fnmatch.fnmatch(full, pat) or fnmatch.fnmatch(r["file_name"].lower(), pat):
            return True
        region = r["file_name"].split("/")[0]
        return region in anatomy

    recs = [r for r in catalog(manifest, kinds=(args.kind,) if args.kind else KINDS)
            if matches(r)]

    section(f"SEARCH {args.pattern!r} — {mname}")
    if not recs:
        print("  no matches. Patterns are globbed and case-insensitive;")
        print("  a bare word is wrapped as *word*. Try `list` to browse.")
        return 1
    print(f"  {'directory':30s} {'file_name':40s} {'size':>11s}  state")
    print(f"  {'-' * 30} {'-' * 40} {'-' * 11:>11s}  {'-' * 7}")
    for r in sorted(recs, key=lambda r: (r["directory"], r["file_name"])):
        state, _ = local_state(r)
        desc = region_description(r["file_name"].split("/")[0])
        print(f"  {r['directory']:30s} {r['file_name']:40s} "
              f"{human(r['size']):>11s}  {state}"
              + (f"   {desc}" if desc else ""))
    print(f"\n  {len(recs)} matches, {human(sum(r['size'] for r in recs))} total")
    return 0


def cmd_regions(args):
    """The region table: code, anatomical name, size, what is cached.

    The codes are opaque — nobody guesses that the prefrontal cortex is
    filed as DFC — so this is the command that answers "what regions are
    there and which one do I want".
    """
    _, manifest, mname = open_cache(args.manifest)
    em = manifest["file_listing"][DIRECTORY_10X]["expression_matrices"]

    section(f"REGIONS IN {DIRECTORY_10X} — {mname}")
    print(f"  {'code':9s} {'anatomical name':44s} {'raw':>9s} {'log2':>9s}  state")
    print(f"  {'-' * 9} {'-' * 44} {'-' * 9:>9s} {'-' * 9:>9s}  {'-' * 7}")
    for r in sorted(em):
        raw = em[r]["raw"]["files"]["h5ad"]
        state, _ = local_state({"relative_path": raw["relative_path"],
                                "size": raw["size"]})
        print(f"  {r:9s} {region_description(r)[:44]:44s} "
              f"{human(raw['size']):>9s} "
              f"{human(em[r]['log2']['files']['h5ad']['size']):>9s}  {state}")

    total = sum(em[r]["raw"]["files"]["h5ad"]["size"] for r in em)
    print(f"\n  {len(em)} regions, {human(total)} for all raw matrices")
    print(f"\n  Aliases accepted by --region and `get`: "
          f"{', '.join(sorted(REGION_ALIASES))}")
    print(f"\n  Pull one:  {pathlib.Path(sys.argv[0]).name} get "
          f"{DIRECTORY_10X} DFC-10X/raw")
    print(f"  Analyse:   python scripts/step2_qc_metrics.py --region PFC")
    return 0


def cmd_presets(_args):
    section("PRESETS")
    for name, p in PRESETS.items():
        print(f"  {name:16s} {p['description']}")
        print(f"  {'':16s} {len(p['files'])} files")
    print(f"\n  Usage: {pathlib.Path(sys.argv[0]).name} preset mtg")
    return 0


def expand_region_alias(pat):
    """Rewrite a leading region alias or anatomical word into its code.

    `get ... PFC/raw` -> `DFC-10X/raw`, and a bare `PFC` becomes
    `DFC-10X/*` so the plan shows both normalizations rather than
    guessing which one was wanted. Patterns whose first segment is not a
    region are returned untouched, so file names like `cell_metadata`
    and globs like `MTG*` behave exactly as before.
    """
    head, sep, tail = pat.partition("/")
    key = head.lower().removesuffix("-10x")

    if head in REGION_INFO:
        code = head
    elif f"{head}-10X" in REGION_INFO:
        code = f"{head}-10X"
    elif key in REGION_ALIASES:
        code = REGION_ALIASES[key]
    else:
        hits = [c for c, (desc, _) in REGION_INFO.items()
                if key and key in desc.lower()]
        code = hits[0] if len(hits) == 1 else None

    if code is None or code == head:
        return pat
    return f"{code}/{tail}" if (sep and tail) else f"{code}/*"


def resolve_selection(manifest, args):
    """Turn command-line selectors into a concrete list of file records."""
    if args.command == "preset":
        preset = PRESETS[args.name]
        wanted = preset["files"]
        index = {(r["directory"], r["file_name"]): r for r in catalog(manifest)}
        recs, missing = [], []
        for key in wanted:
            if key in index:
                recs.append(index[key])
            else:
                missing.append(key)
        if missing:
            raise SystemExit(
                f"\nPreset {args.name!r} names files absent from "
                f"{args.manifest}:\n  "
                + "\n  ".join(f"{d} / {f}" for d, f in missing)
                + "\nThe release may have reorganized. Check with `list`."
            )
        return recs, preset["description"]

    # `get`
    pool = catalog(manifest, directory=args.directory,
                   kinds=(args.kind,) if args.kind else KINDS)
    if not pool:
        known = sorted(manifest["file_listing"])
        raise SystemExit(
            f"\nNo files for directory {args.directory!r}"
            + (f" and kind {args.kind!r}" if args.kind else "")
            + ".\nKnown directories:\n  " + "\n  ".join(known)
        )

    if not args.files:
        if not args.kind:
            raise SystemExit(
                f"\nRefusing to fetch all {len(pool)} files "
                f"({human(sum(r['size'] for r in pool))}) in "
                f"{args.directory} without a narrower selection.\n"
                f"Name files, pass a glob, or restrict with "
                f"--kind {{{','.join(KINDS)}}}.\n"
                f"See what is there:  list {args.directory}"
            )
        return pool, f"all {args.kind} in {args.directory}"

    recs, unmatched = [], []
    for raw_pat in args.files:
        pat = expand_region_alias(raw_pat)
        if pat != raw_pat:
            print(f"  {raw_pat!r} -> {pat!r}  "
                  f"({region_description(pat.split('/')[0])})")
        hits = [r for r in pool if fnmatch.fnmatch(r["file_name"], pat)]
        if hits:
            recs.extend(hits)
        else:
            unmatched.append(raw_pat)
    if unmatched:
        raise SystemExit(
            f"\nNo file in {args.directory} matched: "
            + ", ".join(repr(u) for u in unmatched)
            + f"\nSee what is there:  list {args.directory}"
            f"\nBrain regions by name:  regions"
        )
    # De-duplicate while preserving order (overlapping globs).
    seen, uniq = set(), []
    for r in recs:
        key = (r["directory"], r["file_name"])
        if key not in seen:
            seen.add(key)
            uniq.append(r)
    return uniq, f"{len(uniq)} file(s) from {args.directory}"


def cmd_fetch(args):
    cache, manifest, mname = open_cache(args.manifest)
    recs, label = resolve_selection(manifest, args)

    # ---- Plan, with sizes, before anything transfers. Decision 2. ------
    section("PLAN")
    print(f"  manifest: {mname}")
    print(f"  selection: {label}\n")
    print(f"  {'directory':30s} {'file_name':40s} {'size':>11s}  state")
    print(f"  {'-' * 30} {'-' * 40} {'-' * 11:>11s}  {'-' * 7}")

    to_get, cached_bytes = [], 0
    for r in sorted(recs, key=lambda r: (r["directory"], r["file_name"])):
        state, have = local_state(r)
        desc = region_description(r["file_name"].split("/")[0])
        print(f"  {r['directory']:30s} {r['file_name']:40s} "
              f"{human(r['size']):>11s}  {state}"
              + (f"   {desc}" if desc else ""))
        if state == "cached":
            cached_bytes += have
        else:
            to_get.append(r)

    need = sum(r["size"] for r in to_get)
    print(f"\n  selected:    {len(recs)} files, {human(sum(r['size'] for r in recs))}")
    print(f"  already here: {len(recs) - len(to_get)} files, {human(cached_bytes)}")
    print(f"  TO DOWNLOAD:  {len(to_get)} files, {human(need)}")

    if not to_get:
        print("\n  Nothing to do — everything selected is already cached at")
        print("  the manifest's size. Use --force to re-download.")
        return 0

    if need > args.max_gb * 1e9:
        raise SystemExit(
            f"\nRefusing: {human(need)} exceeds --max-gb {args.max_gb}.\n"
            f"Narrow the selection, or raise the cap deliberately:\n"
            f"  --max-gb {need / 1e9:.0f}"
        )

    if args.dry_run:
        print("\n  --dry-run: nothing was downloaded.")
        return 0

    if need > args.confirm_gb * 1e9 and not args.yes:
        if not sys.stdin.isatty():
            raise SystemExit(
                f"\n{human(need)} needs confirmation but stdin is not a "
                f"terminal. Re-run with --yes."
            )
        print(f"\n  About to download {human(need)}.")
        if input("  Proceed? [y/N] ").strip().lower() not in ("y", "yes"):
            print("\n  Aborted. Nothing was downloaded.")
            return 1

    # ---- Fetch ---------------------------------------------------------
    section("FETCH")
    done = 0
    for i, r in enumerate(to_get, 1):
        print(f"  [{i}/{len(to_get)}] {r['directory']} / {r['file_name']} "
              f"({human(r['size'])}) ...", flush=True)
        # One kind-agnostic accessor (decision 4). Hash checking is left
        # at its default of on (decision 5).
        path = cache.get_file_path(
            directory=r["directory"], file_name=r["file_name"],
            force_download=args.force,
        )
        p = pathlib.Path(path)
        got = p.stat().st_size
        if got != r["size"]:
            raise SystemExit(
                f"\nSize mismatch for {r['file_name']}: got {got:,} bytes, "
                f"manifest says {r['size']:,}. Do not use this file."
            )
        done += got
        print(f"           -> {p}")

    section("DONE")
    print(f"  downloaded: {len(to_get)} files, {human(done)}")
    print(f"  manifest:   {mname}")
    print(f"  cache:      {CACHE_DIR}")
    if args.command == "preset" and args.name == "mtg":
        print("\n  Next:")
        print("    python scripts/step0_discover.py    # provenance (no download)")
        print("    python scripts/step1_inventory.py   # what is in the matrix")
        print("    python scripts/step2_qc_metrics.py  # ~5 min, writes results/")
        print("    python scripts/step3_annotation.py  # joins annotations")
    return 0


# ========================= cli ==========================================

def build_parser():
    ap = argparse.ArgumentParser(
        prog="step0b_fetch.py",
        description="Browse and fetch ABC Atlas data from the Allen S3 bucket.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  step0b_fetch.py preset mtg                        reproduce this project\n"
            "  step0b_fetch.py list                              all directories\n"
            "  step0b_fetch.py list SEA-AD-Multiregion-10X       files in one\n"
            "  step0b_fetch.py regions                           brain regions by name\n"
            "  step0b_fetch.py search PFC                        find files (anatomy aware)\n"
            "  step0b_fetch.py get SEA-AD-Multiregion-10X MTG-10X/raw\n"
            "  step0b_fetch.py get SEA-AD-Multiregion-10X --kind metadata\n"
            "  step0b_fetch.py get SEA-AD-Multiregion-10X 'HIP-10X/*' --dry-run\n"
        ),
    )

    def common(p):
        p.add_argument("--manifest", default=MANIFEST, metavar="REL",
                       help=f"release to use, or 'latest' (default: the project "
                            f"pin, {MANIFEST})")
        p.add_argument("--kind", choices=KINDS,
                       help="restrict to one data kind")

    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("list", help="list directories, or files in one")
    p.add_argument("directory", nargs="?", help="directory to list")
    common(p)
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("search", help="find files across the whole atlas")
    p.add_argument("pattern", help="glob, or a bare word (wrapped as *word*)")
    common(p)
    p.set_defaults(func=cmd_search)

    p = sub.add_parser("regions", help="list brain regions with anatomical names")
    p.add_argument("--manifest", default=MANIFEST, metavar="REL",
                   help=f"release to use, or 'latest' (default: {MANIFEST})")
    p.set_defaults(func=cmd_regions)

    p = sub.add_parser("presets", help="list the named reproducible file sets")
    p.set_defaults(func=cmd_presets)

    def fetch_flags(p):
        common(p)
        p.add_argument("--dry-run", action="store_true",
                       help="print the plan and exit without downloading")
        p.add_argument("--yes", "-y", action="store_true",
                       help="skip the confirmation prompt")
        p.add_argument("--force", action="store_true",
                       help="re-download even if the local hash matches")
        p.add_argument("--max-gb", type=float, default=50.0, metavar="N",
                       help="refuse transfers larger than this (default: 50)")
        p.add_argument("--confirm-gb", type=float, default=1.0, metavar="N",
                       help="prompt above this size (default: 1)")

    p = sub.add_parser("get", help="fetch named files from a directory")
    p.add_argument("directory")
    p.add_argument("files", nargs="*", metavar="FILE",
                   help="file names or globs; omit only with --kind")
    fetch_flags(p)
    p.set_defaults(func=cmd_fetch)

    p = sub.add_parser("preset", help="fetch a named reproducible file set")
    p.add_argument("name", choices=sorted(PRESETS))
    fetch_flags(p)
    p.set_defaults(func=cmd_fetch)

    return ap


def main():
    args = build_parser().parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
