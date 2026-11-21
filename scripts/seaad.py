"""
seaad.py — shared configuration and path resolution for the SEA-AD pipeline.

Not a step. Imported by step0/step1/step2/step3 so that the release pin,
the region selection and the output naming exist in exactly one place.

Why this module exists, given that every step is otherwise a standalone
script with its reasoning in its own docstring:

1. THE RELEASE PIN MUST NOT DRIFT BETWEEN STEPS.
   The pin was duplicated in three scripts and about to be duplicated in
   four. If one copy is updated and another is not, two steps read
   different releases and the mismatch is invisible — results/ would hold
   a table joined from two versions of the same data. One constant, one
   place.

2. REGION IS A PARAMETER, NOT A HARD-CODED PATH.
   The multiregion release covers ten regions that are structurally
   identical: same 36,601 genes, same taxonomy, same donor tables, same
   CSR layout, same manifest version. Only the matrix and the
   feature_matrix_label differ. Every step therefore takes --region and
   defaults to MTG-10X, which is what this project analyses and why the
   rationale in step0_discover.py argues for it.

3. THE VERSION DIRECTORY IS RESOLVED, NEVER TYPED.
   Paths were hard-coded as .../20260630/MTG-10X-raw.h5ad. That string
   silently ties the code to one release, which defeats the pin: change
   MANIFEST and the path still points at the old version directory.
   h5ad_path() reads relative_path out of the pinned manifest instead, so
   the file location follows the pin automatically.

4. OUTPUTS ARE NAMED PER REGION.
   results/step2_qc_metrics.csv would be overwritten by the second region
   analysed, and step3 would silently join one region's QC metrics to
   another region's annotations. Output names carry the region.
"""

import json
import pathlib
import warnings

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

try:
    from abc_atlas_access.abc_atlas_cache.cloud_cache import (
        OutdatedManifestWarning,
    )
    # Staying on an older release is the point of the pin. open_cache()
    # reports the same fact in one line.
    warnings.filterwarnings("ignore", category=OutdatedManifestWarning)
except ImportError:  # pragma: no cover - older abc_atlas_access
    pass

# scripts/ is one level below the project root.
ROOT = pathlib.Path(__file__).resolve().parent.parent
CACHE_DIR = ROOT / "abc_cache"
RESULTS = ROOT / "results"

# ---- The pin. The single source of truth for every step. ---------------
# step0b_fetch.py carries the same value because it must run before this
# module's cache exists; if they ever disagree, step0b_fetch.py wins.
MANIFEST = "releases/20261101/manifest.json"

DIRECTORY = "SEA-AD-Multiregion-10X"
TAXONOMY = "SEA-AD-Multiregion-taxonomy"

DEFAULT_REGION = "MTG-10X"
DEFAULT_NORMALIZATION = "raw"

# Anatomical names for the region codes, so "PFC" is findable without
# knowing that SEA-AD calls it DFC. Taken verbatim from the release's own
# metadata/SEA-AD-Multiregion-10X/value_sets.csv, rows where
# field == 'region_of_interest_label', with their DHBA ontology ids.
# Hard-coded rather than read from that CSV because `step0b_fetch.py
# list` has to work on a cold cache, before any metadata is fetched.
REGION_INFO = {
    "HIP-10X": ("hippocampus (hippocampal formation)", "DHBA:10294"),
    "MEC-10X": ("medial (caudal) entorhinal cortex", "DHBA:10319"),
    "LEC-10X": ("lateral (rostral) entorhinal cortex", "DHBA:10318"),
    "ITG-10X": ("inferior temporal gyrus", "DHBA:12142"),
    "MTG-10X": ("middle temporal gyrus", "DHBA:12141"),
    "FI-10X": ("agranular frontal insular cortex (area FI)", "DHBA:10329"),
    "STG-10X": ("superior temporal gyrus", "DHBA:12140"),
    "DFC-10X": ("dorsolateral prefrontal cortex", "DHBA:10173"),
    "AnG-10X": ("angular gyrus", "DHBA:12136"),
    "V1C-10X": ("primary visual cortex (striate, area 17)", "DHBA:10269"),
    # Head of caudate ships in the separate SEA-AD-CaH-10X directory, not
    # in the multiregion release; listed here so the name resolves.
    "CaH-10X": ("head of caudate", "DHBA:10335"),
}

# Common names that are not the release's own codes. Only unambiguous
# ones belong here: "temporal" and "entorhinal" each match more than one
# region, so they are left to the description search, which reports the
# ambiguity instead of silently picking one.
REGION_ALIASES = {
    "pfc": "DFC-10X",
    "dlpfc": "DFC-10X",
    "prefrontal": "DFC-10X",
    "dpfc": "DFC-10X",
    "hippocampus": "HIP-10X",
    "v1": "V1C-10X",
    "visual": "V1C-10X",
    "insula": "FI-10X",
    "insular": "FI-10X",
    "angular": "AnG-10X",
    "caudate": "CaH-10X",
}


def region_description(region):
    """Anatomical name for a region code, or '' if unknown."""
    return REGION_INFO.get(region, ("", ""))[0]


def section(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


# ========================= manifest / cache =============================

def open_cache(manifest_name=MANIFEST, quiet=False):
    """Return (cache, manifest_dict) with the release pinned.

    from_s3_cache() otherwise loads whichever manifest was last used in
    this cache directory, falling back to the newest on S3 — both mutable
    state that can change the data under an analysis with no code change.
    """
    cache = AbcProjectCache.from_s3_cache(str(CACHE_DIR))
    try:
        cache.load_manifest(manifest_name)
    except Exception as exc:
        raise SystemExit(
            f"\nCould not load manifest {manifest_name!r}: {exc}\n"
            f"Available: {cache.list_manifest_file_names}\n"
            f"Fetch the data first:  python scripts/step0b_fetch.py preset mtg"
        )
    if cache.current_manifest != manifest_name:
        raise SystemExit(
            f"\nPin failed: asked for {manifest_name}, "
            f"got {cache.current_manifest}."
        )

    path = CACHE_DIR / manifest_name
    if not path.exists():
        raise SystemExit(f"Manifest not on disk after load: {path}")
    manifest = json.loads(path.read_text())

    if not quiet:
        print(f"  manifest: {cache.current_manifest}")
        if cache.latest_manifest_file != manifest_name:
            print(f"  (newer release exists: {cache.latest_manifest_file} — "
                  f"staying on the pin)")
    return cache, manifest


def available_regions(manifest):
    """Region labels in the pinned release, e.g. ['AnG-10X', ...]."""
    return sorted(manifest["file_listing"][DIRECTORY]["expression_matrices"])


def _region_table(regions):
    """Available regions with their anatomical names, for error messages."""
    return "\n".join(
        f"  {r:9s} {region_description(r)}" for r in regions)


def resolve_region(manifest, region):
    """Resolve a region code, bare name, alias or anatomical description.

    Accepts, in order: the exact label ('MTG-10X'), the bare code
    ('MTG'), either in any case, a common alias ('PFC' -> 'DFC-10X'),
    or a word from the anatomical name ('prefrontal', 'hippocamp').

    Fails with the full region table rather than a KeyError three frames
    deep: a wrong region here silently analyses the wrong brain area, and
    the codes are opaque enough that this is the likeliest user error.
    An ambiguous description is refused rather than guessed.
    """
    regions = available_regions(manifest)
    if region in regions:
        return region
    candidate = f"{region}-10X"
    if candidate in regions:
        return candidate

    key = region.lower().removesuffix("-10x")
    for r in regions:
        if r.lower().removesuffix("-10x") == key:
            return r

    alias = REGION_ALIASES.get(key)
    if alias:
        if alias in regions:
            return alias
        # A known name, but not in this directory's release.
        raise SystemExit(
            f"\n{region!r} is {alias} ({region_description(alias)}), which "
            f"is not in {DIRECTORY}.\n"
            f"Available ({len(regions)}):\n{_region_table(regions)}"
        )

    # Last resort: match a word against the anatomical names.
    hits = [r for r in regions if key in region_description(r).lower()]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise SystemExit(
            f"\n{region!r} matches more than one region; be specific:\n"
            f"{_region_table(hits)}"
        )

    raise SystemExit(
        f"\nUnknown region {region!r} in {MANIFEST}.\n"
        f"Available ({len(regions)}):\n{_region_table(regions)}\n"
        f"The '-10X' suffix is optional. Aliases: "
        f"{', '.join(sorted(REGION_ALIASES))}"
    )


def matrix_entry(manifest, region, normalization=DEFAULT_NORMALIZATION):
    """The pinned manifest record for one region's expression matrix."""
    em = manifest["file_listing"][DIRECTORY]["expression_matrices"]
    try:
        return em[region][normalization]["files"]["h5ad"]
    except KeyError:
        raise SystemExit(
            f"\nNo {normalization!r} matrix for region {region!r}.\n"
            f"Normalizations available: {sorted(em.get(region, {}))}"
        )


def h5ad_path(manifest, region, normalization=DEFAULT_NORMALIZATION,
              require=True):
    """Local path to a region's matrix, resolved from the pinned manifest.

    The version directory is never typed out — it comes from the
    manifest's relative_path, so the location follows the pin.
    """
    entry = matrix_entry(manifest, region, normalization)
    path = CACHE_DIR / entry["relative_path"]
    if require and not path.exists():
        raise SystemExit(
            f"\nNot found: {path}\n"
            f"({entry['size'] / 1e9:.2f} GB, version {entry['version']})\n"
            f"Fetch it:\n"
            f"  python scripts/step0b_fetch.py get {DIRECTORY} "
            f"{region}/{normalization}"
        )
    if require and path.stat().st_size != entry["size"]:
        raise SystemExit(
            f"\nWrong size: {path}\n"
            f"  on disk:  {path.stat().st_size:,} bytes\n"
            f"  manifest: {entry['size']:,} bytes\n"
            f"Likely a truncated download. Re-fetch:\n"
            f"  python scripts/step0b_fetch.py get {DIRECTORY} "
            f"{region}/{normalization} --force"
        )
    return path


# ========================= outputs ======================================

def out_path(stem, region, ext, mkdir=True):
    """results/<stem>_<region>.<ext> — region-scoped so runs never collide."""
    if mkdir:
        RESULTS.mkdir(parents=True, exist_ok=True)
    return RESULTS / f"{stem}_{region}.{ext}"


def qc_csv(region):
    """step2's output / step3's input."""
    return out_path("step2_qc_metrics", region, "csv")


def annotation_csv(region):
    """step3's output."""
    return out_path("step3_annotation", region, "csv")


# ========================= cli ==========================================

def add_region_args(parser, normalization=False):
    """Attach --region (and optionally --normalization) to a step's parser."""
    parser.add_argument(
        "--region", default=DEFAULT_REGION, metavar="REGION",
        help=f"brain region to analyse; the '-10X' suffix is optional "
             f"(default: {DEFAULT_REGION})",
    )
    if normalization:
        parser.add_argument(
            "--normalization", default=DEFAULT_NORMALIZATION,
            choices=("raw", "log2"),
            help=f"matrix normalization; raw counts are required for QC "
                 f"metrics (default: {DEFAULT_NORMALIZATION})",
        )
    return parser
