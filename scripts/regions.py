"""
regions.py — one brain-region vocabulary shared by every data source.

Not a step. The catalog (and later the fetch and inspect commands) use it
to answer "which brain region is this?" for any name a bucket, the ABC
Atlas or a user might use.

    python scripts/regions.py                 # print the vocabulary
    python scripts/regions.py PFC             # resolve one name
    python scripts/regions.py --exact Caudate_Nucleus

Data lives in catalog/regions.csv (one row per region) and
catalog/region_groups.csv (folder names that cover more than one region).
Editing those files changes the vocabulary; no code change is needed.

Decisions made here, and why:

1. ONE VOCABULARY, CANONICAL CODES WITHOUT THE "-10X" SUFFIX.
   The same region is spelled four ways across the sources: MTG (single-
   cell bucket), MTG-10X (ABC Atlas), middle-temporal-gyrus (spatial and
   neuropathology buckets) and DHBA:12141 (ontology id). If each source
   kept its own spelling, "show me MTG" would silently miss three of
   them. The canonical code is the bare SEA-AD code (MTG, DFC, CaH ...);
   "-10X" names an ABC data product, not a region, so it is stripped.

2. NAMES ARE NORMALISED, NOT LISTED EXHAUSTIVELY.
   Case, '-', '_' and spaces are folded together, and a trailing "-10X"
   is dropped, so 'middle-temporal-gyrus', 'Middle_Temporal_Gyrus' and
   'middle temporal gyrus' are one key. Each region's anatomical name,
   with and without its parenthetical, is automatically an alias. The
   aliases column only holds names that cannot be derived (PFC, v1 ...).

3. TWO LOOKUPS WITH DIFFERENT STRICTNESS.
   codes_for() is EXACT and is what classifies bucket folder names. It
   never fuzzy-matches: a folder called 'temporal-lobe-misc' must come
   back as unknown (and be reviewed), not quietly filed under MTG.
   resolve() is for a person typing a region; it also accepts a word
   from the anatomical name ('hippocamp', 'dorsolateral') but refuses
   anything matching more than one region ('temporal', 'entorhinal').
   A wrong guess here means analysing the wrong brain area.

4. MULTI-REGION FOLDERS RETURN EVERY REGION THEY COVER.
   The neuropathology bucket has folders such as
   middle-temporal-gyrus-and-superior-temporal-gyrus. Forcing that into
   one code would hide STG data from an STG search. These live in
   region_groups.csv and codes_for() returns all their codes. resolve()
   refuses them, because a user choosing "a region" means one region.

5. UNKNOWN STAYS UNKNOWN UNTIL IT CAN BE VERIFIED.
   A region is added only with a real DHBA id. Example: the folder
   primary-visual-cortex-extrastriate-occipital-cortex sat unclassified
   (877 GB) until 2026-11-22, when SEA-AD's own code 'ESOC' (from the
   file V1C-ESOC_all-values_...csv) was confirmed against the Allen DHBA
   ontology as DHBA:10270, 'extrastriate occipital cortex'. Until a name
   can be checked like that, it stays 'unclassified' for a human.

6. THE FILES ARE VALIDATED ON LOAD.
   A key that maps to two different regions, a duplicate code, a
   malformed DHBA id, or a group naming an unknown code raises at load
   time. Collisions are the dangerous case: the last one written would
   win silently.

7. ERRORS ARE EXCEPTIONS, NOT SystemExit.
   seaad.py exits the program on a bad region, which suits a script.
   This module will also sit under a catalog command and possibly a GUI,
   which need to catch the error and show it, so it raises RegionError.
   seaad.py is left unchanged; tests check that the two agree.
"""

import csv
import dataclasses
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
REGIONS_CSV = ROOT / "catalog" / "regions.csv"
GROUPS_CSV = ROOT / "catalog" / "region_groups.csv"

_DHBA = re.compile(r"^DHBA:\d+$")


class RegionError(ValueError):
    """A region name that is unknown, ambiguous, or a bad vocabulary file."""


@dataclasses.dataclass(frozen=True)
class Region:
    code: str
    name: str
    dhba_id: str
    aliases: tuple
    notes: str


def normalise(text):
    """Fold case and separators; drop a trailing '-10X'.

    'MTG-10X' -> 'mtg', 'middle-temporal-gyrus' -> 'middle temporal gyrus',
    'Caudate_Nucleus' -> 'caudate nucleus'.
    """
    s = re.sub(r"[\s_\-]+", " ", str(text).strip().lower()).strip()
    return re.sub(r" 10x$", "", s)


def _strip_parenthetical(name):
    """'hippocampus (hippocampal formation)' -> 'hippocampus'."""
    return re.sub(r"\s*\(.*?\)", "", name).strip()


class Vocabulary:
    """The loaded, validated vocabulary. Build with load()."""

    def __init__(self, regions, groups):
        self.regions = regions          # code -> Region
        self.groups = groups            # normalised alias -> tuple of codes
        self.group_names = {}           # normalised alias -> spelling as written
        self._index = {}                # normalised key -> code
        for r in regions.values():
            keys = {r.code, r.dhba_id, r.name, _strip_parenthetical(r.name),
                    *r.aliases}
            for k in keys:
                self._add(normalise(k), r.code)
        for key, codes in groups.items():
            if key in self._index:
                raise RegionError(
                    f"group alias {key!r} is also a single-region name "
                    f"({self._index[key]})")

    def _add(self, key, code):
        if not key:
            return
        prev = self._index.get(key)
        if prev and prev != code:
            raise RegionError(
                f"name {key!r} maps to both {prev} and {code}; "
                f"fix catalog/regions.csv")
        self._index[key] = code

    # ---- lookups --------------------------------------------------------

    def codes_for(self, text):
        """EXACT lookup for folder/file names. Tuple of codes, or None.

        Never fuzzy. Multi-region folders return every code they cover.
        """
        key = normalise(text)
        if key in self._index:
            return (self._index[key],)
        if key in self.groups:
            return self.groups[key]
        return None

    def resolve(self, text):
        """Lenient lookup for a person choosing ONE region. Returns a code.

        Exact names first, then a unique word match against anatomical
        names. Ambiguous, multi-region and unknown input raise RegionError
        with the candidates listed.
        """
        exact = self.codes_for(text)
        if exact and len(exact) == 1:
            return exact[0]
        if exact:
            raise RegionError(
                f"{text!r} covers more than one region; pick one:\n"
                + self.table(exact))

        key = normalise(text)
        hits = tuple(c for c, r in self.regions.items()
                     if key and key in normalise(r.name))
        if len(hits) == 1:
            return hits[0]
        if hits:
            raise RegionError(
                f"{text!r} matches more than one region; be specific:\n"
                + self.table(hits))
        raise RegionError(
            f"unknown region {text!r}. Known regions:\n" + self.table())

    def table(self, codes=None):
        codes = codes or sorted(self.regions)
        return "\n".join(
            f"  {c:5s} {self.regions[c].name:44s} {self.regions[c].dhba_id}"
            for c in codes)


def load(regions_csv=REGIONS_CSV, groups_csv=GROUPS_CSV):
    """Read and validate the two CSVs."""
    regions = {}
    with open(regions_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            code = row["code"].strip()
            if not code:
                raise RegionError(f"empty code in {regions_csv}")
            if code in regions:
                raise RegionError(f"duplicate code {code!r} in {regions_csv}")
            dhba = row["dhba_id"].strip()
            if not _DHBA.match(dhba):
                raise RegionError(
                    f"{code}: dhba_id {dhba!r} is not of the form DHBA:<digits>")
            aliases = tuple(a.strip() for a in row["aliases"].split(";")
                            if a.strip())
            regions[code] = Region(code, row["name"].strip(), dhba, aliases,
                                   row.get("notes", "").strip())

    groups, group_names = {}, {}
    if pathlib.Path(groups_csv).exists():
        with open(groups_csv, newline="") as fh:
            for row in csv.DictReader(fh):
                codes = tuple(c.strip() for c in row["codes"].split(";")
                              if c.strip())
                unknown = [c for c in codes if c not in regions]
                if unknown:
                    raise RegionError(
                        f"group {row['alias']!r} names unknown code(s) "
                        f"{unknown}")
                if len(codes) < 2:
                    raise RegionError(
                        f"group {row['alias']!r} needs at least two codes; "
                        f"a single region belongs in regions.csv aliases")
                key = normalise(row["alias"])
                if key in groups:
                    raise RegionError(f"duplicate group alias {row['alias']!r}")
                groups[key] = codes
                group_names[key] = row["alias"].strip()
    vocab = Vocabulary(regions, groups)
    vocab.group_names = group_names
    return vocab


_DEFAULT = None


def default():
    """The vocabulary from catalog/, loaded once."""
    global _DEFAULT
    if _DEFAULT is None:
        _DEFAULT = load()
    return _DEFAULT


def main(argv):
    vocab = default()
    if not argv:
        print(f"{len(vocab.regions)} regions  ({REGIONS_CSV.relative_to(ROOT)})")
        print(vocab.table())
        print(f"\n{len(vocab.groups)} multi-region folder names  "
              f"({GROUPS_CSV.relative_to(ROOT)})")
        for alias, codes in vocab.groups.items():
            shown = vocab.group_names.get(alias, alias)
            print(f"  {shown:52s} -> {', '.join(codes)}")
        return 0
    exact = argv[0] == "--exact"
    for text in argv[1:] if exact else argv:
        try:
            if exact:
                codes = vocab.codes_for(text)
                print(f"{text!r} -> {', '.join(codes) if codes else 'UNKNOWN (unclassified)'}")
            else:
                code = vocab.resolve(text)
                print(f"{text!r} -> {code}  ({vocab.regions[code].name})")
        except RegionError as exc:
            print(f"{text!r} -> refused: {exc}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
