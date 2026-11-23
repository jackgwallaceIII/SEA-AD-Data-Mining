"""
test_regions.py — unit tests for scripts/regions.py (the region vocabulary).

No network, no data. The folder names below were taken from read-only
listings of the three open SEA-AD buckets on 2026-11-21, so these tests
pin the classification of the real layout, not an imagined one.

    python -m unittest discover -s tests -t tests -v

Why these tests: the vocabulary decides which brain region every catalog
row is filed under. A wrong answer files data under the wrong region and
nothing downstream would notice, so the tests cover (a) every real folder
name, (b) that exact lookup never guesses, (c) that ambiguity is refused,
and (d) that a broken vocabulary file fails at load instead of silently.
"""

import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import regions  # noqa: E402


def write(dir_, name, text):
    p = pathlib.Path(dir_) / name
    p.write_text(text)
    return p


HEADER = "code,name,dhba_id,aliases,notes\n"
GHEADER = "alias,codes,seen_in,notes\n"


class TestShippedVocabulary(unittest.TestCase):
    """The real catalog/regions.csv and region_groups.csv."""

    @classmethod
    def setUpClass(cls):
        cls.v = regions.load()

    def test_loads_and_has_twelve_regions(self):
        self.assertEqual(
            sorted(self.v.regions),
            ["AnG", "CaH", "DFC", "ESOC", "FI", "HIP", "ITG", "LEC", "MEC",
             "MTG", "STG", "V1C"])

    def test_esoc_is_the_dhba_structure(self):
        # Confirmed 2026-11-22 against the Allen DHBA ontology (graph 16).
        r = self.v.regions["ESOC"]
        self.assertEqual((r.dhba_id, r.name),
                         ("DHBA:10270", "extrastriate occipital cortex"))

    def test_every_alias_resolves_to_its_own_code(self):
        for code, r in self.v.regions.items():
            for name in (code, r.name, r.dhba_id, *r.aliases):
                self.assertEqual(self.v.codes_for(name), (code,), name)

    def test_one_region_four_spellings(self):
        # The core problem: the same region named four ways across sources.
        for name in ("MTG", "mtg", "MTG-10X", "middle-temporal-gyrus",
                     "Middle_Temporal_Gyrus", "middle temporal gyrus",
                     "DHBA:12141"):
            self.assertEqual(self.v.codes_for(name), ("MTG",), name)


class TestRealBucketFolders(unittest.TestCase):
    """Every region-like folder seen in the three open buckets."""

    @classmethod
    def setUpClass(cls):
        cls.v = regions.load()

    def check(self, folders):
        for folder, expected in folders.items():
            self.assertEqual(self.v.codes_for(folder), expected, folder)

    def test_single_cell_profiling(self):
        self.check({
            "AnG": ("AnG",), "FI": ("FI",), "HIP": ("HIP",),
            "ITG": ("ITG",), "LEC": ("LEC",), "MEC": ("MEC",),
            "MTG": ("MTG",), "STG": ("STG",), "V1C": ("V1C",),
            "DLPFC": ("DFC",), "PFC": ("DFC",),
            "Caudate_Nucleus": ("CaH",),
        })

    def test_spatial_transcriptomics(self):
        self.check({
            "caudate_nucleus": ("CaH",),
            "hippocampus": ("HIP",),
            "medial-entorhinal-cortex": ("MEC",),
            "middle-temporal-gyrus": ("MTG",),
        })

    def test_quantitative_neuropathology(self):
        self.check({
            "dorsolateral-prefrontal-cortex": ("DFC",),
            "middle-temporal-gyrus": ("MTG",),
            "middle-temporal-gyrus-and-superior-temporal-gyrus": ("MTG", "STG"),
            "hippocampus-medial-entorhinal-cortex": ("HIP", "MEC"),
        })

    def test_abc_atlas_labels(self):
        for code in ("AnG", "DFC", "FI", "HIP", "ITG", "LEC", "MEC", "MTG",
                     "STG", "V1C", "CaH"):
            self.assertEqual(self.v.codes_for(f"{code}-10X"), (code,))

    def test_non_region_folders_are_unknown(self):
        # Projects, modalities, viewers and scratch space are not regions.
        # Returning None is what makes them 'unclassified' rather than
        # quietly filed under a guessed region.
        for folder in ("DREAM", "Multiregion_2026",
                       "Microglia-and-Immune-for-AAIC", "RNAseq", "ATACseq",
                       "full-tabular-data", "simple_viewer", "temp",
                       "previous_objects", "combined_anndata_object",
                       "H20.33.001"):
            self.assertIsNone(self.v.codes_for(folder), folder)

    def test_extrastriate_folder_is_two_regions(self):
        # Mapping it to V1C alone would hide that half of every slide is
        # extrastriate cortex. ESOC was added once its DHBA id was found.
        self.assertEqual(self.v.codes_for(
            "primary-visual-cortex-extrastriate-occipital-cortex"), ("V1C", "ESOC"))


class TestExactNeverGuesses(unittest.TestCase):

    def setUp(self):
        self.v = regions.load()

    def test_partial_words_are_not_exact_matches(self):
        for word in ("temporal", "hippocamp", "dorsolateral", "cortex",
                     "entorhinal"):
            self.assertIsNone(self.v.codes_for(word), word)


class TestResolveForPeople(unittest.TestCase):

    def setUp(self):
        self.v = regions.load()

    def test_aliases(self):
        for q in ("PFC", "pfc", "DLPFC", "prefrontal", "DFC-10X", "A9", "a9"):
            self.assertEqual(self.v.resolve(q), "DFC", q)

    def test_unique_word_from_anatomical_name(self):
        self.assertEqual(self.v.resolve("dorsolateral"), "DFC")
        self.assertEqual(self.v.resolve("hippocamp"), "HIP")

    def test_ambiguous_word_is_refused_with_candidates(self):
        for q, expected in (("temporal", ("ITG", "MTG", "STG")),
                            ("entorhinal", ("LEC", "MEC"))):
            with self.assertRaises(regions.RegionError) as cm:
                self.v.resolve(q)
            msg = str(cm.exception)
            self.assertIn("more than one", msg)
            for code in expected:
                self.assertIn(code, msg)

    def test_multi_region_folder_is_refused_for_a_single_choice(self):
        with self.assertRaises(regions.RegionError) as cm:
            self.v.resolve("hippocampus-medial-entorhinal-cortex")
        self.assertIn("HIP", str(cm.exception))
        self.assertIn("MEC", str(cm.exception))

    def test_unknown_lists_every_region(self):
        with self.assertRaises(regions.RegionError) as cm:
            self.v.resolve("cerebellum")
        msg = str(cm.exception)
        for code in self.v.regions:
            self.assertIn(code, msg)


class TestBrokenFilesFailAtLoad(unittest.TestCase):
    """A bad vocabulary must fail loudly, never load half-right."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.groups = write(self.tmp.name, "g.csv", GHEADER)

    def tearDown(self):
        self.tmp.cleanup()

    def load(self, body, groups=None):
        r = write(self.tmp.name, "r.csv", HEADER + body)
        g = write(self.tmp.name, "g.csv", GHEADER + groups) if groups else self.groups
        return regions.load(r, g)

    def test_alias_shared_by_two_regions(self):
        with self.assertRaises(regions.RegionError) as cm:
            self.load("AAA,alpha,DHBA:1,shared,\nBBB,beta,DHBA:2,shared,\n")
        self.assertIn("shared", str(cm.exception))

    def test_duplicate_code(self):
        with self.assertRaises(regions.RegionError):
            self.load("AAA,alpha,DHBA:1,,\nAAA,beta,DHBA:2,,\n")

    def test_malformed_dhba(self):
        with self.assertRaises(regions.RegionError):
            self.load("AAA,alpha,12141,,\n")

    def test_group_with_unknown_code(self):
        with self.assertRaises(regions.RegionError):
            self.load("AAA,alpha,DHBA:1,,\n", groups="x-and-y,AAA;ZZZ,,\n")

    def test_group_with_one_code(self):
        with self.assertRaises(regions.RegionError):
            self.load("AAA,alpha,DHBA:1,,\n", groups="x,AAA,,\n")

    def test_group_alias_colliding_with_a_region_name(self):
        with self.assertRaises(regions.RegionError):
            self.load("AAA,alpha,DHBA:1,,\nBBB,beta,DHBA:2,,\n",
                      groups="alpha,AAA;BBB,,\n")

    def test_valid_minimal_file_loads(self):
        v = self.load("AAA,alpha one (extra),DHBA:1,a1,\n")
        self.assertEqual(v.codes_for("Alpha_One"), ("AAA",))
        self.assertEqual(v.codes_for("alpha one (extra)"), ("AAA",))


def seaad_tables():
    """REGION_INFO and REGION_ALIASES read from seaad.py's source text.

    Parsed with ast rather than imported: importing seaad.py pulls in
    abc_atlas_access, so in an environment without it this drift guard
    would error out (or, if skipped, silently stop guarding). Reading the
    two literal dicts needs only the stdlib, so the check always runs.
    """
    import ast
    src = (pathlib.Path(__file__).resolve().parent.parent
           / "scripts" / "seaad.py").read_text()
    found = {}
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            name = getattr(node.targets[0], "id", None)
            if name in ("REGION_INFO", "REGION_ALIASES"):
                found[name] = ast.literal_eval(node.value)
    return found["REGION_INFO"], found["REGION_ALIASES"]


class TestAgreesWithSeaadPy(unittest.TestCase):
    """seaad.py keeps its own region table for the ABC pipeline.

    It is deliberately left unchanged in this step, so this test is the
    guard against the two tables drifting apart.
    """

    def test_codes_names_and_dhba_ids_match(self):
        info, _ = seaad_tables()
        v = regions.load()
        for label, (name, dhba) in info.items():
            code = label.removesuffix("-10X")
            self.assertIn(code, v.regions, label)
            self.assertEqual(v.regions[code].name, name, label)
            self.assertEqual(v.regions[code].dhba_id, dhba, label)

    def test_seaad_aliases_resolve_the_same_way(self):
        _, aliases = seaad_tables()
        v = regions.load()
        for alias, label in aliases.items():
            self.assertEqual(v.resolve(alias), label.removesuffix("-10X"), alias)


if __name__ == "__main__":
    unittest.main(verbosity=2)
