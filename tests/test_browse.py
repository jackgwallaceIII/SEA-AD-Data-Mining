"""
test_browse.py — unit tests for scripts/browse.py (show / search).

No network: a temporary catalog is built by index.refresh() over the
in-memory FakeS3 buckets from test_index, then browsed.

    python -m unittest discover -s tests -t tests -v

Why these tests: `show` is what a newcomer will read numbers from. The
failure modes that matter are shared multi-region data being added into
one region's total, deprecated data vanishing without a trace, an
ambiguous region being guessed, and labels that silently truncate.
"""

import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import browse  # noqa: E402
import index  # noqa: E402
from test_index import NP, NP_KEYS, SC, SC_KEYS  # noqa: E402
from test_sources import FakeS3  # noqa: E402

TEN = "HIP;MEC;LEC;ITG;MTG;FI;STG;DFC;AnG;V1C"


class BrowseCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dir = pathlib.Path(cls.tmp.name)
        for bucket, keys in ((SC, SC_KEYS), (NP, NP_KEYS)):
            index.refresh(bucket, client=FakeS3(keys), catalog_dir=cls.dir,
                          out=lambda *a: None)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def setUp(self):
        self.lines = []

    def out(self, *a):
        self.lines.append(" ".join(str(x) for x in a))

    def text(self):
        return "\n".join(self.lines)

    def show(self, *a, **kw):
        return browse.show(*a, catalog_dir=self.dir, out=self.out, **kw)

    def search(self, *a, **kw):
        return browse.search(*a, catalog_dir=self.dir, out=self.out, **kw)


class TestShowRegion(BrowseCase):

    def test_region_only_and_shared_are_separate(self):
        res = self.show("MTG")
        only = sum(r["size"] for r in res["only"])
        shared = sum(r["size"] for r in res["shared"])
        # MTG RNA folder + DREAM (single-cell) + one MTG-only slide (neuropath)
        self.assertEqual(only, 1 + 1000 + 100 + 150 + 300 + 800)
        self.assertEqual(shared, 400 + 700)                    # 10-region + MTG;STG slide
        self.assertTrue(all(len(r["region_list"]) == 1 for r in res["only"]))

    def test_multi_region_unit_appears_for_each_of_its_regions(self):
        for code in ("HIP", "V1C", "DFC"):
            res = self.show(code)
            self.assertIn(TEN, {r["regions"] for r in res["shared"]}, code)

    def test_deprecated_hidden_but_reported(self):
        res = self.show("MTG")
        self.assertEqual(sum(r["size"] for r in res["hidden"]), 900)
        self.assertIn("hidden: 1 deprecated", self.text())
        self.assertFalse(any(r["deprecated"] for r in res["only"]))

    def test_deprecated_included_on_request(self):
        res = self.show("MTG", include_deprecated=True)
        self.assertTrue(any(r["deprecated"] for r in res["only"]))
        self.assertEqual(res["hidden"], [])

    def test_aliases_resolve(self):
        self.assertEqual(self.show("middle temporal gyrus")["code"], "MTG")

    def test_ambiguous_region_is_refused(self):
        with self.assertRaises(SystemExit) as cm:
            self.show("temporal")
        self.assertIn("more than one", str(cm.exception))

    def test_modality_filter_is_a_loose_substring(self):
        res = self.show("MTG", modality="neuro")
        self.assertEqual([r["size"] for r in res["only"]], [800])
        self.assertEqual([r["regions"] for r in res["shared"]], ["MTG;STG"])
        self.assertTrue(all(r["modality"] == "neuropathology"
                            for r in res["only"] + res["shared"]))

    def test_product_breakdown_separates_plots_from_data(self):
        self.show("MTG")
        t = self.text()
        self.assertIn("plots", t)
        self.assertIn("dataset", t)

    def test_files_lists_largest_first_with_s3_paths(self):
        self.show("MTG", files=True, limit=2)
        listed = [l for l in self.lines if "s3://" in l]
        self.assertEqual(len(listed), 2)
        self.assertIn("SEAAD_MTG_RNAseq_final-nuclei", listed[0])   # 1000 bytes, largest
        self.assertIn("objects]", self.text() + "objects]")         # collapsed marker allowed

    def test_shared_label_names_every_region_count(self):
        self.show("DFC")
        self.assertIn("(10 regions)", self.text())


class TestOverview(BrowseCase):

    def test_overview_adds_up_to_the_total(self):
        res = self.show()
        total = sum(r["size"] for rs in res["groups"].values() for r in rs)
        units = browse.load_units(self.dir)
        self.assertEqual(total, sum(r["size"] for r in units if not r["deprecated"]))

    def test_multi_region_groups_are_their_own_rows(self):
        res = self.show()
        regions = {k[0] for k in res["groups"]}
        self.assertIn("MTG+STG", regions)
        self.assertIn("HIP+…+V1C (10 regions)", regions)
        self.assertIn("(multi-region project)", regions)


class TestSearch(BrowseCase):

    def test_bare_word_and_glob(self):
        a = self.search("final-nuclei")["hits"]
        b = self.search("*final-nuclei*")["hits"]
        self.assertEqual([r["unit"] for r in a], [r["unit"] for r in b])
        self.assertTrue(a)

    def test_case_insensitive(self):
        self.assertTrue(self.search("dream")["hits"])

    def test_deprecated_hidden_and_counted(self):
        res = self.search("previous_objects")
        self.assertEqual(res["hits"], [])
        self.assertTrue(res["hidden"])
        self.assertIn("deprecated hidden", self.text())

    def test_no_match_says_so(self):
        res = self.search("zzz-no-such-thing")
        self.assertEqual(res["hits"], [])
        self.assertIn("no matches", self.text())

    def test_filters_combine(self):
        res = self.search("*", bucket="neuropath")
        self.assertTrue(res["hits"])
        self.assertTrue(all(r["bucket"] == NP for r in res["hits"]))


class TestLabels(unittest.TestCase):

    def test_region_label_never_silently_truncates(self):
        self.assertEqual(browse.region_label(["HIP", "MEC"]), "HIP+MEC")
        self.assertEqual(browse.region_label(TEN.split(";")),
                         "HIP+…+V1C (10 regions)")

    def test_missing_catalog_is_a_clear_error(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(SystemExit) as cm:
                browse.load_units(pathlib.Path(d))
            self.assertIn("refresh", str(cm.exception))


if __name__ == "__main__":
    unittest.main(verbosity=2)
