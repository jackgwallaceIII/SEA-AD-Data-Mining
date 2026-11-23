"""
test_plan.py — unit tests for scripts/plan.py (download preview).

No network and no downloads: a temporary catalog is built over the
in-memory FakeS3 buckets, and a temporary data_cache/ stands in for the
real one, so local-state checks never touch your disk.

    python -m unittest discover -s tests -t tests -v

Why these tests: the plan is the last thing a person reads before a
download. It must refuse a broad selection, never count a file that is
already present, spot a truncated file, enforce --max-gb, and tell real
duplicates (same content) from look-alikes (same name, different content).
"""

import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import index  # noqa: E402
import plan  # noqa: E402
from test_index import SC, SC_KEYS  # noqa: E402
from test_sources import FakeS3  # noqa: E402

MTG_H5AD = "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad"


class PlanCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self.tmp.name)
        self.cat, self.data = base / "catalog", base / "data_cache"
        index.refresh(SC, client=FakeS3(SC_KEYS), catalog_dir=self.cat,
                      out=lambda *a: None)
        self.lines = []

    def tearDown(self):
        self.tmp.cleanup()

    def out(self, *a):
        self.lines.append(" ".join(str(x) for x in a))

    def preview(self, *a, **kw):
        return plan.preview(*a, catalog_dir=self.cat, data_dir=self.data,
                            out=self.out, **kw)

    def text(self):
        return "\n".join(self.lines)


class TestSelection(PlanCase):

    def test_broad_selection_is_refused(self):
        with self.assertRaises(SystemExit) as cm:
            self.preview()
        self.assertIn("Refusing", str(cm.exception))

    def test_region_excludes_shared_units_unless_asked(self):
        p, _ = self.preview("MTG")
        self.assertTrue(all(len(u["region_list"]) == 1 for u in p["units"]))
        self.assertIn("shared with other regions", self.text())
        p2, _ = self.preview("MTG", include_shared=True)
        self.assertGreater(len(p2["units"]), len(p["units"]))

    def test_match_glob_and_word(self):
        a, _ = self.preview(match="final-nuclei")
        b, _ = self.preview(match="*final-nuclei*")
        self.assertEqual([u["unit"] for u in a["units"]], [u["unit"] for u in b["units"]])
        self.assertTrue(a["units"])

    def test_deprecated_left_out_and_reported(self):
        p, _ = self.preview("MTG")
        self.assertFalse(any(u["deprecated"] for u in p["units"]))
        self.assertIn("deprecated", self.text())

    def test_filters_narrow(self):
        p, _ = self.preview("MTG", product="metadata")
        self.assertEqual({u["product"] for u in p["units"]}, {"cell metadata table"})


class TestDestinationAndState(PlanCase):

    def test_destination_mirrors_the_bucket_layout(self):
        p, _ = self.preview(match=MTG_H5AD)
        self.assertEqual(p["units"][0]["dest"], self.data / SC / MTG_H5AD)

    def test_missing_partial_present(self):
        dest = self.data / SC / MTG_H5AD
        p, _ = self.preview(match=MTG_H5AD)
        self.assertEqual(p["units"][0]["state"], "missing")
        dest.parent.mkdir(parents=True)
        dest.write_bytes(b"x" * 10)                       # catalog says 1000
        p, _ = self.preview(match=MTG_H5AD)
        self.assertEqual(p["units"][0]["state"], "partial")
        self.assertEqual(p["need"], 1000)
        dest.write_bytes(b"x" * 1000)
        p, _ = self.preview(match=MTG_H5AD)
        self.assertEqual(p["units"][0]["state"], "present")
        self.assertEqual(p["need"], 0)

    def test_collapsed_unit_state_sums_its_folder(self):
        folder = self.data / SC / "MTG/RNAseq/pseudoprogression-plots"
        (folder / "a").mkdir(parents=True)
        (folder / "a" / "x.jpg").write_bytes(b"x" * 150)  # catalog says 150 total
        p, _ = self.preview(match="pseudoprogression-plots")
        self.assertEqual(p["units"][0]["state"], "present")

    def test_nothing_is_written_to_data_cache(self):
        self.preview("MTG", include_shared=True)
        self.assertFalse(self.data.exists())


class TestGuards(PlanCase):

    def test_max_gb_refuses_and_names_the_needed_cap(self):
        p, code = self.preview("MTG", max_gb=1e-6)        # ~1 kB cap
        self.assertTrue(p["over_cap"])
        self.assertEqual(code, 2)
        self.assertIn("REFUSED", self.text())
        self.assertIn("--max-gb", self.text())

    def test_under_cap_is_accepted(self):
        _, code = self.preview("MTG")
        self.assertEqual(code, 0)

    def test_get_without_dry_run_and_without_yes_does_not_download(self):
        # No terminal to ask and no --yes: refuse. stdin is replaced by a
        # non-terminal so the test behaves the same whether it is launched from
        # a terminal (where the real stdin would be asked y/N, invisibly) or not.
        import contextlib
        import io
        from unittest import mock
        saved = (plan.DATA_DIR, index.CATALOG_DIR)
        plan.DATA_DIR, index.CATALOG_DIR = self.data, self.cat
        try:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), mock.patch("sys.stdin", io.StringIO()):
                args = index.build_parser().parse_args(["get", "--match", MTG_H5AD])
                code = args.func(args)
        finally:
            plan.DATA_DIR, index.CATALOG_DIR = saved
        self.assertEqual(code, 6)
        self.assertIn("needs your confirmation", buf.getvalue())
        self.assertFalse(self.data.exists())


class TestWarnings(unittest.TestCase):

    def row(self, unit, size, etag, kind="file", n=1, **kw):
        base = dict(unit=unit, size=size, etag=etag, kind=kind, n_objects=n,
                    deprecated=False, missing="", format="csv", product="x",
                    bucket=SC, size_is_minimum=False)
        return {**base, **kw}

    def test_same_etag_is_identical_even_with_different_names(self):
        rows = [self.row("a/one.csv", 10, "abc"), self.row("b/two.csv", 10, "abc")]
        identical, likely = plan.duplicate_groups(rows)
        self.assertEqual(len(identical), 1)
        self.assertEqual(likely, [])

    def test_same_name_size_but_different_md5_is_not_a_duplicate(self):
        rows = [self.row("a/x.html", 10, "aaa"), self.row("b/x.html", 10, "bbb")]
        self.assertEqual(plan.duplicate_groups(rows), ([], []))

    def test_same_name_size_multipart_etags_is_probable(self):
        # The real MERFISH case: -95 vs -60 parts.
        rows = [self.row("a/SEAAD_MTG_MERFISH.h5ad", 497737587, "a129-95"),
                self.row("b/SEAAD_MTG_MERFISH.h5ad", 497737587, "0b36-60")]
        identical, likely = plan.duplicate_groups(rows)
        self.assertEqual((len(identical), len(likely)), (0, 1))

    def test_many_object_unit_is_flagged(self):
        rows = [self.row("p/", 10, "", kind="collapsed", n=1_134_632, product="plots")]
        self.assertTrue(any("1,134,632" in w for w in plan.warnings_for(rows)))

    def test_pickle_is_flagged(self):
        rows = [self.row("x/menr.pkl", 10, "e", format="pickle")]
        self.assertTrue(any("pickle" in w for w in plan.warnings_for(rows)))


class TestGitignore(unittest.TestCase):

    def test_data_cache_is_gitignored(self):
        text = (pathlib.Path(__file__).resolve().parent.parent / ".gitignore").read_text()
        self.assertIn("data_cache/", text.splitlines())


if __name__ == "__main__":
    unittest.main(verbosity=2)
