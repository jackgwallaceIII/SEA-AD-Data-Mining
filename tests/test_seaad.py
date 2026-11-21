"""
test_seaad.py — unit tests for scripts/seaad.py (shared config and paths).

Tier 1: no network, no cache, no 15 GB matrix. Every test here runs
against a synthetic manifest built in-process, so a clean clone can run
them immediately after `pip install -r requirements-lock.txt`.

    python -m unittest discover -s tests -v

Why these specific tests: the functions covered here are the ones that
decide WHICH FILE an analysis reads. A silent wrong answer from any of
them produces results computed on the wrong region or the wrong release,
which no downstream assertion would catch. The size check in h5ad_path
is tested because it caught a real truncated download (a dropped S3
connection left a 626 MB stub where 2.05 GB was expected).
"""

import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import seaad  # noqa: E402


def fake_manifest(regions=("MTG-10X", "HIP-10X", "LEC-10X"), version="20260630"):
    """A manifest with the real nesting shape, but no real data.

    Mirrors the structure the live manifest uses: expression matrices nest
    region/normalization/files/ext, metadata nests name/files/ext.
    """
    em = {}
    for r in regions:
        em[r] = {}
        for norm in ("raw", "log2"):
            em[r][norm] = {"files": {"h5ad": {
                "version": version,
                "relative_path":
                    f"expression_matrices/SEA-AD-Multiregion-10X/{version}/"
                    f"{r}-{norm}.h5ad",
                "size": 1000,
                "file_hash": f"hash-{r}-{norm}",
                "url": f"https://example.invalid/{r}-{norm}.h5ad",
            }}}
    return {"file_listing": {
        "SEA-AD-Multiregion-10X": {
            "expression_matrices": em,
            "metadata": {name: {"files": {"csv": {
                "version": "20260801",
                "relative_path": f"metadata/SEA-AD-Multiregion-10X/20260801/{name}.csv",
                "size": 10, "file_hash": f"hash-{name}",
                "url": f"https://example.invalid/{name}.csv",
            }}} for name in ("cell_metadata", "donor", "disease")},
        },
    }}


class TestResolveRegion(unittest.TestCase):
    """Region selection is the main user-facing input; typos must not pass."""

    def setUp(self):
        self.m = fake_manifest()

    def test_exact_label(self):
        self.assertEqual(seaad.resolve_region(self.m, "MTG-10X"), "MTG-10X")

    def test_bare_name_gets_suffix(self):
        self.assertEqual(seaad.resolve_region(self.m, "HIP"), "HIP-10X")

    def test_case_insensitive(self):
        self.assertEqual(seaad.resolve_region(self.m, "lec"), "LEC-10X")
        self.assertEqual(seaad.resolve_region(self.m, "mtg-10x"), "MTG-10X")

    def test_unknown_region_lists_the_valid_ones(self):
        with self.assertRaises(SystemExit) as cm:
            seaad.resolve_region(self.m, "cerebellum")
        msg = str(cm.exception)
        self.assertIn("cerebellum", msg)
        # The error must name every valid option, not just say "invalid".
        for r in ("MTG-10X", "HIP-10X", "LEC-10X"):
            self.assertIn(r, msg)

    def test_full_anatomical_name_resolves(self):
        # 'Hippocampus' is an accepted alias, not an error. This was a
        # behaviour change when anatomical names were added, and the
        # previous version of this test caught it.
        self.assertEqual(seaad.resolve_region(self.m, "Hippocampus"), "HIP-10X")

    def test_available_regions_is_sorted(self):
        self.assertEqual(
            seaad.available_regions(self.m), ["HIP-10X", "LEC-10X", "MTG-10X"])


class TestRegionAnatomy(unittest.TestCase):
    """Codes like DFC are opaque; anatomical names must resolve.

    Nobody guesses that the prefrontal cortex is filed as DFC, so
    `--region PFC` has to work or the data is effectively undiscoverable.
    """

    def setUp(self):
        self.m = fake_manifest(
            regions=("MTG-10X", "DFC-10X", "MEC-10X", "LEC-10X", "ITG-10X",
                     "STG-10X"))

    def test_pfc_aliases_resolve_to_dfc(self):
        for q in ("PFC", "pfc", "DLPFC", "dlpfc", "prefrontal", "dpfc"):
            self.assertEqual(seaad.resolve_region(self.m, q), "DFC-10X", q)

    def test_anatomical_word_resolves(self):
        # Not an alias; matched against the description text.
        self.assertEqual(seaad.resolve_region(self.m, "dorsolateral"), "DFC-10X")

    def test_ambiguous_description_is_refused_not_guessed(self):
        # 'temporal' is ITG, MTG and STG. Picking one silently would
        # analyse the wrong brain region.
        for q in ("temporal", "entorhinal"):
            with self.assertRaises(SystemExit) as cm:
                seaad.resolve_region(self.m, q)
            self.assertIn("more than one", str(cm.exception))

    def test_known_region_absent_from_this_directory(self):
        # Head of caudate is a real SEA-AD region, but ships in a
        # different directory. Say so rather than "unknown region".
        with self.assertRaises(SystemExit) as cm:
            seaad.resolve_region(self.m, "caudate")
        msg = str(cm.exception)
        self.assertIn("CaH-10X", msg)
        self.assertIn("not in", msg)

    def test_errors_include_anatomical_names(self):
        with self.assertRaises(SystemExit) as cm:
            seaad.resolve_region(self.m, "nonsense")
        self.assertIn("dorsolateral prefrontal cortex", str(cm.exception))

    def test_every_multiregion_code_has_a_description(self):
        for code in ("AnG-10X", "DFC-10X", "FI-10X", "HIP-10X", "ITG-10X",
                     "LEC-10X", "MEC-10X", "MTG-10X", "STG-10X", "V1C-10X"):
            self.assertTrue(seaad.region_description(code), code)

    def test_aliases_point_at_real_codes(self):
        for alias, code in seaad.REGION_ALIASES.items():
            self.assertIn(code, seaad.REGION_INFO, f"{alias} -> {code}")


class TestMatrixEntry(unittest.TestCase):

    def setUp(self):
        self.m = fake_manifest()

    def test_returns_the_pinned_record(self):
        e = seaad.matrix_entry(self.m, "MTG-10X", "raw")
        self.assertEqual(e["file_hash"], "hash-MTG-10X-raw")
        self.assertEqual(e["version"], "20260630")

    def test_raw_and_log2_are_different_files(self):
        self.assertNotEqual(
            seaad.matrix_entry(self.m, "MTG-10X", "raw")["relative_path"],
            seaad.matrix_entry(self.m, "MTG-10X", "log2")["relative_path"])

    def test_unknown_normalization_is_refused(self):
        with self.assertRaises(SystemExit):
            seaad.matrix_entry(self.m, "MTG-10X", "cpm")


class TestH5adPath(unittest.TestCase):
    """The file-location contract, including the truncation guard."""

    def setUp(self):
        self.m = fake_manifest()
        self._tmp = tempfile.TemporaryDirectory()
        self._saved = seaad.CACHE_DIR
        seaad.CACHE_DIR = pathlib.Path(self._tmp.name)
        self.entry = seaad.matrix_entry(self.m, "MTG-10X", "raw")
        self.target = seaad.CACHE_DIR / self.entry["relative_path"]
        self.target.parent.mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        seaad.CACHE_DIR = self._saved
        self._tmp.cleanup()

    def test_version_directory_comes_from_the_manifest(self):
        # The point of resolving rather than hard-coding: change the
        # manifest's version and the path must follow it.
        other = fake_manifest(version="20991231")
        p = seaad.h5ad_path(other, "MTG-10X", "raw", require=False)
        self.assertIn("20991231", str(p))
        self.assertNotIn("20260630", str(p))

    def test_missing_file_names_the_fetch_command(self):
        with self.assertRaises(SystemExit) as cm:
            seaad.h5ad_path(self.m, "MTG-10X", "raw")
        msg = str(cm.exception)
        self.assertIn("step0b_fetch.py", msg)
        self.assertIn("MTG-10X/raw", msg)

    def test_truncated_file_is_refused(self):
        # The real failure this guards: an S3 read timeout leaves a short
        # file that h5py would either reject obscurely or partially read.
        self.target.write_bytes(b"x" * 500)  # manifest says 1000
        with self.assertRaises(SystemExit) as cm:
            seaad.h5ad_path(self.m, "MTG-10X", "raw")
        msg = str(cm.exception)
        self.assertIn("500", msg)
        self.assertIn("1,000", msg)
        self.assertIn("--force", msg)

    def test_correct_size_is_accepted(self):
        self.target.write_bytes(b"x" * 1000)
        self.assertEqual(seaad.h5ad_path(self.m, "MTG-10X", "raw"), self.target)

    def test_require_false_skips_both_checks(self):
        p = seaad.h5ad_path(self.m, "MTG-10X", "raw", require=False)
        self.assertFalse(p.exists())


class TestOutputNaming(unittest.TestCase):
    """Outputs must be region-scoped or regions overwrite each other."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self._saved = seaad.RESULTS
        seaad.RESULTS = pathlib.Path(self._tmp.name)

    def tearDown(self):
        seaad.RESULTS = self._saved
        self._tmp.cleanup()

    def test_region_is_in_the_filename(self):
        self.assertEqual(seaad.qc_csv("MTG-10X").name,
                         "step2_qc_metrics_MTG-10X.csv")
        self.assertEqual(seaad.annotation_csv("HIP-10X").name,
                         "step3_annotation_HIP-10X.csv")

    def test_two_regions_never_collide(self):
        self.assertNotEqual(seaad.qc_csv("MTG-10X"), seaad.qc_csv("HIP-10X"))
        self.assertNotEqual(seaad.annotation_csv("MTG-10X"),
                            seaad.annotation_csv("HIP-10X"))

    def test_step2_output_is_step3_input(self):
        # step3 reads exactly what step2 wrote; if these ever diverge the
        # pipeline breaks at a distance.
        self.assertEqual(seaad.qc_csv("LEC-10X"), seaad.qc_csv("LEC-10X"))

    def test_results_directory_is_created(self):
        seaad.RESULTS = pathlib.Path(self._tmp.name) / "nested" / "results"
        p = seaad.qc_csv("MTG-10X")
        self.assertTrue(p.parent.is_dir())


class TestRegionArgs(unittest.TestCase):

    def test_default_region_is_mtg(self):
        import argparse
        ap = seaad.add_region_args(argparse.ArgumentParser())
        self.assertEqual(ap.parse_args([]).region, "MTG-10X")

    def test_normalization_defaults_to_raw_when_offered(self):
        import argparse
        ap = seaad.add_region_args(argparse.ArgumentParser(), normalization=True)
        self.assertEqual(ap.parse_args([]).normalization, "raw")
        self.assertEqual(ap.parse_args(["--normalization", "log2"]).normalization,
                         "log2")


class TestPinIsSingleSourced(unittest.TestCase):
    """The pin drifting between scripts is the failure this guards."""

    def test_steps_do_not_redeclare_the_manifest(self):
        scripts = pathlib.Path(__file__).resolve().parent.parent / "scripts"
        for name in ("step0_discover.py", "step1_inventory.py",
                     "step2_qc_metrics.py", "step3_annotation.py"):
            text = (scripts / name).read_text()
            self.assertNotIn(
                'MANIFEST = "releases', text,
                f"{name} declares its own pin; it must import seaad.MANIFEST")

    def test_steps_do_not_hardcode_a_version_directory(self):
        scripts = pathlib.Path(__file__).resolve().parent.parent / "scripts"
        for name in ("step0_discover.py", "step1_inventory.py",
                     "step2_qc_metrics.py", "step3_annotation.py"):
            text = (scripts / name).read_text()
            self.assertNotIn(
                "/20260630/", text,
                f"{name} hard-codes a version directory; resolve it instead")


if __name__ == "__main__":
    unittest.main(verbosity=2)
