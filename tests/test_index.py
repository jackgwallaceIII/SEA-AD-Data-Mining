"""
test_index.py — unit tests for scripts/index.py (catalog refresh).

No network: the walker runs over test_sources.FakeS3, an in-memory
bucket with real S3 paging behaviour. Every test writes into its own
temporary catalog directory, never into the real catalog/.

    python -m unittest discover -s tests -t tests -v

Why these tests: the catalog is what people will read numbers from. The
failure modes that matter are summary totals that disagree with the
unit table, double counting (multi-region units, deprecated data),
skipped folders reported as zero, a capped count shown as exact,
refreshing one bucket wiping another's rows, and a refresh diff that
misses a change.
"""

import csv
import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import index  # noqa: E402
from test_sources import FakeS3, MTG_RNA, NEURO, PARQUET, PLOTS  # noqa: E402

SC = "sea-ad-single-cell-profiling"
NP = "sea-ad-quantitative-neuropathology"

SC_KEYS = {
    **MTG_RNA,
    **PARQUET,
    "README.md": 6,
    "DREAM/SEAAD_MTG_RNAseq_DREAM.2025-07-15.h5ad": 300,
    "Multiregion_2026/subclass_objects/SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_"
    "STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad": 400,
    "Multiregion_2026/model_outputs/gpboost_differential_expression_modeling/"
    "CPS_Local/Astro_1_DE_results_corrected.csv": 50,
    "Caudate_Nucleus/CaH_all_nuclei.2025-12-29.h5ad": 2000,
}
NP_KEYS = {
    **NEURO,
    "middle-temporal-gyrus-and-superior-temporal-gyrus/H19.33.004/"
    "H19.33.004-A06-NeuN-primary/H19.33.004-A06-NeuN-primary.svs": 700,
}


class CatalogCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = pathlib.Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def refresh(self, bucket, keys, **kw):
        return index.refresh(bucket, client=FakeS3(keys), catalog_dir=self.dir,
                             out=lambda *a: None, **kw)

    def summary(self):
        with open(self.dir / "catalog.csv", newline="") as fh:
            return list(csv.DictReader(fh))

    def row(self, res, region, modality=None):
        hits = [s for s in res["summary"] if s["region"] == region
                and (modality is None or s["modality"] == modality)]
        self.assertEqual(len(hits), 1, f"{region} {modality}: {hits}")
        return hits[0]


class TestTotalsAgree(CatalogCase):

    def test_summary_equals_sum_of_units(self):
        res = self.refresh(SC, SC_KEYS)
        units = [r for r in res["rows"] if r["kind"] != "skipped"]
        cur = sum(r["size"] for r in units if not r["deprecated"])
        dep = sum(r["size"] for r in units if r["deprecated"])
        self.assertEqual(sum(s["total_size"] for s in res["summary"]), cur)
        self.assertEqual(sum(s["deprecated_size"] for s in res["summary"]), dep)

    def test_every_byte_in_the_bucket_is_accounted_for(self):
        res = self.refresh(SC, SC_KEYS)
        self.assertEqual(res["current"] + res["deprecated"], sum(SC_KEYS.values()))

    def test_no_duplicate_units(self):
        res = self.refresh(SC, SC_KEYS)
        units = [r["unit"] for r in res["rows"]]
        self.assertEqual(len(units), len(set(units)))


class TestGrouping(CatalogCase):

    def test_multi_region_unit_is_one_row_not_ten(self):
        res = self.refresh(SC, SC_KEYS)
        ten = "HIP;MEC;LEC;ITG;MTG;FI;STG;DFC;AnG;V1C"
        self.assertEqual(self.row(res, ten)["total_size"], 400)
        # And it is NOT added into the MTG row.
        mtg = self.row(res, "MTG", "snRNA-seq")
        self.assertEqual(mtg["total_size"], 1 + 1000 + 100 + 150 + 300)

    def test_project_scoped_outputs_get_their_own_row(self):
        res = self.refresh(SC, SC_KEYS)
        r = self.row(res, "multi-region (project)")
        self.assertEqual(r["total_size"], 50)

    def test_two_region_folder_is_one_row(self):
        res = self.refresh(NP, NP_KEYS)
        self.assertEqual(self.row(res, "MTG;STG")["total_size"], 700)

    def test_region_aliases_and_names_filled(self):
        res = self.refresh(SC, SC_KEYS)
        r = self.row(res, "MTG", "snRNA-seq")
        self.assertEqual(r["region_name"], "middle temporal gyrus")
        self.assertIn("DHBA:12141", r["region_aliases"])
        self.assertEqual(r["access"], "open")

    def test_unclassified_modality_is_visible_and_counted_for_review(self):
        res = self.refresh(SC, SC_KEYS)
        r = self.row(res, "CaH", "unclassified")
        self.assertEqual((r["total_size"], r["needs_review"]), (2000, 1))


class TestDeprecatedAndSkipped(CatalogCase):

    def test_deprecated_not_in_current_totals(self):
        res = self.refresh(SC, SC_KEYS)
        mtg = self.row(res, "MTG", "snRNA-seq")
        self.assertEqual((mtg["deprecated_units"], mtg["deprecated_size"]), (1, 900))

    def test_skipped_folders_have_unknown_size_and_no_summary_row(self):
        res = self.refresh(NP, NP_KEYS)
        skipped = [r for r in res["rows"] if r["kind"] == "skipped"]
        self.assertEqual({r["unit"] for r in skipped}, {"simple_viewer/", "temp/"})
        self.assertTrue(all(r["size"] == "" for r in skipped))
        self.assertNotIn("(none)", {s["region"] for s in res["summary"]
                                    if s["total_size"] == 0})

    def test_cap_marks_the_summary_row_as_minimum(self):
        res = self.refresh(SC, {**PLOTS, "MTG/RNAseq/x_final-nuclei.h5ad": 5},
                           max_keys_per_unit=2)
        self.assertTrue(self.row(res, "MTG", "snRNA-seq")["size_is_minimum"])
        self.assertTrue(res["minimum"])

    def test_no_cap_means_exact(self):
        res = self.refresh(SC, SC_KEYS)
        self.assertFalse(res["minimum"])


class TestFilesOnDisk(CatalogCase):

    def test_files_written_and_readable(self):
        self.refresh(SC, SC_KEYS)
        units = index.read_units(index.units_path(SC, self.dir))
        self.assertEqual(set(units[0]), set(index.UNIT_FIELDS))
        self.assertTrue(self.summary())
        self.assertTrue((self.dir / "refresh_log.csv").exists())

    def test_dry_run_writes_nothing(self):
        self.refresh(SC, SC_KEYS, dry_run=True)
        self.assertEqual(list(self.dir.iterdir()), [])

    def test_refreshing_one_bucket_keeps_the_others_rows(self):
        self.refresh(SC, SC_KEYS)
        self.refresh(NP, NP_KEYS)
        buckets = {r["bucket"] for r in self.summary()}
        self.assertEqual(buckets, {SC, NP})
        self.refresh(SC, SC_KEYS)                      # again
        self.assertEqual({r["bucket"] for r in self.summary()}, {SC, NP})

    def test_identical_refresh_writes_identical_bytes(self):
        # No spurious git diff when nothing changed.
        self.refresh(SC, SC_KEYS)
        a = index.units_path(SC, self.dir).read_bytes()
        self.refresh(SC, SC_KEYS)
        self.assertEqual(a, index.units_path(SC, self.dir).read_bytes())

    def test_log_appends_one_line_per_refresh(self):
        self.refresh(SC, SC_KEYS)
        self.refresh(SC, SC_KEYS)
        with open(self.dir / "refresh_log.csv") as fh:
            self.assertEqual(len(fh.read().strip().splitlines()), 3)  # header + 2

    def test_unknown_bucket_refused(self):
        with self.assertRaises(SystemExit):
            self.refresh("some-other-bucket", SC_KEYS)


class TestReclassify(CatalogCase):
    """Re-apply rules to the saved catalog, with no S3 access at all."""

    def test_makes_no_s3_calls_and_preserves_sizes(self):
        first = self.refresh(SC, SC_KEYS)
        res = index.reclassify(SC, catalog_dir=self.dir)      # no client given
        self.assertEqual(res["current"], first["current"])
        self.assertEqual(res["deprecated"], first["deprecated"])
        self.assertEqual(res["relabelled"], 0)

    def test_picks_up_a_rule_change(self):
        self.refresh(SC, SC_KEYS)
        # Simulate a catalog written under older rules: blank a label.
        upath = index.units_path(SC, self.dir)
        rows = index.read_units(upath)
        for r in rows:
            if r["unit"].startswith("DREAM/"):
                r["regions"], r["missing"] = "", "region"
        index.write_units(upath, rows)
        res = index.reclassify(SC, catalog_dir=self.dir)
        self.assertEqual(res["relabelled"], 1)
        mtg = [s for s in res["summary"] if s["region"] == "MTG"
               and s["modality"] == "snRNA-seq"][0]
        self.assertEqual(mtg["total_size"], 1 + 1000 + 100 + 150 + 300)

    def test_collapsed_units_keep_their_product(self):
        self.refresh(SC, SC_KEYS)
        res = index.reclassify(SC, catalog_dir=self.dir)
        plots = [r for r in res["rows"] if r["unit"].endswith("pseudoprogression-plots/")]
        self.assertEqual(plots[0]["product"], "plots")

    def test_log_records_the_mode(self):
        self.refresh(SC, SC_KEYS)
        index.reclassify(SC, catalog_dir=self.dir)
        with open(self.dir / "refresh_log.csv", newline="") as fh:
            modes = [r["mode"] for r in csv.DictReader(fh)]
        self.assertEqual(modes, ["refresh", "reclassify"])

    def test_old_log_without_mode_column_is_upgraded(self):
        (self.dir / "refresh_log.csv").write_text(
            "bucket,started,finished,seconds,units,current_size,deprecated_size,"
            "skipped_folders,size_is_minimum,added,removed,changed\n"
            f"{SC},a,b,1,1,1,0,0,False,,,\n")
        self.refresh(SC, SC_KEYS)
        with open(self.dir / "refresh_log.csv", newline="") as fh:
            modes = [r["mode"] for r in csv.DictReader(fh)]
        self.assertEqual(modes, ["refresh", "refresh"])

    def test_requires_an_existing_catalog(self):
        with self.assertRaises(SystemExit):
            index.reclassify(SC, catalog_dir=self.dir)


class TestDiff(CatalogCase):

    def test_first_refresh_has_no_diff(self):
        self.assertIsNone(self.refresh(SC, SC_KEYS)["diff"])

    def test_added_removed_changed(self):
        self.refresh(SC, SC_KEYS)
        keys = dict(SC_KEYS)
        del keys["DREAM/SEAAD_MTG_RNAseq_DREAM.2025-07-15.h5ad"]           # removed
        keys["MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad"] = 1234  # changed
        keys["MTG/RNAseq/pseudoprogression-plots/NEW/overview.jpg"] = 1    # unit changed
        keys["STG/RNAseq/SEAAD_STG_RNAseq_final-nuclei.2026-06-22.h5ad"] = 9   # added
        added, removed, changed = self.refresh(SC, keys)["diff"]
        self.assertEqual(added, ["STG/RNAseq/SEAAD_STG_RNAseq_final-nuclei.2026-06-22.h5ad"])
        self.assertEqual(removed, ["DREAM/SEAAD_MTG_RNAseq_DREAM.2025-07-15.h5ad"])
        self.assertEqual(set(changed), {
            "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad",
            "MTG/RNAseq/pseudoprogression-plots/"})

    def test_unchanged_refresh_reports_nothing(self):
        self.refresh(SC, SC_KEYS)
        self.assertEqual(self.refresh(SC, SC_KEYS)["diff"], ([], [], []))


if __name__ == "__main__":
    unittest.main(verbosity=2)
