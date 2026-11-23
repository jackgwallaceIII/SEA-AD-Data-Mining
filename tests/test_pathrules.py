"""
test_pathrules.py — unit tests for scripts/pathrules.py (key classification).

No network, no data. Keys are real, copied from read-only listings of
the three open SEA-AD buckets on 2026-11-21, EXCEPT where a test is
about an edge case that the buckets do not contain (an unknown
extension, an unknown bucket, a bucket-specific skip rule); those keys
are obviously synthetic ('thing.weird', 'a/b/c.h5ad', 'temp/x.csv' in
the single-cell bucket). The tests pin how the real layout is
classified, so a rule change that re-files real data shows up here.

    python -m unittest discover -s tests -t tests -v

Why these tests: the classifier decides what a catalog row claims a file
is. The failure modes that matter are (a) filing data under the wrong
region or modality, (b) exploding one dataset into thousands of rows,
(c) silently dropping something, and (d) guessing instead of saying
'unclassified'.
"""

import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import pathrules  # noqa: E402

SC = "sea-ad-single-cell-profiling"
ST = "sea-ad-spatial-transcriptomics"
NP = "sea-ad-quantitative-neuropathology"

TEN = ("HIP", "MEC", "LEC", "ITG", "MTG", "FI", "STG", "DFC", "AnG", "V1C")


def c(bucket, key):
    return pathrules.classify(bucket, key)


class TestSingleCellDatasets(unittest.TestCase):

    def test_region_dataset(self):
        r = c(SC, "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad")
        self.assertEqual(r.regions, ("MTG",))
        self.assertEqual(r.region_from, "folder")
        self.assertEqual((r.modality, r.modality_from), ("snRNA-seq", "path"))
        self.assertEqual((r.product, r.format), ("dataset", "h5ad"))
        self.assertFalse(r.deprecated or r.collapsed or r.per_donor)
        self.assertEqual(r.missing, ())

    def test_metadata_table(self):
        r = c(SC, "MTG/RNAseq/SEAAD_MTG_RNAseq_all-nuclei_metadata.2026-06-22.csv")
        self.assertEqual((r.product, r.format), ("cell metadata table", "csv"))

    def test_pfc_folder_holds_dfc_files(self):
        # The folder says PFC, the file says DFC; both are the same code.
        r = c(SC, "PFC/RNAseq/SEAAD_DFC_RNAseq_final-nuclei.2026-06-22.h5ad")
        self.assertEqual(r.regions, ("DFC",))

    def test_atac(self):
        r = c(SC, "MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad")
        self.assertEqual((r.modality, r.product), ("snATAC-seq", "dataset"))

    def test_atac_bigwig(self):
        r = c(SC, "MTG/ATACseq/bw/ADNC2Micro-PVM.bw")
        self.assertEqual((r.modality, r.product, r.format),
                         ("snATAC-seq", "coverage track", "bigwig"))

    def test_per_donor_object(self):
        r = c(SC, "AnG/RNAseq/donors_objects/"
                  "H19.33.004_SEAAD_AnG_RNAseq_final-nuclei.2026-06-22.h5ad")
        self.assertEqual(r.product, "per-donor object")
        self.assertTrue(r.per_donor)
        # The ATAC folder is spelled donor_objects (no s); both must match.
        r = c(SC, "MTG/ATACseq/donor_objects/"
                  "H19.30.002_SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad")
        self.assertEqual(r.product, "per-donor object")

    def test_rds_reference(self):
        r = c(SC, "MTG/RNAseq/Reference_MTG_RNAseq_all-nuclei.2022-06-07.rds")
        self.assertEqual((r.format, r.product), ("rds", "dataset"))

    def test_external_annotation(self):
        r = c(SC, "PFC/RNAseq/public_datasets/Mathys_2019-cell-annotation.2024-03-27.csv")
        self.assertEqual(r.product, "external dataset annotation")

    def test_pickle_has_a_format(self):
        r = c(SC, "MTG/RNAseq/Supplementary Information/SCENIC Results/"
                  "Excitatory-Neurons/motifs/menr.pkl")
        self.assertEqual(r.format, "pickle")
        self.assertEqual(r.product, "supplementary results")


class TestDeprecated(unittest.TestCase):

    def test_previous_objects_flagged_not_dropped(self):
        r = c(SC, "MTG/RNAseq/previous_objects/"
                  "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
        self.assertTrue(r.deprecated)
        self.assertEqual(r.skip, "")          # still catalogued
        self.assertEqual(r.product, "dataset")

    def test_nested_previous_objects(self):
        r = c(SC, "PFC/RNAseq/previous_objects/donor_objects/"
                  "H21.33.046_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad")
        self.assertTrue(r.deprecated)
        self.assertEqual(r.regions, ("DFC",))   # the folder wins over 'A9'


class TestProjectFolders(unittest.TestCase):

    def test_dream_region_from_filename(self):
        r = c(SC, "DREAM/SEAAD_MTG_RNAseq_DREAM.2025-07-15.h5ad")
        self.assertEqual((r.regions, r.region_from), (("MTG",), "filename"))
        self.assertEqual(r.modality_from, "filename")
        self.assertEqual(r.product, "dataset")

    def test_dream_a9_is_dfc(self):
        # 'A9' (Brodmann area 9) was added as a DFC alias on 2026-11-22,
        # a human decision based on the bucket README; see regions.csv notes.
        r = c(SC, "DREAM/SEAAD_A9_RNAseq_DREAM.2025-07-15.h5ad")
        self.assertEqual((r.regions, r.region_from), (("DFC",), "filename"))
        self.assertEqual(r.missing, ())

    def test_multiregion_file_lists_all_ten_regions(self):
        r = c(SC, "Multiregion_2026/subclass_objects/SEAAD_Immune_HIP_MEC_LEC_"
                  "ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad")
        self.assertEqual(r.regions, TEN)
        self.assertEqual(r.product, "per-subclass object")

    def test_multiregion_model_output_is_project_scoped(self):
        r = c(SC, "Multiregion_2026/model_outputs/gpboost_differential_expression_"
                  "modeling/CPS_Local/Astro_1_DE_results_corrected.csv")
        self.assertEqual(r.regions, ())
        self.assertTrue(r.region_from.startswith("project"))
        self.assertNotIn("region", r.missing)
        self.assertEqual((r.modality, r.modality_from),
                         ("snRNA-seq", "project default"))

    def test_region_subfolder_inside_a_project(self):
        r = c(SC, "Multiregion_2026/model_outputs/subclass_annotation/AnG/"
                  "iterative_scANVI_results.2024-06-12.csv")
        self.assertEqual(r.regions, ("AnG",))

    def test_caudate_modality_is_not_assumed(self):
        # The file name says nothing about RNA vs ATAC.
        r = c(SC, "Caudate_Nucleus/CaH_all_nuclei.2025-12-29.h5ad")
        self.assertEqual(r.regions, ("CaH",))
        self.assertEqual(r.modality, "unclassified")
        self.assertIn("modality", r.missing)


class TestCollapsing(unittest.TestCase):
    """Thousands of keys that are one dataset become one unit."""

    def test_plot_folder(self):
        a = c(SC, "MTG/RNAseq/pseudoprogression-plots/A1BG-AS1/overview.jpg")
        b = c(SC, "MTG/RNAseq/pseudoprogression-plots/AC008840.1/supertype-VLMC.jpg")
        self.assertTrue(a.collapsed)
        self.assertEqual(a.unit, "MTG/RNAseq/pseudoprogression-plots/")
        self.assertEqual(a.unit, b.unit)
        self.assertEqual(a.product, "plots")

    def test_zarr_store_including_parquet_inside_it(self):
        base = "caudate_nucleus/H19.33.004/1417642334/combined_spatialdata_zarr/"
        keys = [base + ".zgroup",
                base + "images/morphology_focus/0/0/0/0",
                base + "points/transcripts/points.parquet/part.0.parquet",
                base + "tables/mapped_table/X/data/0"]
        units = {c(ST, k).unit for k in keys}
        self.assertEqual(units, {base})
        self.assertEqual(c(ST, keys[2]).format, "zarr")   # outermost wins

    def test_xenium_bundle(self):
        k = ("caudate_nucleus/H19.33.004/1417642334/xenium_output_and_alignment/"
             "output-XETG00044__0031423__1417642334__20250319__192442/"
             "cell_feature_matrix.h5")
        r = c(ST, k)
        self.assertTrue(r.collapsed)
        self.assertTrue(r.unit.endswith("__192442/"))
        self.assertEqual(r.product, "instrument output (Xenium)")

    def test_deep_zoom_tiles(self):
        k = ("middle-temporal-gyrus-and-superior-temporal-gyrus/H19.33.004/"
             "H19.33.004-A06-NeuN-analysis/H19.33.004-A06-NeuN_analysis_files/10/0_0.jpeg")
        r = c(NP, k)
        self.assertEqual(r.product, "deep-zoom tiles")
        self.assertTrue(r.unit.endswith("_analysis_files/"))
        self.assertEqual(r.regions, ("MTG", "STG"))

    def test_hive_partitioned_parquet(self):
        a = c(SC, "Multiregion_2026/previous_objects/parquet/Multiregion_2026/"
                  "counts_by_gene/gene=A1BG/data.parquet")
        b = c(SC, "Multiregion_2026/previous_objects/parquet/Multiregion_2026/"
                  "counts_by_gene/gene=A1CF/data.parquet")
        self.assertEqual(a.unit, b.unit)
        self.assertTrue(a.unit.endswith("counts_by_gene/"))
        self.assertTrue(a.deprecated)

    def test_unit_is_always_a_prefix_of_the_key(self):
        for b, k in ((SC, "MTG/RNAseq/pseudoprogression-plots/A1BG-AS1/overview.svg"),
                     (SC, "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad")):
            r = c(b, k)
            self.assertTrue(k.startswith(r.unit))

    def test_ordinary_files_are_not_collapsed(self):
        r = c(NP, "middle-temporal-gyrus/H19.33.004/H19.33.004-A06-NEUN/"
                  "H19.33.004-A06-NEUN.svs")
        self.assertFalse(r.collapsed)
        self.assertEqual(r.unit, r.key)


class TestSpatial(unittest.TestCase):

    def test_platform_is_not_asserted(self):
        # Nothing in this path proves MERFISH or Xenium.
        r = c(ST, "middle-temporal-gyrus/H20.33.001/1194111462/"
                  "cellpose-detected_transcripts.csv")
        self.assertEqual((r.modality, r.modality_from),
                         ("spatial transcriptomics", "bucket default"))
        self.assertEqual(r.product, "transcript table")
        self.assertEqual(r.regions, ("MTG",))

    def test_max_projection_image(self):
        r = c(ST, "middle-temporal-gyrus/H20.33.001/1194111462/DAPI_Max.tif")
        self.assertEqual((r.product, r.format), ("image (max projection)", "tiff"))

    def test_combined_anndata(self):
        r = c(ST, "hippocampus/combined_anndata_object/1444211893_HPF_mapped.h5ad")
        self.assertEqual((r.regions, r.product), (("HIP",), "cell-by-gene table"))


class TestNeuropathology(unittest.TestCase):

    def test_whole_slide_image(self):
        r = c(NP, "dorsolateral-prefrontal-cortex/H19.33.004/H19.33.004-A01-AT/"
                  "H19.33.004-A01-AT.svs")
        self.assertEqual((r.regions, r.product, r.format),
                         (("DFC",), "whole-slide image", "svs"))
        self.assertEqual(r.modality, "neuropathology")

    def test_analysis_mask_and_annotation(self):
        base = ("hippocampus-medial-entorhinal-cortex/H19.33.004/H19.33.004-A12-ASYN/"
                "analysis_masks/H19.33.004-A12-ASYN_CA1-L_analysis")
        for ext, fmt in ((".tif", "tiff"), (".annotations", "halo annotations")):
            r = c(NP, base + ext)
            self.assertEqual((r.product, r.format), ("analysis mask", fmt))
            self.assertEqual(r.regions, ("HIP", "MEC"))

    def test_tabular_region_from_filename(self):
        r = c(NP, "full-tabular-data/CaH_all-values_aibs-ids_20250922.csv")
        self.assertEqual((r.regions, r.region_from), (("CaH",), "filename"))
        self.assertEqual(r.product, "quantification table")

    def test_extrastriate_folder_is_v1c_and_esoc(self):
        r = c(NP, "primary-visual-cortex-extrastriate-occipital-cortex/H19.33.004/"
                  "H19.33.004-A8-ASYN/H19.33.004-A8-ASYN.svs")
        self.assertEqual((r.regions, r.region_from), (("V1C", "ESOC"), "folder"))
        self.assertEqual(r.missing, ())

    def test_v1c_esoc_table(self):
        r = c(NP, "full-tabular-data/V1C-ESOC_all-values_aibs-ids_20250922.csv")
        self.assertEqual((r.regions, r.region_from), (("V1C", "ESOC"), "filename"))


class TestSkipped(unittest.TestCase):
    """Not data — skipped, but always with a reason."""

    def test_viewer_and_scratch(self):
        self.assertEqual(
            c(NP, "simple_viewer/openseadragon/images/button_grouphover.png").skip,
            "web viewer assets")
        self.assertEqual(c(NP, "temp/images/ABAT.svg").skip, "scratch folder")

    def test_checkpoint_copies(self):
        r = c(NP, "temp/images/.ipynb_checkpoints/ABAT-checkpoint.svg")
        self.assertTrue(r.skip)

    def test_folder_marker(self):
        self.assertEqual(c(SC, "MTG/RNAseq/").skip, "folder marker")

    def test_skip_rules_are_bucket_specific(self):
        # 'temp' is scratch in the neuropathology bucket only.
        self.assertEqual(c(SC, "temp/x.csv").skip, "")

    def test_skipped_records_report_nothing_missing(self):
        self.assertEqual(
            c(NP, "simple_viewer/dorsolateral-prefrontal-cortex/H19.33.004.html").missing,
            ())


class TestNeuropathologyFullListing(unittest.TestCase):
    """Patterns first seen in the full neuropathology refresh, 2026-11-22."""

    MS = ("middle-temporal-gyrus-and-superior-temporal-gyrus/H19.33.004/"
          "H19.33.004-A06-NeuN-analysis/")

    def test_layer_masks_and_annotations(self):
        for k in (self.MS + "subregion_annotations/H19.33.004-A06-NeuN_Layer1_analysis.annotations",
                  self.MS + "halo_masked_subregions/H19.33.004-A06-NeuN_Layer1_analysis.tif"):
            r = c(NP, k)
            self.assertEqual(r.product, "analysis mask", k)
            self.assertEqual(r.regions, ("MTG", "STG"))
            self.assertEqual(r.missing, ())

    def test_annotation_overlay(self):
        self.assertEqual(c(NP, self.MS + "annotation.svg").product, "annotation overlay")

    def test_per_slide_table(self):
        r = c(NP, "dorsolateral-prefrontal-cortex/H19.33.004/H19.33.004-A01-AT/"
                  "downsampled/H19.33.004-A01-AT.csv")
        self.assertEqual(r.product, "per-slide table")

    def test_slide_annotation(self):
        r = c(NP, "middle-temporal-gyrus/H19.33.004/H19.33.004-A06-NEUN/"
                  "H19.33.004-A06-NEUN.annotations")
        self.assertEqual(r.product, "slide annotation")

    def test_per_donor_json(self):
        r = c(NP, "middle-temporal-gyrus-and-superior-temporal-gyrus/H19.33.004/"
                  "H19.33.004.json")
        self.assertEqual(r.product, "per-donor json")

    def test_hyphenated_region_pair_in_filename(self):
        r = c(NP, "full-tabular-data/Hip-MEC_all-values_aibs-ids_20250922.csv")
        self.assertEqual((r.regions, r.region_from), (("HIP", "MEC"), "filename"))

    def test_hyphen_rule_needs_every_part_to_be_a_region(self):
        # 'Visual' alone is a V1C alias; 'Higher-Order-Visual' must not be.
        # SYNTHETIC key: the real file sits in a region folder, which would
        # decide the region before the filename is ever read, so a
        # non-region folder is needed to exercise the filename rule.
        r = c(NP, "full-tabular-data/H19.33.004-A8-ASYN_Higher-Order-Visual_analysis.csv")
        self.assertEqual(r.regions, ())
        self.assertIn("region", r.missing)

    def test_os_junk_is_skipped(self):
        for k in ("hippocampus-medial-entorhinal-cortex/H20.33.001/H20.33.001-A12-GFAP/"
                  "analysis_masks/Thumbs.db",
                  "middle-temporal-gyrus-and-superior-temporal-gyrus/H20.33.038/"
                  "H20.33.038-A6-NeuN-analysis/.DS_Store"):
            self.assertEqual(c(NP, k).skip, "operating-system junk file", k)


class TestSpatialPlatform(unittest.TestCase):

    def test_merfish_named_in_file(self):
        r = c(ST, "middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad")
        self.assertEqual((r.modality, r.modality_from),
                         ("spatial transcriptomics (MERFISH)", "filename"))
        self.assertEqual(r.product, "cell-by-gene table")
        self.assertEqual(r.missing, ())

    def test_xenium_named_in_folder(self):
        r = c(ST, "caudate_nucleus/H19.33.004/1417642334/xenium_output_and_alignment/"
                  "output-XETG00044__0031423__1417642334__20250319__192442/"
                  "cell_feature_matrix.h5")
        self.assertEqual((r.modality, r.modality_from),
                         ("spatial transcriptomics (Xenium)", "path"))

    def test_unnamed_platform_stays_generic(self):
        r = c(ST, "middle-temporal-gyrus/H20.33.001/1194111462/DAPI_Max.tif")
        self.assertEqual(r.modality, "spatial transcriptomics")


class TestTestData(unittest.TestCase):
    """A mouse dev table must not inherit 'human multi-region snRNA-seq'."""

    UNIT = ("Multiregion_2026/previous_objects/parquet/TEST_dev_mouse/"
            "counts_by_gene/")

    def test_collapsed_test_unit_gets_no_project_defaults(self):
        r = pathrules.classify_unit(SC, self.UNIT, "partitioned parquet",
                                    "partitioned table")
        self.assertEqual(r.product, "test/development data")
        self.assertEqual((r.regions, r.region_from), ((), ""))
        self.assertEqual(r.modality, "")
        self.assertEqual(r.missing, ())
        self.assertTrue(r.deprecated)

    def test_real_project_unit_still_gets_defaults(self):
        r = pathrules.classify_unit(
            SC, "Multiregion_2026/previous_objects/parquet/Multiregion_2026/"
                "counts_by_gene/", "partitioned parquet", "partitioned table")
        self.assertEqual(r.modality_from, "project default")
        self.assertTrue(r.region_from.startswith("project"))


class TestUnclassifiedIsHonest(unittest.TestCase):

    def test_unknown_extension(self):
        r = c(SC, "MTG/RNAseq/thing.weird")
        self.assertEqual(r.format, "unknown:weird")
        self.assertIn("format", r.missing)

    def test_unknown_product(self):
        r = c(SC, "DREAM/test.csv")
        self.assertEqual(r.product, "unclassified")
        self.assertIn("product", r.missing)

    def test_documentation_needs_no_region(self):
        for k in ("README.md", "MTG/RNAseq/Changelog.html", "index.html"):
            r = c(SC, k)
            self.assertEqual(r.product, "documentation", k)
            self.assertEqual(r.missing, (), k)

    def test_every_record_has_every_field(self):
        for b, k in ((SC, "a/b/c.h5ad"), (ST, "x"), (NP, "y/z.svs"),
                     ("some-other-bucket", "q/r.csv")):
            r = c(b, k)
            self.assertTrue(r.format)
            self.assertTrue(r.product)
            self.assertTrue(r.modality)


if __name__ == "__main__":
    unittest.main(verbosity=2)
