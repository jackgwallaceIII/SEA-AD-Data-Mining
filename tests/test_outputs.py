"""
test_outputs.py — regression tests on the pipeline's actual outputs.

Tier 3: reads results/*.csv. Every test is SKIPPED when the file is not
there, so a clean clone still passes; they only run once the pipeline
has produced output. Nothing here downloads or recomputes anything.

    python -m unittest discover -s tests -v
    python -m unittest tests.test_outputs -v          # just these

Why these tests: the facts asserted here were established by inspecting
the data and are the ones a refactor could silently break. In particular
the row count, the completeness of the join, and the exact agreement
between Allen's umi_count and step2's independently computed
total_counts — that last one is the only check that would catch row
misalignment between the h5ad and the metadata CSV, which would corrupt
every result downstream while looking entirely normal.

These are regression tests, not correctness proofs: they pin the current
behaviour so a change that alters it has to be deliberate. If a number
here changes because the pinned release changed, update it and say so in
the commit.
"""

import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import seaad  # noqa: E402

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    pd = None

REGION = "MTG-10X"
QC = seaad.qc_csv(REGION)
ANN = seaad.annotation_csv(REGION)

# Established facts for the MTG-10X matrix under releases/20261101.
N_NUCLEI = 1_178_694
N_METADATA_ROWS = 1_861_334      # before the in-matrix restriction
N_DROPPED = 682_640              # metadata rows with no expression data
N_DONORS = 84
N_IMMUNE = 36_686                # the microglia + PVM subclass
N_COLUMNS = 27

ANN_COLUMNS = {
    "cell_label", "doublet_score", "umi_count", "donor_label", "cluster_alias",
    "class", "subclass", "supertype", "donor_sex", "donor_age_value", "PMI",
    "Brain pH", "Overall AD neuropathological Change", "Thal", "Braak",
    "CERAD score", "Cognitive Status", "APOE Genotype", "CPS_Global",
    "Severely Affected Donor", "LATE", "total_counts", "n_genes_detected",
    "mito_counts", "ribo_counts", "pct_mito", "pct_ribo",
}

CHUNK = 400_000


def chunks(path, usecols):
    return pd.read_csv(path, usecols=usecols, chunksize=CHUNK)


@unittest.skipIf(pd is None, "pandas not installed")
@unittest.skipUnless(QC.exists(), f"step2 output not present ({QC.name})")
class TestQcMetrics(unittest.TestCase):
    """step2: per-nucleus QC, computed by streaming the CSR matrix."""

    def test_one_row_per_nucleus(self):
        n = sum(len(c) for c in chunks(QC, ["cell_label"]))
        self.assertEqual(n, N_NUCLEI)

    def test_cell_label_is_unique(self):
        seen = set()
        for c in chunks(QC, ["cell_label"]):
            seen.update(c["cell_label"])
        self.assertEqual(len(seen), N_NUCLEI)

    def test_no_all_zero_nuclei(self):
        for c in chunks(QC, ["total_counts"]):
            self.assertEqual((c["total_counts"] == 0).sum(), 0)

    def test_allen_prefilter_is_visible(self):
        # Allen filtered at <5% mito and >=500 genes before release. If
        # these bounds are ever violated the matrix is not what we think.
        for c in chunks(QC, ["pct_mito", "n_genes_detected"]):
            self.assertLess(c["pct_mito"].max(), 5.0)
            self.assertGreaterEqual(c["n_genes_detected"].min(), 500)

    def test_percentages_are_in_range(self):
        for c in chunks(QC, ["pct_mito", "pct_ribo"]):
            self.assertGreaterEqual(c["pct_mito"].min(), 0.0)
            self.assertGreaterEqual(c["pct_ribo"].min(), 0.0)
            self.assertLessEqual(c["pct_ribo"].max(), 100.0)


@unittest.skipIf(pd is None, "pandas not installed")
@unittest.skipUnless(ANN.exists(), f"step3 output not present ({ANN.name})")
class TestAnnotation(unittest.TestCase):
    """step3: the joined per-nucleus table."""

    def test_restricted_to_nuclei_in_the_matrix(self):
        # The bug this pins: the table once held 1,861,334 rows because
        # the assertions were baselined on cell_metadata rather than on
        # the matrix, so 682,640 annotation-free rows passed unnoticed.
        n = sum(len(c) for c in chunks(ANN, ["cell_label"]))
        self.assertEqual(n, N_NUCLEI)
        self.assertEqual(N_METADATA_ROWS - N_NUCLEI, N_DROPPED)

    def test_column_set(self):
        head = pd.read_csv(ANN, nrows=1)
        self.assertEqual(set(head.columns), ANN_COLUMNS)
        self.assertEqual(len(head.columns), N_COLUMNS)

    def test_no_missing_values_anywhere(self):
        # After the restriction every row must be fully annotated; a NaN
        # here means a join silently failed.
        for c in pd.read_csv(ANN, chunksize=CHUNK):
            na = c.isna().sum()
            self.assertEqual(na[na > 0].to_dict(), {})

    def test_umi_count_matches_independently_computed_total(self):
        # The row-alignment check. Allen's umi_count and step2's
        # total_counts are computed independently; exact agreement is the
        # only evidence that h5ad row i corresponds to metadata row i.
        worst = 0
        for c in chunks(ANN, ["umi_count", "total_counts"]):
            worst = max(worst, (c["umi_count"] - c["total_counts"]).abs().max())
        self.assertEqual(worst, 0)

    def test_donor_count(self):
        donors = set()
        for c in chunks(ANN, ["donor_label"]):
            donors.update(c["donor_label"])
        self.assertEqual(len(donors), N_DONORS)

    def test_immune_subclass_population(self):
        # The project's cell type of interest; a change here changes
        # every microglial result.
        n = sum((c["subclass"] == "Immune").sum()
                for c in chunks(ANN, ["subclass"]))
        self.assertEqual(n, N_IMMUNE)

    def test_subclass_counts_sum_to_the_total(self):
        counts = {}
        for c in chunks(ANN, ["subclass"]):
            for k, v in c["subclass"].value_counts().items():
                counts[k] = counts.get(k, 0) + v
        self.assertEqual(sum(counts.values()), N_NUCLEI)
        self.assertEqual(len(counts), 24)

    def test_cps_global_is_a_fraction(self):
        for c in chunks(ANN, ["CPS_Global"]):
            self.assertGreaterEqual(c["CPS_Global"].min(), 0.0)
            self.assertLessEqual(c["CPS_Global"].max(), 1.0)


@unittest.skipIf(pd is None, "pandas not installed")
@unittest.skipUnless(QC.exists() and ANN.exists(), "pipeline outputs not present")
class TestStepsAgree(unittest.TestCase):
    """step3 must carry step2's numbers through unchanged."""

    def test_same_nucleus_count(self):
        a = sum(len(c) for c in chunks(QC, ["cell_label"]))
        b = sum(len(c) for c in chunks(ANN, ["cell_label"]))
        self.assertEqual(a, b)

    def test_qc_columns_survive_the_join(self):
        qc_cols = set(pd.read_csv(QC, nrows=1).columns)
        ann_cols = set(pd.read_csv(ANN, nrows=1).columns)
        self.assertTrue(qc_cols.issubset(ann_cols),
                        f"lost in the join: {qc_cols - ann_cols}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
