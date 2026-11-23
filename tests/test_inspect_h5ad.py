"""
test_inspect_h5ad.py — unit tests for the HDF5 / h5ad part of inspect_data.py.

No network, no real data: tiny AnnData files are written into a temporary
folder with the anndata library, so the on-disk layout is the real one.
Values that must never appear in a report are planted on purpose (a donor
ID, cell labels, a gene symbol, an expression value, a uns value) and the
tests assert they are absent from both the .md and the .json.

    python -m unittest discover -s tests -t tests -v

Skipped automatically if anndata, h5py or scipy is not installed.

Why these tests: an h5ad report has to answer "which matrix holds raw
counts?" correctly (it decides every downstream analysis), describe obs
and var with the same withholding rules as a CSV, and never read the
values of uns or the matrices into the report.
"""

import json
import pathlib
import sys
import tempfile
import unittest
import warnings

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import inspect_data as ins  # noqa: E402

try:
    import anndata as ad
    import h5py
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    HAVE = True
except ImportError:                                  # pragma: no cover
    HAVE = False

SC = "sea-ad-single-cell-profiling"
SECRET_DONOR = "H21.33.999"
SECRET_VALUE = 31415.9265
SECRETS = (SECRET_DONOR, "SECRETCELL", "SECRETGENE", "SECRETUNS", "31415",
           "1234.5", "991.0")
N_OBS, N_VAR = 60, 40


def write(adata, path):
    # anndata warns about pandas' future string handling while writing;
    # that is about the test fixture, not about the code under test.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata.write_h5ad(path)


def make_adata():
    rng = np.random.default_rng(0)
    counts = rng.poisson(2.0, size=(N_OBS, N_VAR)).astype(np.float32)
    logged = np.log1p(counts)
    logged[0, 0] = SECRET_VALUE                     # must never be printed
    obs = pd.DataFrame({
        # 6 donors: above the per-donor floor, so levels may show
        "Donor ID": pd.Categorical([SECRET_DONOR if i % 6 == 1 else f"H20.33.{i % 6:03d}"
                                    for i in range(N_OBS)]),
        "Fresh Brain Weight": pd.Categorical(["1234.5", "991.0"] * (N_OBS // 2)),
        "specify other race": pd.Categorical(["nan"] * N_OBS),
        "Secondary Study Name": pd.Categorical(["nan", "SLS"] * (N_OBS // 2)),
        "anything": pd.Categorical(["ACGTACGTACGTACGT-1", "TTGGCCAATTGGCCAA-1"] * (N_OBS // 2)),
        "Braak": pd.Categorical(["Braak II", "Braak IV", "Braak VI"] * (N_OBS // 3)),
        "Severely Affected Donor": pd.Categorical(["N", "Y"] * (N_OBS // 2)),
        "Supertype": pd.Categorical([f"type_{i % 40}" for i in range(N_OBS)]),
        "PMI": [np.nan if i < 5 else float(4 + i % 5) for i in range(N_OBS)],
        "Used in analysis": [True] * N_OBS,
        "Genes detected": np.arange(N_OBS, dtype=np.int64),
        "maybe": pd.array([None if i < 3 else i % 4 for i in range(N_OBS)], dtype="Int64"),
    }, index=[f"SECRETCELL{i:04d}" for i in range(N_OBS)])
    var = pd.DataFrame({"gene_symbol": [f"SECRETGENE{i}" for i in range(N_VAR)]},
                       index=[f"ENSG{i:011d}" for i in range(N_VAR)])
    a = ad.AnnData(X=sp.csr_matrix(logged), obs=obs, var=var)
    a.layers["UMIs"] = sp.csr_matrix(counts)
    a.obsm["X_umap"] = rng.normal(size=(N_OBS, 2)).astype(np.float32)
    a.obsp["connectivities"] = sp.random(N_OBS, N_OBS, density=0.05, format="csr",
                                         random_state=0)
    a.uns["Braak_colors"] = np.array(["#000000", "#111111", "#222222"])
    a.uns["params"] = {"method": "SECRETUNS-value"}
    a.raw = ad.AnnData(X=sp.csr_matrix(counts), var=var)
    return a


@unittest.skipUnless(HAVE, "anndata / h5py / scipy not installed")
class H5adCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(cls.tmp.name)
        cls.data, cls.out = base / "data_cache", base / "results"
        cls.path = cls.data / SC / "LEC/RNAseq/donors_objects/test.h5ad"
        cls.path.parent.mkdir(parents=True)
        write(make_adata(), cls.path)
        cls.r, cls.md = ins.inspect(cls.path, data_dir=cls.data, out_dir=cls.out)
        cls.a = cls.r["anndata"]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def matrix(self, name):
        return next(m for m in self.a["matrices"] if m["name"] == name)

    def col(self, frame, name):
        return next(c for c in self.a[frame]["columns"] if c["name"] == name)


class TestShapeAndMatrices(H5adCase):

    def test_recognised_as_anndata(self):
        self.assertEqual(self.r["detected_format"], "hdf5")
        self.assertEqual(self.r["status"], "inspected (AnnData)")
        self.assertEqual((self.a["n_obs"], self.a["n_vars"]), (N_OBS, N_VAR))

    def test_x_is_sparse_and_judged_transformed(self):
        x = self.matrix("X")
        self.assertEqual((x["storage"], x["shape"]), ("csr sparse", (N_OBS, N_VAR)))
        self.assertIn("non-whole", x["values"])

    def test_counts_layer_judged_raw_counts(self):
        self.assertIn("raw counts", self.matrix("layers/UMIs")["values"])
        self.assertIn("raw counts", self.matrix("raw/X")["values"])

    def test_embeddings_and_graphs_listed_without_values(self):
        names = {m["name"]: m for m in self.a["obsm"] + self.a["obsp"]}
        self.assertEqual(names["obsm/X_umap"]["shape"], (N_OBS, 2))
        self.assertEqual(names["obsp/connectivities"]["storage"], "csr sparse")
        self.assertNotIn("values", names["obsm/X_umap"])


class TestObsVar(H5adCase):

    def test_categories_follow_the_csv_rules(self):
        self.assertEqual(self.col("obs", "Braak")["levels"],
                         ["Braak II", "Braak IV", "Braak VI"])
        self.assertEqual(self.col("obs", "Braak")["type"], "category")
        self.assertEqual(self.col("obs", "Severely Affected Donor")["levels"], ["N", "Y"])
        sup = self.col("obs", "Supertype")
        self.assertEqual((sup["distinct"], sup["levels"]), (40, None))

    def test_identifiers_withheld_by_name_and_by_value(self):
        self.assertIsNotNone(self.col("obs", "Donor ID")["levels_withheld"])
        self.assertIn("donor IDs or barcodes", self.col("obs", "anything")["levels_withheld"])

    def test_numeric_missing_and_nullable(self):
        self.assertEqual(self.col("obs", "PMI")["missing"], 5)
        self.assertIsNone(self.col("obs", "PMI")["levels"])
        self.assertEqual(self.col("obs", "maybe")["missing"], 3)
        self.assertEqual(self.col("obs", "Genes detected")["type"], "integer")
        self.assertEqual(self.col("obs", "Used in analysis")["levels"], ["True"])

    def test_indexes_described_not_listed(self):
        self.assertEqual(self.a["obs"]["index"]["entries"], N_OBS)
        self.assertTrue(self.a["obs"]["index"]["unique"])
        self.assertEqual(self.a["var"]["index"]["looks_like"], "Ensembl gene IDs")
        self.assertIsNotNone(self.col("var", "gene_symbol")["levels_withheld"])


class TestQuirks(H5adCase):

    def test_numbers_stored_as_text_are_not_listed(self):
        c = self.col("obs", "Fresh Brain Weight")
        self.assertEqual(c["type"], "category (numbers as text)")
        self.assertIsNone(c["levels"])

    def test_nan_text_counts_as_missing(self):
        c = self.col("obs", "specify other race")
        self.assertEqual((c["missing"], c["empty"], c["type"]), (N_OBS, True, "empty"))
        c = self.col("obs", "Secondary Study Name")
        self.assertEqual((c["missing"], c["levels"]), (N_OBS // 2, ["SLS"]))
        self.assertIn("'nan'", c["note"])
        self.assertIn("Missing values stored as the text 'nan'", self.md)

    def test_many_donors_no_floor(self):
        self.assertIsNone(self.a["donor_floor"])


@unittest.skipUnless(HAVE, "anndata / h5py / scipy not installed")
class TestPerDonorObject(unittest.TestCase):

    def test_single_donor_object_hides_personal_levels(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            obs = pd.DataFrame({
                "Donor ID": pd.Categorical([SECRET_DONOR] * 12),
                "APOE Genotype": pd.Categorical(["SECRET-APOE"] * 12),
                "Braak": pd.Categorical(["SECRET-BRAAK"] * 12),
                "Subclass": pd.Categorical(["Astro", "Oligo", "Vip"] * 4),
                # the real ATAC case: a category column holding only the word 'nan'
                "Supertype (non-expanded)": pd.Categorical(["nan"] * 12),
            }, index=[f"c{i}" for i in range(12)])
            p = base / "H20.33.036_SEAAD_LEC_RNAseq_final-nuclei.h5ad"
            write(ad.AnnData(X=sp.csr_matrix(np.ones((12, 3), np.float32)), obs=obs), p)
            r, md = ins.inspect(p, data_dir=base, out_dir=base / "out")
        a = r["anndata"]
        self.assertEqual(a["donor_floor"]["donors"], 1)
        levels = {c["name"]: c["levels"] for c in a["obs"]["columns"]}
        self.assertIsNone(levels["APOE Genotype"])
        self.assertIsNone(levels["Braak"])
        self.assertEqual(levels["Subclass"], ["Astro", "Oligo", "Vip"])
        blank = next(c for c in a["obs"]["columns"] if c["name"] == "Supertype (non-expanded)")
        self.assertEqual((blank["type"], blank["levels_withheld"]), ("empty", None))
        self.assertNotIn("SECRET", md + json.dumps(r))


class TestNothingLeaks(H5adCase):

    def test_planted_values_appear_in_neither_output(self):
        blob = self.md + json.dumps(self.r)
        for secret in SECRETS:
            self.assertNotIn(secret, blob, secret)

    def test_uns_key_names_listed_values_not(self):
        keys = [u["key"] for u in self.a["uns"]]
        self.assertIn("params/method", keys)
        self.assertIn("Braak_colors", keys)
        self.assertNotIn("#111111", self.md + json.dumps(self.r))

    def test_markdown_has_the_sections(self):
        for heading in ("## AnnData:", "### Matrices", "### obs:", "### var:",
                        "### Embeddings and graphs", "### uns:"):
            self.assertIn(heading, self.md)


@unittest.skipUnless(HAVE, "anndata / h5py / scipy not installed")
class TestOtherLayouts(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = pathlib.Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def run_on(self, name):
        return ins.inspect(self.base / name, data_dir=self.base, out_dir=self.base / "out")

    def test_dense_scaled_matrix_is_called_scaled(self):
        x = np.random.default_rng(1).normal(size=(30, 10)).astype(np.float32)
        write(ad.AnnData(X=x), self.base / "scaled.h5ad")
        r, _ = self.run_on("scaled.h5ad")
        m = r["anndata"]["matrices"][0]
        self.assertEqual(m["storage"], "dense")
        self.assertIn("scaled or centred", m["values"])

    def test_integer_matrix_is_called_counts(self):
        x = sp.csr_matrix(np.random.default_rng(2).poisson(3, size=(30, 10)).astype(np.int32))
        write(ad.AnnData(X=x), self.base / "ints.h5ad")
        r, _ = self.run_on("ints.h5ad")
        self.assertIn("integer type", r["anndata"]["matrices"][0]["values"])

    def test_atac_peaks_and_binarised_matrix(self):
        x = sp.csr_matrix((np.random.default_rng(3).random((30, 12)) > 0.7).astype(np.float32))
        var = pd.DataFrame(index=[f"chr{1 + i % 3}:{1000 * i}-{1000 * i + 500}"
                                  for i in range(12)])
        write(ad.AnnData(X=x, var=var), self.base / "atac.h5ad")
        r, _ = self.run_on("atac.h5ad")
        a = r["anndata"]
        self.assertEqual(a["var"]["index"]["looks_like"], "genomic intervals (peaks), not genes")
        self.assertIn("binarised", a["matrices"][0]["values"])

    def test_plain_hdf5_gets_an_object_listing(self):
        with h5py.File(self.base / "tenx.h5", "w") as f:
            f.create_dataset("matrix/data", data=np.arange(5, dtype=np.int32))
            f.create_dataset("matrix/barcodes", data=np.array([b"SECRETCELL"] * 5))
        r, md = self.run_on("tenx.h5")
        self.assertIn("not AnnData", r["status"])
        paths = {o["path"] for o in r["hdf5_tree"]}
        self.assertIn("matrix/barcodes", paths)
        self.assertNotIn("SECRETCELL", md + json.dumps(r))


if __name__ == "__main__":
    unittest.main(verbosity=2)
