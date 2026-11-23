"""
test_stores.py — unit tests for scripts/descriptors.py and scripts/stores.py.

No network, no real data. The Zarr fixture is hand-written consolidated
metadata shaped like the real SEA-AD SpatialData store (same groups and
element names as a listing of caudate_nucleus/.../combined_spatialdata_zarr/
shows), so every number in the expected report is known. Run identifiers
and gene names are planted and must not reach a report.

    python -m unittest discover -s tests -t tests -v

Why these tests: `get --descriptors` must fetch only descriptor files and
nothing else from a multi-gigabyte folder, through the normal verified
downloader; the store report must get element shapes, coordinate systems
and pixel scales right and must never contain a category level or a run
identifier.
"""

import json
import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import descriptors  # noqa: E402
import index  # noqa: E402
import inspect_data as ins  # noqa: E402
import stores  # noqa: E402
from test_images import FakeObjectStore  # noqa: E402

ST = "sea-ad-spatial-transcriptomics"
STORE = "caudate_nucleus/H00.00.000/1/combined_spatialdata_zarr/"
BUNDLE = "caudate_nucleus/H00.00.000/1/xenium_output_and_alignment/output-XETG00000__1/"


def arr(shape, dtype="<f4", chunks=None, filters=None):
    return {"shape": shape, "dtype": dtype, "chunks": chunks or shape, "zarr_format": 2,
            "compressor": {"id": "blosc", "cname": "zstd"}, "filters": filters}


GLOBAL = {"name": "global", "axes": [{"name": "x", "type": "space", "unit": "unit"}]}
META = {
    ".zgroup": {"zarr_format": 2},
    ".zattrs": {"spatialdata_attrs": {"version": "0.1"}},
    "images/.zgroup": {"zarr_format": 2},
    "images/morphology_focus/.zattrs": {
        "multiscales": [{"version": "0.4", "axes": [{"name": "c", "type": "channel"},
                                                    {"name": "y", "type": "space"},
                                                    {"name": "x", "type": "space"}],
                         "datasets": [{"path": "0", "coordinateTransformations": [
                                          {"type": "scale", "scale": [1.0, 1.0, 1.0]}]},
                                      {"path": "1", "coordinateTransformations": [
                                          {"type": "scale", "scale": [1.0, 2.0, 2.0]}]}],
                         "coordinateTransformations": [
                             {"type": "identity", "input": {"name": "cyx"}, "output": GLOBAL}]}],
        "omero": {"channels": [{"label": "DAPI"}, {"label": "ATP1A1/CD45/E-Cadherin"}]}},
    "images/morphology_focus/0/.zarray": arr([2, 40000, 30000], "<u2", [1, 4096, 4096]),
    "images/morphology_focus/1/.zarray": arr([2, 20000, 15000], "<u2", [1, 4096, 4096]),
    "labels/.zgroup": {"zarr_format": 2},
    "labels/cell_labels/.zattrs": {
        "multiscales": [{"version": "0.4", "axes": [{"name": "y"}, {"name": "x"}],
                         "datasets": [{"path": "0"}],
                         "coordinateTransformations": [
                             {"type": "affine", "affine": [[0.2, 0, 1], [0, 0.2, 2], [0, 0, 1]],
                              "input": {"name": "yx"}, "output": GLOBAL}]}]},
    "labels/cell_labels/0/.zarray": arr([40000, 30000], "<u4"),
    "points/.zgroup": {"zarr_format": 2},
    "points/transcripts/.zattrs": {
        "coordinateTransformations": [{"type": "scale", "scale": [1 / 0.2125, 1 / 0.2125, 1.0],
                                       "input": {"name": "xyz"}, "output": GLOBAL}],
        "spatialdata_attrs": {"feature_key": "feature_name", "instance_key": "cell_id"}},
    "shapes/.zgroup": {"zarr_format": 2},
    "shapes/cell_boundaries/.zattrs": {
        "coordinateTransformations": [{"type": "identity", "output": GLOBAL}],
        "encoding-type": "ngff:shapes"},
    "tables/.zgroup": {"zarr_format": 2},
    "tables/mapped_table/.zattrs": {"encoding-type": "anndata"},
    "tables/mapped_table/X/.zattrs": {"encoding-type": "csr_matrix", "shape": [5000, 300]},
    "tables/mapped_table/X/data/.zarray": arr([123456], "<f4"),
    "tables/mapped_table/layers/counts/.zattrs": {"encoding-type": "csr_matrix",
                                                  "shape": [5000, 300]},
    "tables/mapped_table/layers/counts/data/.zarray": arr([123456], "<u4"),
    "tables/mapped_table/obs/.zattrs": {"_index": "_index", "encoding-type": "dataframe",
                                        "column-order": ["CDM_Subclass_name", "cell_area",
                                                         "cell_id"]},
    "tables/mapped_table/obs/_index/.zarray": arr([5000], "|O", filters=[{"id": "vlen-utf8"}]),
    "tables/mapped_table/obs/CDM_Subclass_name/.zattrs": {"encoding-type": "categorical"},
    "tables/mapped_table/obs/CDM_Subclass_name/categories/.zarray":
        arr([24], "|O", filters=[{"id": "vlen-utf8"}]),
    "tables/mapped_table/obs/CDM_Subclass_name/codes/.zarray": arr([5000], "|i1"),
    "tables/mapped_table/obs/cell_area/.zarray": arr([5000], "<f8"),
    "tables/mapped_table/obs/cell_id/.zarray": arr([5000], "|O", filters=[{"id": "vlen-utf8"}]),
    "tables/mapped_table/var/.zattrs": {"_index": "_index", "column-order": ["gene_ids"]},
    "tables/mapped_table/var/_index/.zarray": arr([300], "|O", filters=[{"id": "vlen-utf8"}]),
    "tables/mapped_table/var/gene_ids/.zarray": arr([300], "|O", filters=[{"id": "vlen-utf8"}]),
    "tables/mapped_table/obsm/spatial/.zarray": arr([5000, 2], "<f8"),
    "tables/mapped_table/uns/mapping_config/.zgroup": {"zarr_format": 2},
    "tables/mapped_table/uns/taxonomy_ROI/.zarray": arr([], "|O"),
}

EXPERIMENT = {"major_version": 5, "minor_version": 0, "run_name": "SECRETRUN",
              "region_name": "SECRETREGION", "instrument_sn": "SECRETSN",
              "slide_id": "0036828", "run_start_time": "2025-03-05T20:39:43Z",
              "panel_name": "Xenium Prime 5K Human Pan Tissue", "panel_organism": "Human",
              "pixel_size": 0.2125, "num_cells": 81234, "transcripts_per_cell": 456,
              "analysis_sw_version": "xenium-3.1", "experiment_uuid": "SECRETUUID",
              "xenium_ranger": {"version": "xenium-ranger-3.1", "command_line": "SECRETCMD",
                                "run_id": "SECRETRUNID"},
              "images": {"morphology_filepath": "morphology.ome.tif",
                         "morphology_focus_filepath": "morphology_focus/"}}
PANEL = {"metadata": {"sha256_checksum": "x", "ui_props": {"visible": True}},
         "payload": {
             "chemistry": {"decoder": "d", "version": "v1"},
             "customer": {"email": "SECRET.customer@example.org", "name": "SECRET Customer"},
             "designer": {"date_created": "2024-01-01", "version": "2",
                          "designed_by": {"email": "SECRET.designer@example.org",
                                          "name": "SECRET Designer"}},
             "panel": {"description": "d", "identity": {"design_id": "SECRETDESIGN",
                                                        "name": "hBrain_464g", "version": "1"},
                       "num_gene_targets": 5, "species": "Human", "tissue": "brain",
                       "type": {"descriptor": "custom"}},
             "spec_version": "3.0",
             "targets": [{"type": {"descriptor": "gene", "data": {"name": f"SECRETGENE{i}"}}}
                         for i in range(5)]
             + [{"type": {"descriptor": "negative_control"}}] * 2}}


class TestPick(unittest.TestCase):

    def test_consolidated_file_alone(self):
        objs = [("zmetadata", 92402, "e"), (".zgroup", 24, "e"), ("images/0/0.0.0", 9, "e"),
                ("points/transcripts/points.parquet/part.0.parquet", 10 ** 8, "e")]
        self.assertEqual([o[0] for o in descriptors.pick("zarr", objs)], ["zmetadata"])

    def test_metadata_files_when_not_consolidated(self):
        objs = [(".zgroup", 24, "e"), ("images/a/.zattrs", 50, "e"), ("images/a/0/.zarray", 30, "e"),
                ("images/a/0/0.0.0", 10 ** 6, "e"), ("shapes/s/shapes.parquet", 10 ** 6, "e")]
        self.assertEqual([o[0] for o in descriptors.pick("zarr", objs)],
                         [".zgroup", "images/a/.zattrs", "images/a/0/.zarray"])

    def test_xenium_two_files_never_data(self):
        objs = [("experiment.xenium", 2692, "e"), ("gene_panel.json", 233598, "e"),
                ("metrics_summary.csv", 2053, "e"), ("transcripts.parquet", 10 ** 8, "e"),
                ("morphology.ome.tif", 10 ** 9, "e")]
        self.assertEqual(sorted(o[0] for o in descriptors.pick("xenium bundle", objs)),
                         ["experiment.xenium", "gene_panel.json"])


class TestZarrReport(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.z = stores.describe_zarr(META, "test")
        cls.md = "\n".join(stores.zarr_markdown(cls.z))

    def test_elements_found(self):
        self.assertEqual({g: sorted(e) for g, e in self.z["elements"].items()},
                         {"images": ["morphology_focus"], "labels": ["cell_labels"],
                          "points": ["transcripts"], "shapes": ["cell_boundaries"],
                          "tables": ["mapped_table"]})
        self.assertEqual(self.z["coordinate_systems"], ["global"])
        self.assertEqual(self.z["spatialdata_version"], "0.1")

    def test_image_levels_channels_and_pixel_scale(self):
        im = self.z["elements"]["images"]["morphology_focus"]
        self.assertEqual([lv["shape"] for lv in im["levels"]],
                         [[2, 40000, 30000], [2, 20000, 15000]])
        self.assertEqual(im["levels"][0]["dtype"], "uint16")
        self.assertEqual(im["channels"], ["DAPI", "ATP1A1/CD45/E-Cadherin"])
        self.assertEqual(im["transforms"], [("global", "identity", None)])

    def test_units_explained_from_identity_image_and_scaled_points(self):
        self.assertIn("`global` is measured in image pixels", self.md)
        self.assertIn("1 ÷ 4.706 = 0.2125", self.md)
        self.assertIn("micrometres", self.md)

    def test_scales_rounded(self):
        z = stores.describe_zarr({**META, "images/morphology_focus/.zattrs": {
            **META["images/morphology_focus/.zattrs"],
            "multiscales": [{**META["images/morphology_focus/.zattrs"]["multiscales"][0],
                             "datasets": [{"path": "0", "coordinateTransformations": [
                                 {"type": "scale", "scale": [1.0, 2.000083843380565, 2.0]}]}]}]}},
            "t")
        md = "\n".join(stores.zarr_markdown(z))
        self.assertIn("[1, 2, 2]", md)
        self.assertNotIn("2.000083", md)

    def test_affine_named_not_printed(self):
        lab = self.z["elements"]["labels"]["cell_labels"]
        self.assertEqual(lab["transforms"], [("global", "affine", None)])
        self.assertNotIn("[[0.2", self.md)

    def test_table_shape_types_and_category_counts(self):
        t = self.z["elements"]["tables"]["mapped_table"]
        self.assertEqual((t["n_obs"], t["n_vars"]), (5000, 300))
        self.assertEqual(t["X"]["storage"], "csr sparse")
        self.assertEqual(t["layers"]["counts"]["dtype"], "uint32")
        types = {c["name"]: c["type"] for c in t["obs"]["columns"]}
        self.assertEqual(types, {"CDM_Subclass_name": "category (24 levels)",
                                 "cell_area": "float64", "cell_id": "text"})
        self.assertEqual(t["obsm"], {"spatial": [5000, 2]})
        self.assertIn("5,000 cells × 300 genes", self.md)

    def test_folder_of_metadata_files_gives_the_same_answer(self):
        with tempfile.TemporaryDirectory() as tmp:
            for rel, content in META.items():
                p = pathlib.Path(tmp) / rel
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(json.dumps(content))
            again = stores.describe_zarr(stores.metadata_from_folder(tmp), "test")
        self.assertEqual(again["elements"], self.z["elements"])


class TestXenium(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.data = pathlib.Path(self.tmp.name) / "data_cache"
        base = self.data / ST / BUNDLE
        base.mkdir(parents=True)
        (base / "experiment.xenium").write_text(json.dumps(EXPERIMENT))
        (base / "gene_panel.json").write_text(json.dumps(PANEL))
        self.base = base

    def tearDown(self):
        self.tmp.cleanup()

    def inspect(self, name):
        return ins.inspect(self.base / name, data_dir=self.data,
                           out_dir=pathlib.Path(self.tmp.name) / "out")

    def test_whitelisted_fields_only(self):
        r, md = self.inspect("experiment.xenium")
        x = r["xenium"]
        self.assertEqual(x["fields"]["pixel_size"], 0.2125)
        self.assertEqual(x["fields"]["num_cells"], 81234)
        self.assertIn("instrument_sn", x["withheld_fields"])
        self.assertEqual(x["file_groups"]["images"],
                         ["morphology_filepath", "morphology_focus_filepath"])
        blob = md + json.dumps(r)
        for secret in ("SECRETRUN", "SECRETREGION", "SECRETSN", "SECRETUUID", "0036828",
                       "2025-03-05"):
            self.assertNotIn(secret, blob, secret)

    def test_panel_counts_targets_not_names(self):
        r, md = self.inspect("gene_panel.json")
        p = r["panel"]
        self.assertEqual((p["n_targets"], p["targets_by_type"]),
                         (7, {"gene": 5, "negative_control": 2}))
        self.assertEqual((p["panel_name"], p["organism"], p["tissue"]),
                         ("hBrain_464g", "Human", "brain"))
        self.assertIn("matches the declared gene count", md)
        blob = md + json.dumps(r)
        for secret in ("SECRETGENE", "SECRET.customer", "SECRET.designer", "SECRET Customer",
                       "SECRET Designer", "SECRETDESIGN"):
            self.assertNotIn(secret, blob, secret)

    def test_xenium_ranger_version_only(self):
        r, md = self.inspect("experiment.xenium")
        self.assertEqual(r["xenium"]["fields"]["xenium_ranger_version"], "xenium-ranger-3.1")
        self.assertNotIn("SECRETCMD", md + json.dumps(r))
        self.assertNotIn("SECRETRUNID", md + json.dumps(r))


class TestGetDescriptors(unittest.TestCase):
    """End to end: catalog -> folder listing -> only descriptors fetched, verified."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self.tmp.name)
        self.cat, self.data = base / "catalog", base / "data_cache"
        zmeta = json.dumps({"metadata": META, "zarr_consolidated_format": 1}).encode()
        self.objects = {
            STORE + "zmetadata": zmeta, STORE + ".zgroup": b'{"zarr_format": 2}',
            STORE + "images/morphology_focus/0/0.0.0": b"SECRETPIXELS" * 100,
            STORE + "points/transcripts/points.parquet/part.0.parquet": b"PAR1" * 100,
            BUNDLE + "experiment.xenium": json.dumps(EXPERIMENT).encode(),
            BUNDLE + "gene_panel.json": json.dumps(PANEL).encode(),
            BUNDLE + "transcripts.parquet": b"PAR1" * 1000,
            BUNDLE + "metrics_summary.csv": b"num_cells\n1\n",
            "caudate_nucleus/H00.00.000/1/IF_images/plain.tif": b"II*\x00" + b"\x00" * 60,
        }
        self.store = FakeObjectStore(self.objects)
        self.store.etag = lambda key: __import__("hashlib").md5(self.objects[key]).hexdigest()
        self.store.list_objects_v2 = self._list(self.store.list_objects_v2)
        index.refresh(ST, client=self.store, catalog_dir=self.cat, out=lambda *a: None)
        self.lines = []

    def _list(self, inner):
        """FakeS3 lists 'etag-<n>'; give the real MD5 so verification can pass."""
        import hashlib

        def wrapped(**kw):
            out = inner(**kw)
            for o in out["Contents"]:
                o["ETag"] = f'"{hashlib.md5(self.objects[o["Key"]]).hexdigest()}"'
            return out
        return wrapped

    def tearDown(self):
        self.tmp.cleanup()

    def get_object(self, Bucket, Key, Range=None, IfMatch=None):
        import io
        return {"Body": io.BytesIO(self.objects[Key])}

    def run_get(self, match, **kw):
        self.store.get_object = self.get_object
        return descriptors.run(match, client=self.store, catalog_dir=self.cat,
                               data_dir=self.data, out=self.lines.append, **kw)

    def test_only_descriptors_land_on_disk(self):
        code = self.run_get("*/H00.00.000/1/*/", yes=True)
        self.assertEqual(code, 0, "\n".join(self.lines))
        got = sorted(p.relative_to(self.data / ST).as_posix()
                     for p in (self.data / ST).rglob("*") if p.is_file())
        self.assertEqual(got, sorted([BUNDLE + "experiment.xenium", BUNDLE + "gene_panel.json",
                                      STORE + "zmetadata"]))
        ledger = (self.data / "_ledger.csv").read_text()
        self.assertEqual(ledger.count("md5 match"), 3)

    def test_dry_run_and_preview(self):
        self.assertEqual(self.run_get("combined_spatialdata_zarr", dry_run=True), 0)
        text = "\n".join(self.lines)
        self.assertIn("only these 1 file(s)", text)
        self.assertFalse(self.data.exists())

    def test_inspect_the_fetched_zmetadata(self):
        self.run_get("combined_spatialdata_zarr", yes=True)
        r, md = ins.inspect(self.data / ST / STORE / "zmetadata", data_dir=self.data,
                            out_dir=pathlib.Path(self.tmp.name) / "out")
        self.assertEqual(r["zarr"]["coordinate_systems"], ["global"])
        self.assertIn("### tables/mapped_table (AnnData)", md)
        self.assertIn("(no extension)", md)

    def test_plain_files_are_refused(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_get("plain.tif")
        self.assertIn("Zarr stores and Xenium bundles only", str(cm.exception))


if __name__ == "__main__":
    unittest.main(verbosity=2)
