"""
test_sources.py — unit tests for scripts/sources.py (read-only S3 listing).

No network. Two kinds of fake S3 are used:

- botocore's Stubber, wrapped around a real anonymous client. It checks
  the exact API calls and parameters, and FAILS on any call it was not
  told to expect — which is how "this module only lists" is enforced.
- FakeS3, a few lines that answer list_objects_v2 from an in-memory key
  list with real Prefix/Delimiter/pagination behaviour (page size 3, so
  every test crosses page boundaries). Used for whole-walk behaviour.

    python -m unittest discover -s tests -t tests -v

Why these tests: the walker decides what the catalog sees. The ways it
can be wrong are (a) dropping keys at a page boundary, (b) counting a key
twice, (c) listing inside a folder that should be skipped, (d) reporting
a capped count as exact, and (e) sending credentials.
"""

import datetime
import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import pathrules  # noqa: E402
import sources  # noqa: E402

SC = "sea-ad-single-cell-profiling"
ST = "sea-ad-spatial-transcriptomics"
NP = "sea-ad-quantitative-neuropathology"
T0 = datetime.datetime(2026, 6, 22, tzinfo=datetime.timezone.utc)


class FakeS3:
    """Answers list_objects_v2 like S3 does, from {key: size}."""

    PAGE = 3

    def __init__(self, keys):
        self.keys = dict(keys)
        self.calls = []

    def list_objects_v2(self, Bucket, Prefix="", Delimiter=None,
                        ContinuationToken=None):
        self.calls.append({"Prefix": Prefix, "Delimiter": Delimiter})
        entries = []                       # ('P', prefix) or ('K', key)
        seen = set()
        for k in sorted(self.keys):
            if not k.startswith(Prefix):
                continue
            rest = k[len(Prefix):]
            if Delimiter and Delimiter in rest:
                p = Prefix + rest.split(Delimiter, 1)[0] + Delimiter
                if p not in seen:
                    seen.add(p)
                    entries.append(("P", p))
            else:
                entries.append(("K", k))
        start = int(ContinuationToken or 0)
        page = entries[start:start + self.PAGE]
        out = {"IsTruncated": start + self.PAGE < len(entries),
               "CommonPrefixes": [{"Prefix": v} for t, v in page if t == "P"],
               "Contents": [{"Key": v, "Size": self.keys[v],
                             "ETag": f'"etag-{len(v)}"', "LastModified": T0}
                            for t, v in page if t == "K"]}
        if out["IsTruncated"]:
            out["NextContinuationToken"] = str(start + self.PAGE)
        return out


def walk(bucket, keys, **kw):
    fake = FakeS3(keys)
    src = sources.S3Source(bucket, client=fake)
    return list(src.walk(**kw)), fake, src


def accounted(records):
    """Total size covered by the records (files + complete units)."""
    return sum(r.size if isinstance(r, sources.FileObject) else r.total_size
               for r in records if not isinstance(r, sources.SkippedFolder))


# ---- real keys from the 2026-11-21 listings ------------------------------

PLOTS = {
    "MTG/RNAseq/pseudoprogression-plots/A1BG-AS1/overview.jpg": 10,
    "MTG/RNAseq/pseudoprogression-plots/A1BG-AS1/overview.svg": 20,
    "MTG/RNAseq/pseudoprogression-plots/A1BG-AS1/overview_subclass.jpg": 30,
    "MTG/RNAseq/pseudoprogression-plots/AC008840.1/supertype-Sst.jpg": 40,
    "MTG/RNAseq/pseudoprogression-plots/AC008840.1/supertype-VLMC.jpg": 50,
}
MTG_RNA = {
    "MTG/RNAseq/": 0,                                                # marker
    "MTG/RNAseq/Changelog.html": 1,
    "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad": 1000,
    "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei_metadata.2026-06-22.csv": 100,
    "MTG/RNAseq/previous_objects/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad": 900,
    **PLOTS,
}
ZARR = "caudate_nucleus/H19.33.004/1417642334/combined_spatialdata_zarr/"
XENIUM = ("caudate_nucleus/H19.33.004/1417642334/xenium_output_and_alignment/"
          "output-XETG00044__0031423__1417642334__20250319__192442/")
SPATIAL = {
    ZARR + ".zgroup": 1,
    ZARR + "images/morphology_focus/0/0/0/0": 500,
    ZARR + "points/transcripts/points.parquet/part.0.parquet": 300,
    ZARR + "tables/mapped_table/X/data/0": 200,
    XENIUM + "cell_feature_matrix.h5": 70,
    XENIUM + "cell_feature_matrix.zarr.zip": 80,
    "caudate_nucleus/H19.33.004/1417642334/IF_images/DAPI.tif": 5000,
}
HIVE = "Multiregion_2026/previous_objects/parquet/Multiregion_2026/counts_by_gene/"
PARQUET = {
    HIVE + "gene=A1BG-AS1/data.parquet": 7,
    HIVE + "gene=A1BG/data.parquet": 8,
    HIVE + "gene=A1CF/data.parquet": 9,
    HIVE + "gene=A2M/data.parquet": 10,
    # SYNTHETIC: a file directly inside the partitioned folder. Not seen in
    # the real listing, but it is the case that would be double-counted
    # if the walker yielded this level's files before collapsing it.
    HIVE + "_common_metadata": 5,
    "Multiregion_2026/previous_objects/parquet/Multiregion_2026/manifest.json": 1,
}
NEURO = {
    "simple_viewer/dorsolateral-prefrontal-cortex/H19.33.004.html": 1,
    "simple_viewer/openseadragon/images/button_grouphover.png": 1,
    "temp/images/ABAT.svg": 1,
    "full-tabular-data/CaH_all-values_aibs-ids_20250922.csv": 40,
    "middle-temporal-gyrus/H19.33.004/H19.33.004-A06-NEUN/H19.33.004-A06-NEUN.svs": 800,
}


class TestAnonymousClient(unittest.TestCase):

    def test_requests_are_unsigned(self):
        from botocore import UNSIGNED
        client = sources.anonymous_client()
        self.assertIs(client.meta.config.signature_version, UNSIGNED)
        self.assertEqual(client.meta.region_name, "us-west-2")


class TestStubbedApiCalls(unittest.TestCase):
    """Exact calls against a real (stubbed) boto3 client."""

    def setUp(self):
        from botocore.stub import Stubber
        self.client = sources.anonymous_client()
        self.stub = Stubber(self.client)
        self.src = sources.S3Source(SC, client=self.client)

    def obj(self, key, size):
        return {"Key": key, "Size": size, "ETag": '"abc123"',
                "LastModified": T0}

    def test_list_level_follows_continuation_tokens(self):
        base = {"Bucket": SC, "Prefix": "MTG/", "Delimiter": "/"}
        self.stub.add_response("list_objects_v2", {
            "IsTruncated": True, "NextContinuationToken": "tok1",
            "CommonPrefixes": [{"Prefix": "MTG/ATACseq/"}],
            "Contents": [self.obj("MTG/", 0)],
        }, base)
        self.stub.add_response("list_objects_v2", {
            "IsTruncated": False,
            "CommonPrefixes": [{"Prefix": "MTG/RNAseq/"}],
            "Contents": [self.obj("MTG/README.txt", 5)],
        }, {**base, "ContinuationToken": "tok1"})
        with self.stub:
            folders, files = self.src.list_level("MTG/")
        self.stub.assert_no_pending_responses()
        self.assertEqual(folders, ["MTG/ATACseq/", "MTG/RNAseq/"])
        self.assertEqual([f.key for f in files], ["MTG/README.txt"])
        self.assertEqual(self.src.folder_markers, 1)

    def test_etag_quotes_are_stripped_and_time_is_iso(self):
        self.stub.add_response("list_objects_v2", {
            "IsTruncated": False, "Contents": [self.obj("a.csv", 3)]},
            {"Bucket": SC, "Prefix": "", "Delimiter": "/"})
        with self.stub:
            _, files = self.src.list_level("")
        self.assertEqual(files[0].etag, "abc123")
        self.assertEqual(files[0].last_modified, "2026-06-22T00:00:00+00:00")

    def test_a_skipped_prefix_makes_no_api_call_at_all(self):
        # No responses queued: any S3 call would raise inside the stubber.
        src = sources.S3Source(NP, client=self.client)
        with self.stub:
            recs = list(src.walk("simple_viewer/"))
        self.assertEqual(len(recs), 1)
        self.assertIsInstance(recs[0], sources.SkippedFolder)

    def test_module_exposes_no_download_or_read(self):
        public = {n for n in dir(sources.S3Source) if not n.startswith("_")}
        self.assertEqual(public, {"list_level", "summarise", "walk"})


class TestWalk(unittest.TestCase):

    def test_every_key_accounted_once_across_page_boundaries(self):
        recs, _, src = walk(SC, MTG_RNA, prefix="MTG/")
        real = {k: v for k, v in MTG_RNA.items() if not k.endswith("/")}
        self.assertEqual(accounted(recs), sum(real.values()))
        files = [r.key for r in recs if isinstance(r, sources.FileObject)]
        self.assertEqual(len(files), len(set(files)))
        self.assertEqual(src.folder_markers, 1)

    def test_deprecated_files_are_listed_not_skipped(self):
        recs, _, _ = walk(SC, MTG_RNA, prefix="MTG/")
        keys = {r.key for r in recs if isinstance(r, sources.FileObject)}
        self.assertIn("MTG/RNAseq/previous_objects/"
                      "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad", keys)

    def test_plot_folder_becomes_one_complete_unit(self):
        recs, _, _ = walk(SC, MTG_RNA, prefix="MTG/")
        units = [r for r in recs if isinstance(r, sources.CollapsedUnit)]
        self.assertEqual(len(units), 1)
        u = units[0]
        self.assertEqual(u.prefix, "MTG/RNAseq/pseudoprogression-plots/")
        self.assertEqual((u.n_objects, u.total_size, u.complete), (5, 150, True))
        self.assertEqual(u.product, "plots")

    def test_cap_makes_the_unit_a_lower_bound(self):
        recs, _, _ = walk(SC, PLOTS, prefix="MTG/", max_keys_per_unit=2)
        u = next(r for r in recs if isinstance(r, sources.CollapsedUnit))
        self.assertFalse(u.complete)
        self.assertEqual(u.n_objects, 2)
        self.assertLess(u.total_size, sum(PLOTS.values()))

    def test_zarr_and_xenium_collapse_and_ordinary_file_does_not(self):
        recs, _, _ = walk(ST, SPATIAL, prefix="caudate_nucleus/")
        units = {r.prefix: r for r in recs if isinstance(r, sources.CollapsedUnit)}
        self.assertEqual(set(units), {ZARR, XENIUM})
        self.assertEqual(units[ZARR].total_size, 1001)     # parquet inside counted once
        self.assertEqual(units[XENIUM].n_objects, 2)
        files = [r.key for r in recs if isinstance(r, sources.FileObject)]
        self.assertEqual(files, ["caudate_nucleus/H19.33.004/1417642334/IF_images/DAPI.tif"])
        self.assertEqual(accounted(recs), sum(SPATIAL.values()))

    def test_never_lists_inside_a_collapsed_folder_by_level(self):
        _, fake, _ = walk(ST, SPATIAL, prefix="caudate_nucleus/")
        level_calls = [c["Prefix"] for c in fake.calls if c["Delimiter"]]
        self.assertFalse(any(p.startswith(ZARR) for p in level_calls))

    def test_hive_partitions_collapse_without_double_counting(self):
        recs, _, _ = walk(SC, PARQUET, prefix="Multiregion_2026/")
        units = [r for r in recs if isinstance(r, sources.CollapsedUnit)]
        self.assertEqual([u.prefix for u in units], [HIVE])
        self.assertEqual(units[0].n_objects, 5)
        self.assertEqual(accounted(recs), sum(PARQUET.values()))

    def test_skipped_folders_are_reported_and_never_listed(self):
        recs, fake, _ = walk(NP, NEURO, prefix="")
        skipped = {r.prefix: r.reason for r in recs
                   if isinstance(r, sources.SkippedFolder)}
        self.assertEqual(skipped, {"simple_viewer/": "web viewer assets",
                                   "temp/": "scratch folder"})
        listed = [c["Prefix"] for c in fake.calls]
        self.assertFalse(any(p.startswith(("simple_viewer/", "temp/"))
                             for p in listed))
        self.assertEqual(accounted(recs), 40 + 800)

    def test_empty_prefix_yields_nothing_rather_than_failing(self):
        recs, _, _ = walk(SC, MTG_RNA, prefix="NOPE/")
        self.assertEqual(recs, [])


class TestFolderRuleAgreesWithClassify(unittest.TestCase):
    """The walker's folder decision and classify()'s unit must match."""

    def test_every_key_in_a_collapsed_unit_classifies_to_that_unit(self):
        for bucket, keys, prefix in ((SC, MTG_RNA, "MTG/"),
                                     (ST, SPATIAL, "caudate_nucleus/"),
                                     (SC, PARQUET, "Multiregion_2026/")):
            recs, _, _ = walk(bucket, keys, prefix=prefix)
            for u in (r for r in recs if isinstance(r, sources.CollapsedUnit)):
                inside = [k for k in keys if k.startswith(u.prefix)
                          and not self.sibling_only(u, k)]
                self.assertTrue(inside)
                for k in inside:
                    self.assertEqual(pathrules.classify(bucket, k).unit,
                                     u.prefix, k)

    @staticmethod
    def sibling_only(unit, key):
        """A file directly inside a Hive table folder.

        Its own path has no 'key=value' segment, so classify() — which
        sees one key at a time — cannot know it belongs to the table.
        Only the walker, which sees the sibling folders, can. The catalog
        therefore takes unit membership from the walker. See the next test.
        """
        return (unit.format == "partitioned parquet"
                and "/" not in key[len(unit.prefix):])

    def test_known_limit_classify_cannot_see_siblings(self):
        # Documented, not hidden: key-only classification leaves this file
        # as its own unit, while the walker puts it in the table.
        k = HIVE + "_common_metadata"
        self.assertEqual(pathrules.classify(SC, k).unit, k)
        recs, _, _ = walk(SC, PARQUET, prefix="Multiregion_2026/")
        files = {r.key for r in recs if isinstance(r, sources.FileObject)}
        self.assertNotIn(k, files)

    def test_files_yielded_by_the_walk_are_not_collapsed_by_classify(self):
        recs, _, _ = walk(ST, SPATIAL, prefix="caudate_nucleus/")
        for r in recs:
            if isinstance(r, sources.FileObject):
                self.assertFalse(pathrules.classify(ST, r.key).collapsed, r.key)


class TestCli(unittest.TestCase):

    def test_whole_bucket_walk_is_refused(self):
        with self.assertRaises(SystemExit) as cm:
            sources.main(["walk", SC])
        self.assertIn("Refusing", str(cm.exception))


if __name__ == "__main__":
    unittest.main(verbosity=2)
