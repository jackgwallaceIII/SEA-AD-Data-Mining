"""
test_fetch.py — unit tests for scripts/step0b_fetch.py (browse and fetch).

Tier 1 (no network) for the selection logic, plus a small Tier 2 block
that checks the presets against the real pinned manifest and is skipped
when the cache is absent.

    python -m unittest discover -s tests -v

Why these tests: this script decides what gets downloaded. Its two real
hazards are fetching 157 GB by accident and silently fetching the wrong
thing, so the guards and the glob matching are what is covered here.
"""

import json
import pathlib
import sys
import tempfile
import types
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import step0b_fetch as fetch  # noqa: E402
from test_seaad import fake_manifest  # noqa: E402


def args(**kw):
    """An argparse-like namespace with the fetch defaults filled in."""
    base = dict(command="get", directory="SEA-AD-Multiregion-10X", files=[],
                kind=None, name=None, manifest=fetch.MANIFEST,
                dry_run=False, yes=True, force=False,
                max_gb=50.0, confirm_gb=1.0)
    base.update(kw)
    return types.SimpleNamespace(**base)


class TestWalk(unittest.TestCase):
    """Nesting depth differs by kind; the walker must handle both."""

    def test_expression_matrices_join_two_levels(self):
        m = fake_manifest()
        tree = m["file_listing"]["SEA-AD-Multiregion-10X"]["expression_matrices"]
        names = sorted(n for n, _ in fetch.walk(tree, []))
        self.assertIn("MTG-10X/raw", names)
        self.assertIn("MTG-10X/log2", names)

    def test_metadata_is_a_flat_name(self):
        m = fake_manifest()
        tree = m["file_listing"]["SEA-AD-Multiregion-10X"]["metadata"]
        names = sorted(n for n, _ in fetch.walk(tree, []))
        self.assertEqual(names, ["cell_metadata", "disease", "donor"])

    def test_every_leaf_is_reached(self):
        m = fake_manifest(regions=("A-10X", "B-10X"))
        tree = m["file_listing"]["SEA-AD-Multiregion-10X"]["expression_matrices"]
        self.assertEqual(len(list(fetch.walk(tree, []))), 4)  # 2 regions x 2 norms


class TestCatalog(unittest.TestCase):

    def setUp(self):
        self.m = fake_manifest()

    def test_flattens_every_kind(self):
        recs = fetch.catalog(self.m)
        kinds = {r["kind"] for r in recs}
        self.assertEqual(kinds, {"expression_matrices", "metadata"})

    def test_kind_filter(self):
        recs = fetch.catalog(self.m, kinds=("metadata",))
        self.assertTrue(all(r["kind"] == "metadata" for r in recs))
        self.assertEqual(len(recs), 3)

    def test_unknown_directory_yields_nothing(self):
        self.assertEqual(fetch.catalog(self.m, directory="NOPE"), [])

    def test_records_carry_size_and_hash(self):
        rec = next(r for r in fetch.catalog(self.m)
                   if r["file_name"] == "MTG-10X/raw")
        self.assertEqual(rec["size"], 1000)
        self.assertEqual(rec["md5"], "hash-MTG-10X-raw")


class TestLocalState(unittest.TestCase):
    """missing / partial / cached drives both the plan and the skip logic."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self._saved = fetch.CACHE_DIR
        fetch.CACHE_DIR = pathlib.Path(self._tmp.name)
        self.rec = {"relative_path": "a/b.h5ad", "size": 100}
        self.path = fetch.CACHE_DIR / "a/b.h5ad"
        self.path.parent.mkdir(parents=True)

    def tearDown(self):
        fetch.CACHE_DIR = self._saved
        self._tmp.cleanup()

    def test_missing(self):
        self.assertEqual(fetch.local_state(self.rec)[0], "missing")

    def test_partial(self):
        self.path.write_bytes(b"x" * 50)
        self.assertEqual(fetch.local_state(self.rec)[0], "partial")

    def test_cached(self):
        self.path.write_bytes(b"x" * 100)
        self.assertEqual(fetch.local_state(self.rec), ("cached", 100))


class TestSelectionGuards(unittest.TestCase):
    """The expensive mistakes: too much, or the wrong thing."""

    def setUp(self):
        self.m = fake_manifest()

    def test_bare_directory_without_kind_is_refused(self):
        # Against the real manifest this would be 157 GB.
        with self.assertRaises(SystemExit) as cm:
            fetch.resolve_selection(self.m, args(files=[]))
        self.assertIn("Refusing", str(cm.exception))

    def test_kind_alone_is_an_acceptable_narrowing(self):
        recs, _ = fetch.resolve_selection(self.m, args(files=[], kind="metadata"))
        self.assertEqual(len(recs), 3)

    def test_named_file(self):
        recs, _ = fetch.resolve_selection(self.m, args(files=["MTG-10X/raw"]))
        self.assertEqual([r["file_name"] for r in recs], ["MTG-10X/raw"])

    def test_glob(self):
        recs, _ = fetch.resolve_selection(self.m, args(files=["MTG-10X/*"]))
        self.assertEqual(sorted(r["file_name"] for r in recs),
                         ["MTG-10X/log2", "MTG-10X/raw"])

    def test_overlapping_globs_do_not_duplicate(self):
        recs, _ = fetch.resolve_selection(
            self.m, args(files=["MTG-10X/*", "MTG-10X/raw", "*/raw"]))
        names = [r["file_name"] for r in recs]
        self.assertEqual(len(names), len(set(names)))

    def test_unmatched_glob_is_an_error_not_a_silent_noop(self):
        # Fetching nothing because of a typo must not look like success.
        with self.assertRaises(SystemExit) as cm:
            fetch.resolve_selection(self.m, args(files=["ZZZ*"]))
        self.assertIn("ZZZ*", str(cm.exception))

    def test_unknown_directory_lists_the_known_ones(self):
        with self.assertRaises(SystemExit) as cm:
            fetch.resolve_selection(self.m, args(directory="NOPE", files=["x"]))
        self.assertIn("SEA-AD-Multiregion-10X", str(cm.exception))


class TestPresetResolution(unittest.TestCase):

    def test_preset_missing_files_fails_loudly(self):
        # A release that reorganized must not silently fetch a subset.
        thin = fake_manifest(regions=("MTG-10X",))
        del thin["file_listing"]["SEA-AD-Multiregion-10X"]["metadata"]["donor"]
        with self.assertRaises(SystemExit) as cm:
            fetch.resolve_selection(thin, args(command="preset", name="mtg"))
        self.assertIn("donor", str(cm.exception))

    def test_presets_are_explicit_file_lists(self):
        # Not a directory sweep: every entry names one file.
        for name, p in fetch.PRESETS.items():
            for item in p["files"]:
                self.assertEqual(len(item), 2, f"{name} entry is not (dir, file)")


class TestHuman(unittest.TestCase):

    def test_units(self):
        self.assertEqual(fetch.human(512), "512 B")
        self.assertEqual(fetch.human(2_000_000), "2.00 MB")
        self.assertEqual(fetch.human(15_190_788_656), "15.19 GB")


# ---- Tier 2: needs the real pinned manifest on disk --------------------

REAL = fetch.CACHE_DIR / fetch.MANIFEST


@unittest.skipUnless(REAL.exists(), f"pinned manifest not cached ({REAL})")
class TestAgainstRealManifest(unittest.TestCase):
    """Checks the presets actually resolve in the pinned release."""

    @classmethod
    def setUpClass(cls):
        cls.m = json.loads(REAL.read_text())

    def test_every_preset_resolves(self):
        for name in fetch.PRESETS:
            recs, _ = fetch.resolve_selection(
                self.m, args(command="preset", name=name))
            self.assertEqual(len(recs), len(fetch.PRESETS[name]["files"]))

    def test_mtg_preset_matches_the_recorded_matrix(self):
        recs, _ = fetch.resolve_selection(
            self.m, args(command="preset", name="mtg"))
        matrix = next(r for r in recs if r["file_name"] == "MTG-10X/raw")
        self.assertEqual(matrix["md5"], "affce583f5b2db357022d289a9e9892a")
        self.assertEqual(matrix["size"], 15_190_788_656)

    def test_metadata_preset_excludes_the_matrix(self):
        recs, _ = fetch.resolve_selection(
            self.m, args(command="preset", name="mtg-metadata"))
        self.assertTrue(all(r["kind"] == "metadata" for r in recs))

    def test_all_ten_regions_present(self):
        em = self.m["file_listing"]["SEA-AD-Multiregion-10X"]["expression_matrices"]
        self.assertEqual(len(em), 10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
