"""
test_download.py — unit tests for scripts/download.py (verified download).

No network. FakeObjects serves real bytes with ETags computed the way S3
computes them (MD5 for single-part; MD5-of-part-MD5s + '-N' for multipart),
independently of download.ETagHasher, so the checksum code is checked
against the real formula rather than against itself. Every test writes
into its own temporary data_cache/.

    python -m unittest discover -s tests -t tests -v

Why these tests: this is the code that puts bytes on disk. The failures
that matter are a corrupted or truncated file left under its real name,
a file that changed on S3 being fetched as if it were the previewed one,
a download starting without confirmation, and a ledger that claims a
stronger check than was done.
"""

import csv
import hashlib
import io
import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import download  # noqa: E402

SC = "sea-ad-single-cell-profiling"


def s3_etag(data, part_size=None):
    """The ETag S3 would give these bytes."""
    if part_size is None:
        return hashlib.md5(data).hexdigest()
    parts = [data[i:i + part_size] for i in range(0, len(data), part_size)]
    joined = b"".join(hashlib.md5(p).digest() for p in parts)
    return f"{hashlib.md5(joined).hexdigest()}-{len(parts)}"


class Body:
    def __init__(self, data, truncate_at=None):
        self.buf = io.BytesIO(data[:truncate_at] if truncate_at else data)

    def read(self, n=-1):
        return self.buf.read(n)


class FakeObjects:
    """head_object / get_object over in-memory bytes. Records every call."""

    def __init__(self):
        self.objects = {}           # key -> (bytes, etag, part_size)
        self.calls = []
        self.corrupt = set()
        self.truncate = {}
        self.no_part_head = False

    def put(self, key, data, part_size=None):
        self.objects[key] = (data, s3_etag(data, part_size), part_size)
        return self.objects[key][1]

    def head_object(self, Bucket, Key, PartNumber=None):
        self.calls.append(("head", Key, PartNumber))
        data, etag, ps = self.objects[Key]
        if PartNumber is not None:
            if self.no_part_head or ps is None:
                raise RuntimeError("part head unavailable")
            return {"ContentLength": min(ps, len(data)), "ETag": f'"{etag}"'}
        return {"ContentLength": len(data), "ETag": f'"{etag}"',
                "LastModified": "2026-06-22T00:00:00+00:00"}

    def get_object(self, Bucket, Key):
        self.calls.append(("get", Key, None))
        data = self.objects[Key][0]
        if Key in self.corrupt:
            data = b"X" + data[1:]
        return {"Body": Body(data, self.truncate.get(Key))}


class DownloadCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.data = pathlib.Path(self.tmp.name) / "data_cache"
        self.s3 = FakeObjects()
        self.lines = []

    def tearDown(self):
        self.tmp.cleanup()

    def out(self, *a):
        self.lines.append(" ".join(str(x) for x in a))

    def unit(self, key, data, part_size=None, state="missing", kind="file",
             catalog_size=None, catalog_etag=None):
        etag = self.s3.put(key, data, part_size)
        return {"bucket": SC, "unit": key, "kind": kind,
                "size": catalog_size if catalog_size is not None else len(data),
                "etag": catalog_etag if catalog_etag is not None else etag,
                "state": state, "dest": self.data / SC / key}

    def run_plan(self, units, **kw):
        kw.setdefault("yes", True)
        plan = {"units": units}
        return download.run(plan, self.data, client=self.s3, out=self.out, **kw)

    def ledger(self):
        with open(self.data / "_ledger.csv", newline="") as fh:
            return list(csv.DictReader(fh))


class TestHappyPath(DownloadCase):

    def test_single_part_md5_verified_and_ledgered(self):
        data = b"cell_label,donor\n" * 1000
        u = self.unit("LEC/RNAseq/meta.csv", data)
        self.assertEqual(self.run_plan([u]), 0)
        self.assertEqual(u["dest"].read_bytes(), data)
        self.assertFalse(u["dest"].with_name("meta.csv.part").exists())
        row = self.ledger()[0]
        self.assertEqual(row["verified"], "md5 match")
        self.assertEqual(row["sha256"], hashlib.sha256(data).hexdigest())
        self.assertEqual(int(row["size"]), len(data))

    def test_multipart_etag_recomputed_exactly(self):
        data = bytes(range(256)) * 41          # 10,496 bytes
        u = self.unit("MTG/x.h5ad", data, part_size=1000)   # 11 parts
        self.assertTrue(u["etag"].endswith("-11"))
        self.assertEqual(self.run_plan([u]), 0)
        self.assertEqual(self.ledger()[0]["verified"], "multipart etag match")

    def test_multipart_part_boundaries_not_aligned_to_read_chunks(self):
        saved = download.CHUNK
        download.CHUNK = 333                    # forces parts to straddle reads
        try:
            data = b"abcdefghij" * 997
            u = self.unit("MTG/y.h5ad", data, part_size=1024)
            self.assertEqual(self.run_plan([u]), 0)
        finally:
            download.CHUNK = saved

    def test_size_only_when_part_size_cannot_be_learned(self):
        self.s3.no_part_head = True
        data = b"z" * 5000
        u = self.unit("MTG/z.h5ad", data, part_size=1000)
        self.assertEqual(self.run_plan([u]), 0)
        self.assertEqual(self.ledger()[0]["verified"], "size only")

    def test_ledger_appends(self):
        a = self.unit("a.csv", b"1" * 10)
        b = self.unit("b.csv", b"2" * 20)
        self.run_plan([a])
        self.run_plan([b])
        self.assertEqual([r["key"] for r in self.ledger()], ["a.csv", "b.csv"])


class TestFailuresLeaveNothingBehind(DownloadCase):

    def assert_clean(self, u):
        self.assertFalse(u["dest"].exists())
        self.assertFalse(u["dest"].with_name(u["dest"].name + ".part").exists())
        self.assertFalse((self.data / "_ledger.csv").exists())

    def test_changed_on_s3_since_catalog_fetches_nothing(self):
        u = self.unit("a.csv", b"new bytes", catalog_size=5)
        self.assertEqual(self.run_plan([u]), 7)
        self.assertNotIn("get", [c[0] for c in self.s3.calls])
        self.assertIn("changed on S3", "\n".join(self.lines))
        self.assert_clean(u)

    def test_changed_etag_same_size_fetches_nothing(self):
        u = self.unit("a.csv", b"abcdef", catalog_etag="deadbeef")
        self.assertEqual(self.run_plan([u]), 7)
        self.assertNotIn("get", [c[0] for c in self.s3.calls])
        self.assert_clean(u)

    def test_corrupted_bytes_are_deleted(self):
        u = self.unit("a.csv", b"hello world" * 100)
        self.s3.corrupt.add("a.csv")
        self.assertEqual(self.run_plan([u]), 7)
        self.assertIn("checksum mismatch", "\n".join(self.lines))
        self.assert_clean(u)

    def test_truncated_stream_is_deleted(self):
        u = self.unit("a.csv", b"x" * 1000)
        self.s3.truncate["a.csv"] = 400
        self.assertEqual(self.run_plan([u]), 7)
        self.assertIn("received 400", "\n".join(self.lines))
        self.assert_clean(u)

    def test_failure_stops_the_batch_and_reports_progress(self):
        a = self.unit("a.csv", b"a" * 10)
        b = self.unit("b.csv", b"b" * 10, catalog_size=3)
        c = self.unit("c.csv", b"c" * 10)
        self.assertEqual(self.run_plan([a, b, c]), 7)
        self.assertTrue(a["dest"].exists())
        self.assertFalse(c["dest"].exists())
        self.assertIn("1 of 3", "\n".join(self.lines))


class TestLedgerPaths(DownloadCase):
    """Ledger paths are relative to data_cache/, so moving the project keeps them valid."""

    def test_written_relative_to_data_cache(self):
        u = self.unit("LEC/RNAseq/meta.csv", b"x" * 10)
        self.assertEqual(self.run_plan([u]), 0)
        self.assertEqual(self.ledger()[0]["path"], f"{SC}/LEC/RNAseq/meta.csv")
        self.assertNotIn(str(self.data), (self.data / "_ledger.csv").read_text())
        self.assertEqual(download.read_ledger(self.data)[0]["file"], u["dest"])

    def old_row(self, path):
        download.append_ledger(self.data / "_ledger.csv", {
            "downloaded_at": "2026-11-20", "bucket": SC, "key": "a.csv", "size": 1,
            "etag": "e", "s3_last_modified": "", "sha256": "0" * 64,
            "verified": "md5 match", "path": str(path)})
        return download.read_ledger(self.data)[0]["file"]

    def test_old_absolute_row_still_read(self):
        f = self.data / SC / "a.csv"
        f.parent.mkdir(parents=True)
        f.write_bytes(b"x")
        self.assertEqual(self.old_row(f.resolve()), f.resolve())

    def test_old_absolute_row_found_after_the_folder_moved(self):
        f = self.data / SC / "a.csv"
        f.parent.mkdir(parents=True)
        f.write_bytes(b"x")
        found = self.old_row(pathlib.Path("/no/such/old-folder/data_cache") / SC / "a.csv")
        self.assertEqual(found, self.data / SC / "a.csv")

    def test_missing_unrelated_path_is_left_as_is(self):
        self.assertEqual(self.old_row("/no/such/place/a.csv"),
                         pathlib.Path("/no/such/place/a.csv"))

    def test_no_ledger_means_no_rows(self):
        self.assertEqual(download.read_ledger(self.data / "nope"), [])


class TestProgressMeter(DownloadCase):
    """The percentage meter is for a person at a terminal; tests and logs stay clean."""

    def stderr_of(self, **kw):
        import contextlib
        u = self.unit("a.csv", b"x" * 5000)
        class TerminalStderr(io.StringIO):     # a stderr that says a person is watching
            def isatty(self):
                return True
        err = TerminalStderr()
        with contextlib.redirect_stderr(err):
            self.assertEqual(self.run_plan([u], **kw), 0)
        return err.getvalue()

    def test_quiet_when_output_is_redirected(self):
        self.assertEqual(self.stderr_of(), "")

    def test_shown_when_asked(self):
        self.assertIn("100.0%", self.stderr_of(show_progress=True))


class TestGuards(DownloadCase):

    def test_no_terminal_and_no_yes_refuses(self):
        u = self.unit("a.csv", b"x" * 10)
        code = self.run_plan([u], yes=False, isatty=False)
        self.assertEqual(code, 6)
        self.assertEqual(self.s3.calls, [])

    def test_confirmation_no_aborts(self):
        u = self.unit("a.csv", b"x" * 10)
        code = self.run_plan([u], yes=False, isatty=True, ask=lambda _: "n")
        self.assertEqual(code, 1)
        self.assertEqual(self.s3.calls, [])

    def test_confirmation_yes_proceeds(self):
        u = self.unit("a.csv", b"x" * 10)
        code = self.run_plan([u], yes=False, isatty=True, ask=lambda _: "y")
        self.assertEqual(code, 0)

    def test_default_asks_even_for_tiny_downloads(self):
        asked = []
        u = self.unit("a.csv", b"x")
        self.run_plan([u], yes=False, isatty=True,
                      ask=lambda q: asked.append(q) or "n")
        self.assertEqual(len(asked), 1)

    def test_present_files_are_skipped_unless_forced(self):
        u = self.unit("a.csv", b"x" * 10, state="present")
        self.assertEqual(self.run_plan([u]), 0)
        self.assertEqual(self.s3.calls, [])
        self.assertEqual(self.run_plan([u], force=True), 0)
        self.assertIn("get", [c[0] for c in self.s3.calls])

    def test_collapsed_units_are_refused(self):
        u = self.unit("MTG/RNAseq/pseudoprogression-plots/", b"", kind="collapsed")
        self.assertEqual(self.run_plan([u]), 4)
        self.assertEqual(self.s3.calls, [])

    def test_unknown_bucket_is_refused(self):
        u = self.unit("a.csv", b"x")
        u["bucket"] = "someone-elses-bucket"
        self.assertEqual(self.run_plan([u]), 4)
        self.assertEqual(self.s3.calls, [])

    def test_not_enough_disk_space_is_refused(self):
        u = self.unit("a.csv", b"x" * 100)

        class Usage:
            free = 50
        code = self.run_plan([u], disk_usage=lambda p: Usage)
        self.assertEqual(code, 5)
        self.assertEqual(self.s3.calls, [])


class TestETagHasher(unittest.TestCase):

    def test_matches_s3_formula(self):
        data = b"q" * 12345
        for ps in (None, 1000, 4096, 12345, 20000):
            h = download.ETagHasher(ps)
            h.update(data[:5000])
            h.update(data[5000:])
            self.assertEqual(h.hexdigest(), s3_etag(data, ps), ps)

    def test_multipart_count(self):
        self.assertEqual(download.multipart_count("abc-95"), 95)
        self.assertIsNone(download.multipart_count("0123abc"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
