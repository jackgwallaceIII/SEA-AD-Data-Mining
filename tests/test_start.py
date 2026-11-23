"""
test_start.py — unit tests for scripts/start.py (the front door).

No network: a temporary catalog over the in-memory FakeS3 bucket, and a
fake package finder to simulate a broken environment.

    python -m unittest discover -s tests -t tests -v

Why these tests: `start` is the first command a newcomer runs. It must
name a missing package and say how to fix it, recognise the conda 'base'
mistake, cope with no catalog at all, show the five workflow steps with
real commands for the region asked about, and contact nothing.
"""

import contextlib
import io
import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import index  # noqa: E402
import start  # noqa: E402
from test_index import SC, SC_KEYS  # noqa: E402
from test_sources import FakeS3  # noqa: E402


def env(missing=(), conda="seaad", version=(3, 12, 1)):
    return start.environment(find_spec=lambda m: None if m in missing else object(),
                             version=version, env={"CONDA_DEFAULT_ENV": conda})


class StartCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self.tmp.name)
        self.cat, self.data, self.rep = base / "catalog", base / "data_cache", base / "reports"
        index.refresh(SC, client=FakeS3(SC_KEYS), catalog_dir=self.cat, out=lambda *a: None)
        self.lines = []

    def tearDown(self):
        self.tmp.cleanup()

    def run_start(self, region=None, **kw):
        kw.setdefault("env_check", env())
        code = start.run(region, catalog_dir=self.cat, data_dir=self.data, report_dir=self.rep,
                         out=self.lines.append, **kw)
        return code, "\n".join(self.lines)


class TestSetup(StartCase):

    def test_all_good(self):
        code, text = self.run_start()
        self.assertEqual(code, 0)
        self.assertIn("ok  boto3", text)

    def test_missing_package_named_with_the_fix(self):
        code, text = self.run_start(env_check=env(missing={"boto3"}))
        self.assertEqual(code, 3)
        self.assertIn("NO  boto3", text)
        self.assertIn("conda activate seaad", text)
        self.assertIn("pip install -r requirements-lock.txt", text)

    def test_conda_base_is_called_out(self):
        _, text = self.run_start(env_check=env(missing={"h5py"}, conda="base"))
        self.assertIn("You are in conda 'base'", text)

    def test_old_python(self):
        code, text = self.run_start(env_check=env(version=(3, 8, 0)))
        self.assertEqual(code, 3)
        self.assertIn("needs 3.10 or newer", text)

    def test_optional_package_does_not_fail(self):
        code, text = self.run_start(env_check=env(missing={"pyarrow"}))
        self.assertEqual(code, 0)
        self.assertIn("--  pyarrow", text)


class TestDownloadCount(StartCase):

    def test_relative_and_old_absolute_ledger_rows_both_counted(self):
        import download
        a, b = self.data / SC / "a.csv", self.data / SC / "b.csv"
        a.parent.mkdir(parents=True)
        a.write_bytes(b"x" * 10)
        b.write_bytes(b"y" * 20)
        for path, size in ((f"{SC}/a.csv", 10), (str(b.resolve()), 20)):
            download.append_ledger(self.data / "_ledger.csv", {
                "downloaded_at": "2026-11-23", "bucket": SC, "key": "k", "size": size,
                "etag": "e", "s3_last_modified": "", "sha256": "0" * 64,
                "verified": "md5 match", "path": path})
        st = start.local_state(self.cat, self.data, self.rep)
        self.assertEqual((st["downloads"], st["download_size"]), (2, 30))


class TestWorkflow(StartCase):

    def test_five_steps_with_region_commands(self):
        _, text = self.run_start("MTG")
        for step in ("1. LOOK", "2. READ", "3. PREVIEW", "4. DOWNLOAD", "5. DESCRIBE"):
            self.assertIn(step, text)
        self.assertIn("python scripts/index.py report MTG", text)
        self.assertIn("--dry-run", text)
        self.assertIn("conda activate seaad", text)

    def test_example_is_a_small_real_file_of_the_region(self):
        _, text = self.run_start("MTG")
        self.assertIn('get --match "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei_metadata.2026-06-22.csv"',
                      text)

    def test_commands_have_no_comments(self):
        _, text = self.run_start("MTG")
        cmds = [ln.strip() for ln in text.splitlines() if ln.strip().startswith("python scripts")]
        self.assertTrue(cmds)
        self.assertFalse([c for c in cmds if "#" in c])

    def test_unknown_region_lists_the_regions(self):
        code, text = self.run_start("Atlantis")
        self.assertEqual(code, 2)
        self.assertIn("LEC (", text)

    def test_promises_stated(self):
        _, text = self.run_start()
        self.assertIn("Nothing downloads without a size preview and your yes", text)
        self.assertIn("never downloaded", text)


class TestNoCatalog(unittest.TestCase):

    def test_explains_how_to_build_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            lines = []
            base = pathlib.Path(tmp)
            code = start.run(None, catalog_dir=base / "catalog", data_dir=base / "d",
                             report_dir=base / "r", env_check=env(), out=lines.append)
        text = "\n".join(lines)
        self.assertEqual(code, 0)
        self.assertIn("catalog:   none yet", text)
        self.assertIn("index.py refresh --bucket", text)


class TestEntryPoint(unittest.TestCase):

    def test_no_command_runs_start(self):
        calls = []
        saved = start.run
        start.run = lambda *a, **k: calls.append(a) or 0
        try:
            self.assertEqual(index.main([]), 0)
        finally:
            start.run = saved
        self.assertEqual(calls, [()])

    def test_start_is_listed_first_in_help(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), self.assertRaises(SystemExit):
            index.main(["--help"])
        text = " ".join(buf.getvalue().split())          # argparse wraps lines
        self.assertIn("python scripts/index.py start", text)
        self.assertIn("START HERE", text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
