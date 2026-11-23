"""
test_issues.py — unit tests for scripts/issues.py and where issues are shown.

No network: a temporary catalog is built over the in-memory FakeS3 buckets
and a temporary known_issues.csv is written next to it.

    python -m unittest discover -s tests -t tests -v

Why these tests: a known issue is only useful if it reaches the person
using the file. These check that the CSV is validated (a typo must not
silently hide an issue), that patterns match only what they claim, and
that show, search, get and inspect each print the issues for their units.
The last test checks the real catalog/known_issues.csv is valid.
"""

import contextlib
import csv
import io
import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import browse  # noqa: E402
import index  # noqa: E402
import issues  # noqa: E402
import plan  # noqa: E402
from test_index import SC, SC_KEYS  # noqa: E402
from test_sources import FakeS3  # noqa: E402

MTG_H5AD = "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad"


def issue(**kw):
    base = dict(id="KI-900", bucket=SC, applies_to="MTG/RNAseq/*final-nuclei*.h5ad",
                severity="caution", summary="synthetic test issue", detail="d",
                evidence="e", found="2026-11-22", status="open")
    return {**base, **kw}


def write_issues(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=issues.FIELDS)
        w.writeheader()
        w.writerows(rows)


class IssueCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self.tmp.name)
        self.cat, self.data = base / "catalog", base / "data_cache"
        index.refresh(SC, client=FakeS3(SC_KEYS), catalog_dir=self.cat, out=lambda *a: None)
        self.file = self.cat / "known_issues.csv"
        self.lines = []

    def tearDown(self):
        self.tmp.cleanup()

    def out(self, *a):
        self.lines.append(" ".join(str(x) for x in a))

    def text(self):
        return "\n".join(self.lines)


class TestValidation(IssueCase):

    def test_good_file_loads(self):
        write_issues(self.file, [issue()])
        self.assertEqual([i["id"] for i in issues.load(self.file)], ["KI-900"])

    def test_bad_rows_are_refused_with_the_line(self):
        for bad in (dict(severity="urgent"), dict(bucket="some-other-bucket"),
                    dict(found="22/11/2026"), dict(status="maybe"), dict(applies_to="")):
            write_issues(self.file, [issue(**bad)])
            with self.assertRaises(issues.IssueError, msg=str(bad)) as cm:
                issues.load(self.file)
            self.assertIn("line 2", str(cm.exception))

    def test_web_source_accepted(self):
        write_issues(self.file, [issue(bucket="web", applies_to="cdn.example.org/*.xlsx")])
        found = issues.load(self.file)
        self.assertTrue(issues.applies(found[0], "web", "cdn.example.org/a.xlsx"))

    def test_duplicate_id_refused(self):
        write_issues(self.file, [issue(), issue()])
        with self.assertRaises(issues.IssueError):
            issues.load(self.file)

    def test_resolved_hidden_unless_asked(self):
        write_issues(self.file, [issue(status="resolved")])
        self.assertEqual(issues.load(self.file), [])
        self.assertEqual(len(issues.load(self.file, include_resolved=True)), 1)

    def test_missing_file_means_no_issues(self):
        self.assertEqual(issues.load(self.cat / "nope.csv"), [])


class TestMatching(unittest.TestCase):

    def test_glob_ignores_case_and_crosses_folders(self):
        i = issue(applies_to="mtg/*.h5ad")
        self.assertTrue(issues.applies(i, SC, MTG_H5AD))
        self.assertFalse(issues.applies(i, SC, "LEC/RNAseq/x.h5ad"))
        self.assertFalse(issues.applies(i, "sea-ad-spatial-transcriptomics", MTG_H5AD))

    def test_several_patterns(self):
        i = issue(applies_to="LEC/*.csv; MTG/*.h5ad")
        self.assertTrue(issues.applies(i, SC, MTG_H5AD))
        self.assertTrue(issues.applies(i, SC, "LEC/RNAseq/m.csv"))
        self.assertFalse(issues.applies(i, SC, "LEC/RNAseq/m.h5ad"))

    def test_most_severe_first(self):
        rows = [{"bucket": SC, "unit": MTG_H5AD}]
        found = issues.matching([issue(id="A", severity="note"),
                                 issue(id="B", severity="problem")], rows)
        self.assertEqual([i["id"] for i, _ in found], ["B", "A"])


class TestShownWhereUsed(IssueCase):

    def test_show_region(self):
        write_issues(self.file, [issue()])
        r = browse.show("MTG", catalog_dir=self.cat, out=self.out)
        self.assertEqual(len(r["issues"]), 1)
        self.assertIn("KI-900", self.text())

    def test_show_other_region_stays_quiet(self):
        write_issues(self.file, [issue()])
        r = browse.show("LEC", catalog_dir=self.cat, out=self.out)
        self.assertEqual(r["issues"], [])
        self.assertNotIn("KI-900", self.text())

    def test_overview_counts_them(self):
        write_issues(self.file, [issue()])
        browse.show(catalog_dir=self.cat, out=self.out)
        self.assertIn("1 known upstream issue", self.text())

    def test_search(self):
        write_issues(self.file, [issue()])
        r = browse.search("final-nuclei", catalog_dir=self.cat, out=self.out)
        self.assertEqual(len(r["issues"]), 1)

    def test_get_preview_warns_but_does_not_block(self):
        write_issues(self.file, [issue(severity="problem")])
        p, code = plan.preview(match=MTG_H5AD, catalog_dir=self.cat, data_dir=self.data,
                               out=self.out)
        self.assertEqual(code, 0)
        self.assertEqual(len(p["issues"]), 1)
        self.assertIn("known issue KI-900 (problem)", self.text())

    def test_inspect_report_carries_them(self):
        import inspect_data
        write_issues(self.file, [issue(applies_to="MTG/RNAseq/*.csv")])
        f = self.data / SC / "MTG/RNAseq/meta.csv"
        f.parent.mkdir(parents=True)
        f.write_text("a,b\n1,2\n")
        r, md = inspect_data.inspect(f, data_dir=self.data, out_dir=self.data / "out",
                                     issues_path=self.file)
        self.assertEqual([i["id"] for i in r["known_issues"]], ["KI-900"])
        self.assertIn("## Known upstream issues for this file", md)

    def test_issues_command_lists_everything(self):
        write_issues(self.file, [issue(), issue(id="KI-901", status="resolved")])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            shown = issues.report(path=self.file)
        self.assertEqual([i["id"] for i in shown], ["KI-900"])
        self.assertIn("synthetic test issue", buf.getvalue())


class TestRealFile(unittest.TestCase):

    def test_committed_issues_file_is_valid(self):
        found = issues.load(include_resolved=True)
        self.assertTrue(found, "catalog/known_issues.csv is missing or empty")
        for i in found:
            self.assertTrue(i["evidence"], f"{i['id']} has no evidence")


if __name__ == "__main__":
    unittest.main(verbosity=2)
