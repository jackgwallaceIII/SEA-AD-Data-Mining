"""
test_report.py — unit tests for scripts/report.py (per-region data-structure reports).

No network, no data: a temporary catalog is built over the in-memory FakeS3
bucket, and hand-written inspect-report JSONs stand in for files that were
inspected. One of them carries a category level that must NOT be copied
into the regional report.

    python -m unittest discover -s tests -t tests -v

Why these tests: the regional report is what the next person reads first.
It must keep region-only and shared units apart, say how many units of a
kind were actually checked, list kinds never checked, attach files inside
folder units to the right folder, carry the known issues, repeat no
category level, and give identical output when nothing changed.
"""

import json
import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import index  # noqa: E402
import report  # noqa: E402
from test_index import SC, SC_KEYS  # noqa: E402
from test_issues import issue, write_issues  # noqa: E402
from test_sources import FakeS3  # noqa: E402

H5AD = "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad"
META = "MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei_metadata.2026-06-22.csv"


def h5ad_report():
    return {"bucket": SC, "key": H5AD, "status": "inspected (AnnData)",
            "ledger": "matches ledger (md5 match, 2026-11-20T10:00:00+00:00)",
            "anndata": {"n_obs": 1234, "n_vars": 36601,
                        "var": {"index": {"looks_like": "gene symbols (not Ensembl IDs)"}},
                        "matrices": [
                            {"name": "X", "storage": "csr sparse",
                             "values": "non-whole numbers (sample of 9): transformed"},
                            {"name": "layers/UMIs", "storage": "csr sparse",
                             "values": "whole numbers only (sample of 9): looks like raw"}]}}


def csv_report():
    return {"bucket": SC, "key": META, "status": "inspected", "rows": 50,
            "columns": [{"name": "Braak", "levels": ["SECRETLEVEL"], "levels_withheld": None,
                         "empty": False},
                        {"name": "x", "levels": None, "levels_withheld": "id", "empty": True}]}


def plot_report():
    """A file INSIDE the collapsed plots folder."""
    return {"bucket": SC, "key": "MTG/RNAseq/pseudoprogression-plots/a/overview.jpg",
            "status": "detected binary; no structural reader for it yet"}


class ReportCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self.tmp.name)
        self.cat, self.ins, self.out = base / "catalog", base / "inspect", base / "reports"
        index.refresh(SC, client=FakeS3(SC_KEYS), catalog_dir=self.cat, out=lambda *a: None)
        for rep in (h5ad_report(), csv_report(), plot_report()):
            p = self.ins / rep["bucket"] / (rep["key"] + ".json")
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps(rep))
        write_issues(self.cat / "known_issues.csv",
                     [issue(id="KI-900", applies_to="MTG/RNAseq/*.h5ad")])
        self.lines = []

    def tearDown(self):
        self.tmp.cleanup()

    def run_report(self, *a, **kw):
        return report.run(*a, catalog_dir=self.cat, inspect_dir=self.ins, report_dir=self.out,
                          out=self.lines.append, **kw)

    def text(self, code="MTG"):
        return (self.out / f"{code}.md").read_text()


class TestRegionPage(ReportCase):

    def setUp(self):
        super().setUp()
        self.run_report("MTG")
        self.md = self.text()

    def test_sections_present(self):
        for heading in ("# MTG — middle temporal gyrus", "## At a glance", "### MTG only",
                        "### Shared with other regions", "## What was checked inside the files",
                        "## Known upstream issues", "## How to get it",
                        "## How to read this report"):
            self.assertIn(heading, self.md)

    def test_verdicts_in_plain_words(self):
        self.assertIn("1,234 cells × 36,601 genes; X: normalised/transformed; "
                      "layers/UMIs: raw counts", self.md)
        self.assertIn("downloaded, 2026-11-20", self.md)

    def test_coverage_and_not_yet_inspected(self):
        self.assertIn("1 of 1", self.md)              # the region's one dataset was checked
        self.assertIn("**Not yet inspected**", self.md)
        not_yet = self.md.split("**Not yet inspected**")[1].split("\n")[0]
        self.assertNotIn("cell metadata table", not_yet)

    def test_file_inside_a_folder_counts_for_that_folder(self):
        self.assertIn("`a/overview.jpg`", self.md)
        not_yet = self.md.split("**Not yet inspected**")[1].split("\n")[0]
        self.assertNotIn("plots", not_yet)

    def test_shared_units_kept_apart(self):
        only = self.md.split("### MTG only")[1].split("### Shared")[0]
        self.assertNotIn("per-subclass object", only)
        shared = self.md.split("### Shared with other regions")[1].split("## What was")[0]
        self.assertIn("per-subclass object", shared)

    def test_known_issue_carried(self):
        self.assertIn("KI-900", self.md)

    def test_no_category_level_repeated(self):
        self.assertNotIn("SECRETLEVEL", self.md)

    def test_links_point_at_the_per_file_reports(self):
        self.assertIn("(../inspect/" + SC + "/" + H5AD + ".md)", self.md)

    def test_example_command_uses_an_inspected_file(self):
        self.assertIn(f'get --match "{META}" --dry-run', self.md)


class TestRunModes(ReportCase):

    def test_same_inputs_same_output(self):
        self.run_report("MTG")
        first = self.text()
        self.run_report("MTG")
        self.assertEqual(first, self.text())

    def test_all_regions_and_index(self):
        summaries = self.run_report(all_regions=True)
        idx = (self.out / "README.md").read_text()
        self.assertIn("[MTG](MTG.md)", idx)
        self.assertIn("[CaH](CaH.md)", idx)
        self.assertEqual({s["code"] for s in summaries},
                         {p.stem for p in self.out.glob("*.md")} - {"README"})

    def test_index_totals_add_up_to_the_catalog(self):
        import browse
        self.run_report(all_regions=True)
        idx = (self.out / "README.md").read_text()
        rows, hidden = browse.filter_units(browse.load_units(self.cat), None, None, None, False)
        total = next(ln for ln in idx.splitlines() if ln.startswith("| **total current data**"))
        self.assertIn(f"| {len(rows):,} |", total)
        parts = {name: int(ln.split("|")[3].strip().replace(",", ""))
                 for ln in idx.splitlines() for name in ("all regions, region-only",
                                                         "shared by several regions",
                                                         "no single region")
                 if ln.startswith(f"| **{name}**")}
        self.assertEqual(sum(parts.values()), len(rows))

    def test_shared_inspections_not_counted_as_region_only(self):
        summaries = {s["code"]: s for s in self.run_report(all_regions=True)}
        mtg = summaries["MTG"]
        self.assertEqual(mtg["inspected_only"] + mtg["inspected_shared"], mtg["inspected"])
        self.assertEqual(mtg["inspected_shared"], 0)    # all three fixtures are MTG-only

    def test_unknown_region_refused(self):
        with self.assertRaises(SystemExit):
            self.run_report("Atlantis")

    def test_region_with_nothing_inspected_says_so(self):
        self.run_report("CaH")
        self.assertIn("No file of this region has been inspected yet.", self.text("CaH"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
