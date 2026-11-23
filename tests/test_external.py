"""
test_external.py — unit tests for scripts/external.py (web files, controlled studies)
and scripts/xlsx.py (Excel workbooks read with the standard library).

No network: a fake URL opener answers HEAD and GET from memory. The
workbook is assembled from its XML parts with zipfile, so every cell is
known; donor IDs and free text are planted and must not reach a report.

    python -m unittest discover -s tests -t tests -v

Why these tests: controlled-access data must be impossible to download;
open web files must get the same preview / yes / size check / ledger as
S3 files; and an Excel sheet must be described by the same column rules
as a CSV, including counting blanks in columns that start partway down.
"""

import csv
import hashlib
import io
import json
import pathlib
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import external  # noqa: E402
import inspect_data as ins  # noqa: E402
import xlsx  # noqa: E402

URL = "https://cdn.example.org/files/abc_sea-ad_cohort_donor_metadata.xlsx"


class FakeResponse:
    def __init__(self, body, headers):
        self.body, self.headers, self.pos = body, headers, 0

    def read(self, n=-1):
        out = self.body[self.pos:] if n < 0 else self.body[self.pos:self.pos + n]
        self.pos += len(out)
        return out

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class FakeWeb:
    """urlopen stand-in: {url: bytes}; records every request method."""

    def __init__(self, files, etag="md5", announce=None):
        self.files, self.etag, self.announce, self.calls = files, etag, announce, []

    def __call__(self, req, timeout=None):
        body = self.files[req.full_url]
        self.calls.append((req.get_method(), req.full_url, req.headers.get("User-agent")))
        tag = hashlib.md5(body).hexdigest() if self.etag == "md5" else self.etag
        headers = {"Content-Length": str(self.announce or len(body)), "ETag": f'"{tag}"',
                   "Last-Modified": "Tue, 01 Jul 2025 00:00:00 GMT",
                   "Content-Type": "application/vnd.openxmlformats"}
        return FakeResponse(b"" if req.get_method() == "HEAD" else body, headers)


def sources_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=external.FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({**{k: "" for k in external.FIELDS}, **r})


OPEN = {"id": "web-donor-metadata", "title": "donors", "access": "open", "url": URL,
        "format": "xlsx", "modality": "donor metadata", "listed_on": "https://example.org"}
CTRL = {"id": "syn26223298", "title": "raw", "access": "controlled",
        "url": "https://adknowledgeportal.synapse.org/x?Study=syn26223298"}


class SourceCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = pathlib.Path(self.tmp.name)
        self.src, self.status = self.base / "sources.csv", self.base / "status.csv"
        sources_csv(self.src, [OPEN, CTRL])
        self.data = self.base / "data_cache"
        self.lines = []

    def tearDown(self):
        self.tmp.cleanup()


class TestList(SourceCase):

    def test_real_list_is_valid_and_has_the_confirmed_studies(self):
        rows = {r["id"]: r for r in external.load()}
        for syn in ("syn26223298", "syn64410371"):
            self.assertEqual(rows[syn]["access"], "controlled")
        self.assertTrue(all(r["url"].startswith("https://") for r in rows.values()
                            if r["access"] == "open"))
        self.assertEqual(sum(r["access"] == "open" for r in rows.values()), 7)

    def test_bad_rows_refused(self):
        for bad in ({"access": "public"}, {"url": "ftp://x"},
                    {"url": "http://cdn.example.org/a.xlsx"}):
            sources_csv(self.src, [{**OPEN, **bad}])
            with self.assertRaises(external.SourceError, msg=str(bad)):
                external.load(self.src)
        sources_csv(self.src, [{**CTRL, "id": "syn12ab"}])
        with self.assertRaises(external.SourceError):
            external.load(self.src)
        sources_csv(self.src, [OPEN, OPEN])
        with self.assertRaises(external.SourceError):
            external.load(self.src)


class TestCheck(SourceCase):

    def test_head_only_and_status_written(self):
        web = FakeWeb({URL: b"x" * 1234})
        code = external.check(self.src, self.status, opener=web, out=self.lines.append)
        self.assertEqual(code, 0)
        self.assertEqual([m for m, _, _ in web.calls], ["HEAD"])     # never GET
        self.assertEqual(web.calls[0][2], external.USER_AGENT)
        st = external.load_status(self.status)["web-donor-metadata"]
        self.assertEqual(st["size"], "1234")

    def test_change_is_reported(self):
        external.check(self.src, self.status, opener=FakeWeb({URL: b"x" * 10}), out=lambda *a: 0)
        external.check(self.src, self.status, opener=FakeWeb({URL: b"y" * 20}),
                       out=self.lines.append)
        self.assertIn("CHANGED since last check", "\n".join(self.lines))


class TestGet(SourceCase):

    def get(self, pattern="web-donor-metadata", **kw):
        kw.setdefault("yes", True)
        return external.run_get(pattern, data_dir=self.data, sources_path=self.src,
                                out=self.lines.append, **kw)

    def test_controlled_is_always_refused(self):
        for kw in ({}, {"yes": True, "force": True}, {"dry_run": True}):
            with self.assertRaises(SystemExit) as cm:
                self.get("syn26223298", opener=FakeWeb({}), **kw)
            self.assertIn("never downloads", str(cm.exception))
            self.assertIn("adknowledgeportal", str(cm.exception))

    def test_glob_including_controlled_is_refused(self):
        with self.assertRaises(SystemExit):
            self.get("*", opener=FakeWeb({URL: b"x"}))

    def test_dry_run_only_heads(self):
        web = FakeWeb({URL: b"x" * 100})
        self.assertEqual(self.get(dry_run=True, opener=web), 0)
        self.assertEqual({m for m, _, _ in web.calls}, {"HEAD"})
        self.assertFalse(self.data.exists())

    def test_download_verified_and_ledgered(self):
        body = b"PK\x03\x04" + b"z" * 5000
        self.assertEqual(self.get(opener=FakeWeb({URL: body})), 0)
        dest = self.data / "web" / "cdn.example.org" / "abc_sea-ad_cohort_donor_metadata.xlsx"
        self.assertEqual(dest.read_bytes(), body)
        led = list(csv.DictReader((self.data / "_ledger.csv").read_text().splitlines()))
        self.assertEqual(led[0]["verified"], "md5 matches ETag")
        self.assertEqual(led[0]["sha256"], hashlib.sha256(body).hexdigest())
        self.assertEqual(led[0]["bucket"], "web")
        self.assertEqual(led[0]["path"],                         # relative to data_cache/
                         "web/cdn.example.org/abc_sea-ad_cohort_donor_metadata.xlsx")

    def test_etag_that_is_not_md5_says_size_only(self):
        self.get(opener=FakeWeb({URL: b"abc"}, etag="W/opaque-123"))
        led = list(csv.DictReader((self.data / "_ledger.csv").read_text().splitlines()))
        self.assertEqual(led[0]["verified"], "size only")

    def test_short_transfer_deletes_the_part_file(self):
        code = self.get(opener=FakeWeb({URL: b"abc"}, announce=999))
        self.assertEqual(code, 7)
        self.assertEqual(list((self.data / "web").rglob("*")),
                         [self.data / "web" / "cdn.example.org"])

    def test_needs_a_yes(self):
        code = self.get(yes=False, isatty=False, opener=FakeWeb({URL: b"abc"}))
        self.assertEqual(code, 6)
        self.assertFalse((self.data / "web").exists())


# ---- xlsx ----------------------------------------------------------------------

M = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


def cell(ref, value, shared):
    if isinstance(value, bool):
        return f'<c r="{ref}" t="b"><v>{int(value)}</v></c>'
    if isinstance(value, (int, float)):
        return f'<c r="{ref}"><v>{value}</v></c>'
    if value not in shared:
        shared.append(value)
    return f'<c r="{ref}" t="s"><v>{shared.index(value)}</v></c>'


def workbook(sheets):
    """sheets: [(name, hidden, rows)] with rows as lists (None = empty cell)."""
    shared, buf = [], io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        wb = (f'<workbook xmlns="{M}" xmlns:r="{R}"><sheets>'
              + "".join(f'<sheet name="{n}" sheetId="{i + 1}" r:id="rId{i + 1}"'
                        + (' state="hidden"' if h else "") + "/>"
                        for i, (n, h, _) in enumerate(sheets)) + "</sheets></workbook>")
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/'
                   'relationships">' + "".join(
                       f'<Relationship Id="rId{i + 1}" Target="worksheets/sheet{i + 1}.xml" '
                       f'Type="x"/>' for i in range(len(sheets))) + "</Relationships>")
        for i, (_, _, rows) in enumerate(sheets):
            xml = "".join(f'<row r="{r + 1}">' + "".join(
                cell(f"{chr(65 + j)}{r + 1}", v, shared) for j, v in enumerate(row)
                if v is not None) + "</row>" for r, row in enumerate(rows))
            z.writestr(f"xl/worksheets/sheet{i + 1}.xml",
                       f'<worksheet xmlns="{M}"><sheetData>{xml}</sheetData></worksheet>')
        z.writestr("xl/sharedStrings.xml", f'<sst xmlns="{M}">' + "".join(
            f"<si><t>{s}</t></si>" for s in shared) + "</sst>")
    return buf.getvalue()


DONORS = [["Donor ID", "Sex", "Age at Death", "Cognitive Status", "Note"]] + [
    [f"H20.33.{i:03d}" if i != 3 else "H21.33.999", ["Female", "Male"][i % 2], 70 + i,
     ["Dementia", "No dementia"][i % 2], "SECRET free text" if i >= 4 else None]
    for i in range(8)]


class TestWorkbook(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.data = pathlib.Path(self.tmp.name) / "data_cache"
        self.path = self.data / "web" / "cdn.example.org" / "donors.xlsx"
        self.path.parent.mkdir(parents=True)

    def tearDown(self):
        self.tmp.cleanup()

    def inspect(self, data):
        self.path.write_bytes(data)
        return ins.inspect(self.path, data_dir=self.data,
                           out_dir=pathlib.Path(self.tmp.name) / "out")

    def col(self, sheet, name):
        return next(c for c in sheet["columns"] if c["name"] == name)

    def test_sheet_rows_types_levels_and_withholding(self):
        r, md = self.inspect(workbook([("Donors", False, DONORS)]))
        sh = r["workbook"]["sheets"][0]
        self.assertEqual((sh["rows"], sh["header_row"], len(sh["columns"])), (8, 1, 5))
        self.assertEqual(self.col(sh, "Sex")["levels"], ["Female", "Male"])
        self.assertEqual(self.col(sh, "Age at Death")["type"], "integer")
        self.assertIsNone(self.col(sh, "Age at Death")["levels"])
        self.assertIsNotNone(self.col(sh, "Donor ID")["levels_withheld"])
        self.assertIn("Column names taken from row 1", md)
        blob = md + json.dumps(r)
        for secret in ("H21.33.999", "H20.33.001"):
            self.assertNotIn(secret, blob)

    def test_group_label_row_above_the_header_is_skipped(self):
        rows = [[None, "RIPA Buffer Tissue extractions", None, None, "GuHCl extractions"],
                ["Donor ID", "abeta40", "abeta42", "tTau", "abeta40 (GuHCl)"]] + [
               [f"H20.33.{i:03d}", 1.5 + i, 2.0 + i, 3 + i, 0.25 * i] for i in range(8)]
        r, md = self.inspect(workbook([("MTG_Luminex", False, rows)]))
        sh = r["workbook"]["sheets"][0]
        self.assertEqual((sh["header_row"], sh["rows_above_header"], sh["rows"]), (2, 1, 8))
        self.assertEqual([c["name"] for c in sh["columns"]][:3], ["Donor ID", "abeta40", "abeta42"])
        self.assertEqual(self.col(sh, "abeta40")["type"], "float")
        self.assertIn("1 title/group-label row(s) above it were skipped", md)

    def test_repeated_names_get_their_group_label_and_blank_rows_skip(self):
        rows = [[None, "RIPA", None, "GuHCl", None],
                ["Donor ID", "ABeta40 pg/ug", "tTAU pg/ug", "ABeta40 pg/ug", "tTAU pg/ug"]] + [
               [f"H20.33.{i:03d}", 1.5, 2.5, 3.5, 4.5] for i in range(6)] + [[]]
        r, _ = self.inspect(workbook([("Lx", False, rows)]))
        sh = r["workbook"]["sheets"][0]
        self.assertEqual([c["name"] for c in sh["columns"]],
                         ["Donor ID", "ABeta40 pg/ug [RIPA]", "tTAU pg/ug [RIPA]",
                          "ABeta40 pg/ug [GuHCl]", "tTAU pg/ug [GuHCl]"])
        self.assertEqual(sh["rows"], 6)
        self.assertEqual({c["missing"] for c in sh["columns"]}, {0})

    def test_source_line_is_the_web_address(self):
        _, md = self.inspect(workbook([("Donors", False, DONORS)]))
        self.assertIn("**Source:** `https://cdn.example.org/donors.xlsx`", md)
        self.assertNotIn("s3://web", md)

    def test_column_starting_partway_counts_its_blanks(self):
        r, _ = self.inspect(workbook([("Donors", False, DONORS)]))
        self.assertEqual(self.col(r["workbook"]["sheets"][0], "Note")["missing"], 4)

    def test_hidden_sheet_named(self):
        r, md = self.inspect(workbook([("Donors", False, DONORS),
                                       ("old", True, [["a", "b"], [1, 2]])]))
        self.assertTrue(r["workbook"]["sheets"][1]["hidden"])
        self.assertIn("## Sheet 'old' (hidden)", md)

    def test_few_donor_sheet_hides_levels(self):
        r, _ = self.inspect(workbook([("One", False, DONORS[:1] + DONORS[1:2] * 3)]))
        self.assertIsNone(self.col(r["workbook"]["sheets"][0], "Sex")["levels"])

    def test_zip_that_is_not_a_workbook(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr("readme.txt", "hi")
        r, _ = self.inspect(buf.getvalue())
        self.assertIn("could not read it", r["status"])

    def test_oversized_archive_refused(self):
        saved = xlsx.MAX_UNPACKED
        xlsx.MAX_UNPACKED = 100
        try:
            r, _ = self.inspect(workbook([("Donors", False, DONORS)]))
        finally:
            xlsx.MAX_UNPACKED = saved
        self.assertIn("unpacks to", r["warnings"][0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
