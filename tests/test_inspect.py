"""
test_inspect.py — unit tests for scripts/inspect_data.py (structure reports).

No network, no real data: every file is a tiny synthetic one written into
a temporary folder. Values that must never appear in a report are planted
on purpose (a donor ID, a barcode, measurement numbers) and the tests
assert they are absent from both the .md and the .json.

    python -m unittest discover -s tests -t tests -v

Why these tests: the inspect report is meant to be shared with the next
person. Its hazards are leaking identifiers or values, trusting a file
name over its bytes, opening an unsafe pickle, and missing that a file
changed after it was downloaded.
"""

import csv
import gzip
import json
import pathlib
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import inspect_data as ins  # noqa: E402

SC = "sea-ad-single-cell-profiling"
SECRET_DONOR = "H21.33.999"
SECRET_BARCODE = "ACGTACGTACGTACGT-1"
SECRET_NUMBER = "31415.9265"


def rows(n=50):
    out = []
    for i in range(n):
        out.append({
            "cell_label": f"cell{i:05d}",                 # unique per row
            # 6 donors: above the per-donor floor, so levels may show
            "Donor ID": SECRET_DONOR if i % 6 == 1 else f"H20.33.{i % 6:03d}",
            "anything": SECRET_BARCODE if i == 0 else "ACGTTTGGCCAAGGTT-1",
            "Braak": ["Braak II", "Braak IV"][i % 2],
            "Sex": ["Female", "Male"][i % 2],
            "Severely Affected Donor": ["N", "Y"][i % 2],
            "PMI": SECRET_NUMBER if i == 3 else str(4 + i % 5),
            "Supertype": f"type_{i % 40}",                # 40 levels > cap
            "empty_col": "",
            # one label per library: 25 levels is under the cap, so only
            # the name rule keeps them out (the real LEC file has 34)
            "load_name": f"SECRETLIB{i % 25:02d}",
            "rna_amplification": f"SECRETAMP{i % 25:02d}",
            "rna_amplification_pass_fail": "Pass",
        })
    return out


class InspectCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self.tmp.name)
        self.data, self.out = base / "data_cache", base / "results"
        self.path = self.data / SC / "LEC/RNAseq/meta.csv"
        self.path.parent.mkdir(parents=True)

    def tearDown(self):
        self.tmp.cleanup()

    def write_csv(self, path=None, compress=False):
        path = path or self.path
        data = rows()
        opener = (lambda p: gzip.open(p, "wt", newline="")) if compress else \
                 (lambda p: open(p, "w", newline=""))
        with opener(path) as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
        return path

    def inspect(self, path=None):
        return ins.inspect(path or self.path, data_dir=self.data, out_dir=self.out)

    def col(self, report, name):
        return next(c for c in report["columns"] if c["name"] == name)


class TestCsvStructure(InspectCase):

    def test_shape_types_and_missing(self):
        self.write_csv()
        r, _ = self.inspect()
        self.assertEqual(r["rows"], 50)
        self.assertEqual(len(r["columns"]), 12)
        self.assertEqual(self.col(r, "PMI")["type"], "float")
        self.assertTrue(self.col(r, "empty_col")["empty"])
        self.assertEqual(self.col(r, "empty_col")["type"], "empty")   # not pandas' float
        self.assertEqual(self.col(r, "empty_col")["missing"], 50)

    def test_few_level_categories_show_levels(self):
        self.write_csv()
        r, _ = self.inspect()
        self.assertEqual(self.col(r, "Braak")["levels"], ["Braak II", "Braak IV"])
        self.assertEqual(self.col(r, "Sex")["levels"], ["Female", "Male"])

    def test_more_than_cap_levels_are_counted_not_listed(self):
        self.write_csv()
        r, _ = self.inspect()
        c = self.col(r, "Supertype")
        self.assertEqual((c["distinct"], c["levels"]), (40, None))

    def test_yes_no_flag_shown_despite_donor_in_name(self):
        self.write_csv()
        r, _ = self.inspect()
        self.assertEqual(self.col(r, "Severely Affected Donor")["levels"], ["N", "Y"])

    def test_gzipped_csv(self):
        p = self.write_csv(self.path.with_name("meta.csv.gz"), compress=True)
        r, _ = self.inspect(p)
        self.assertEqual((r["detected_format"], r["rows"]), ("gzip", 50))


class TestNothingLeaks(InspectCase):

    def test_identifier_columns_are_withheld(self):
        self.write_csv()
        r, _ = self.inspect()
        self.assertIsNone(self.col(r, "Donor ID")["levels"])
        self.assertIsNotNone(self.col(r, "Donor ID")["levels_withheld"])
        self.assertIsNotNone(self.col(r, "cell_label")["levels_withheld"])

    def test_identifier_values_under_an_innocent_name_are_withheld(self):
        self.write_csv()
        r, _ = self.inspect()
        self.assertIn("donor IDs or barcodes", self.col(r, "anything")["levels_withheld"])

    def test_planted_values_appear_in_neither_output(self):
        self.write_csv()
        r, md = self.inspect()
        blob = md + json.dumps(r)
        for secret in (SECRET_DONOR, SECRET_BARCODE, SECRET_NUMBER, "cell00007"):
            self.assertNotIn(secret, blob, secret)

    def test_per_library_labels_withheld_even_under_the_cap(self):
        self.write_csv()
        r, md = self.inspect()
        for name in ("load_name", "rna_amplification"):
            self.assertEqual(self.col(r, name)["distinct"], 25)
            self.assertIsNone(self.col(r, name)["levels"], name)
            self.assertIsNotNone(self.col(r, name)["levels_withheld"], name)
        self.assertNotIn("SECRETLIB", md + json.dumps(r))
        self.assertNotIn("SECRETAMP", md + json.dumps(r))
        # its pass/fail companion is a flag and stays visible
        self.assertEqual(self.col(r, "rna_amplification_pass_fail")["levels"], ["Pass"])

    def test_numeric_columns_never_list_levels(self):
        self.write_csv()
        r, _ = self.inspect()
        c = self.col(r, "PMI")
        self.assertIsNone(c["levels"])
        self.assertIsNone(c["levels_withheld"])


class TestTextThatIsNotATable(InspectCase):
    """Regression: `inspect --downloaded` crashed on a UCSC trackDB.txt."""

    TRACKDB = ("hub HumanATAC-seq Hub\nshortLabel  ATAC-seq Hub\n\ngenome hg38\n\n"
               "track Lamp5Lhx6\ncontainer multiWig\n\tparent Lamp5Lhx6\n\tcolor 58,37,31\n")

    def test_track_config_is_plain_text_not_a_crash(self):
        p = self.path.with_name("trackDB.txt")
        p.write_text(self.TRACKDB)
        r, md = self.inspect(p)
        self.assertIn("plain text, not a table", r["status"])
        self.assertEqual(r["text_lines"], self.TRACKDB.count("\n"))
        self.assertNotIn("columns", r)
        self.assertNotIn("Lamp5Lhx6", md + json.dumps(r))

    def test_tsv_is_split_on_tabs(self):
        p = self.path.with_name("t.tsv")
        p.write_text("a\tb\tc\n1\t2,5\t3\n4\t5\t6\n")
        r, _ = self.inspect(p)
        self.assertEqual(([c["name"] for c in r["columns"]], r["delimiter"]),
                         (["a", "b", "c"], "tab"))

    def test_tab_separated_txt_is_a_table(self):
        p = self.path.with_name("t.txt")
        p.write_text("gene\tscore\nA\t1\nB\t2\n")
        r, _ = self.inspect(p)
        self.assertEqual(r["rows"], 2)

    def test_ragged_csv_is_reported_not_raised(self):
        p = self.path.with_name("bad.csv")
        p.write_text("a,b\n1,2\n3,4,5,6\n")
        r, md = self.inspect(p)
        self.assertEqual(r["status"], "could not be parsed as a table")
        self.assertTrue(r["warnings"])
        self.assertNotIn("columns", r)

    def test_one_failure_does_not_stop_the_rest(self):
        import contextlib
        import io
        import index
        good = self.write_csv()
        calls = []
        real = ins.inspect

        def flaky(p, *a, **k):
            calls.append(p)
            if str(p).endswith("boom.csv"):
                raise RuntimeError("simulated")
            return real(p, data_dir=self.data, out_dir=self.out)
        ins.inspect = flaky
        try:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                args = index.build_parser().parse_args(
                    ["inspect", str(self.path.with_name("boom.csv")), str(good)])
                code = args.func(args)
        finally:
            ins.inspect = real
        self.assertEqual(code, 1)
        self.assertEqual(len(calls), 2)
        self.assertIn("FAILED to inspect", buf.getvalue())
        self.assertIn("1 of 2 file(s) could not be inspected", buf.getvalue())

    def test_downloaded_with_nothing_downloaded_says_so(self):
        # A fresh clone has no data_cache/: the message must not tell the user
        # to pass --downloaded when they just did.
        import index
        real = ins.downloaded_files
        ins.downloaded_files = lambda *a, **k: []
        try:
            args = index.build_parser().parse_args(["inspect", "--downloaded"])
            with self.assertRaises(SystemExit) as cm:
                args.func(args)
        finally:
            ins.downloaded_files = real
        self.assertIn("Nothing to inspect", str(cm.exception.code))
        self.assertNotIn("Give a file path", str(cm.exception.code))


class TestFreeText(InspectCase):
    """Real case: a donor table's 'If other Consensus dx, describe' column."""

    def test_free_text_withheld_by_name_and_by_length(self):
        data = [{"Donor ID": f"H20.33.{i:03d}", "Sex": ["Female", "Male"][i % 2],
                 "If other Consensus dx, describe": "SECRETDX early onset" if i < 3 else "",
                 "specify other race": "SECRETRACE" if i == 0 else "",
                 "Summary": ("SECRETLONG " + "x" * 50) if i == 1 else "short"}
                for i in range(8)]
        with open(self.path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
        r, md = self.inspect()
        for name in ("If other Consensus dx, describe", "specify other race", "Summary"):
            self.assertIsNone(self.col(r, name)["levels"], name)
            self.assertIn("free text", self.col(r, name)["levels_withheld"], name)
        self.assertEqual(self.col(r, "Sex")["levels"], ["Female", "Male"])
        self.assertNotIn("SECRET", md + json.dumps(r))


class TestFewDonors(InspectCase):
    """A per-donor file must not list that donor's own values."""

    def write_rows(self, path, data):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)

    def one_donor(self):
        return [{"Donor ID": SECRET_DONOR, "Sex": "SECRETSEX", "APOE": "SECRET-2/3",
                 "Subclass": ["Astro", "Oligo", "Vip"][i % 3], "blank": ""}
                for i in range(30)]

    def test_single_donor_hides_donor_level_columns_only(self):
        self.write_rows(self.path, self.one_donor())
        r, md = self.inspect()
        self.assertEqual(r["donor_floor"]["donors"], 1)
        for name in ("Sex", "APOE"):
            self.assertIsNone(self.col(r, name)["levels"], name)
            self.assertIn("only 1 donor", self.col(r, name)["levels_withheld"])
        # varies between cells, so it describes cells, not the person
        self.assertEqual(self.col(r, "Subclass")["levels"], ["Astro", "Oligo", "Vip"])
        self.assertNotIn("SECRET", md + json.dumps(r))
        self.assertIn("Donor-level levels withheld", md)
        # an empty column has nothing to withhold (real ATAC report said otherwise)
        self.assertEqual((self.col(r, "blank")["type"], self.col(r, "blank")["levels_withheld"]),
                         ("empty", None))

    def test_donor_id_in_the_file_name_counts_as_one_donor(self):
        p = self.path.with_name("H20.33.036_meta.csv")
        self.write_rows(p, [{"Sex": "SECRETSEX", "x": str(i)} for i in range(10)])
        r, md = self.inspect(p)
        self.assertEqual(r["donor_floor"]["source"], "the file name names one donor")
        self.assertNotIn("SECRETSEX", md + json.dumps(r))

    def test_many_donors_no_floor(self):
        self.write_csv()
        r, _ = self.inspect()
        self.assertIsNone(r["donor_floor"])


class TestFormatDetection(InspectCase):

    def test_bytes_win_over_the_extension(self):
        fake = self.path.with_name("actually_gzip.csv")
        with gzip.open(fake, "wb") as fh:
            fh.write(b"not,a,real,csv\n")
        r, _ = self.inspect(fake)
        self.assertEqual(r["detected_format"], "gzip")
        self.assertTrue(any("bytes say gzip" in w for w in r["warnings"]))

    def test_pickle_is_never_opened(self):
        p = self.path.with_name("menr.pkl")
        # Protocol-4 header followed by garbage: loading it would fail, so
        # a crash would prove someone tried. The report must not.
        p.write_bytes(b"\x80\x04" + b"\x00" * 64)
        r, md = self.inspect(p)
        self.assertEqual(r["detected_format"], "pickle")
        self.assertIn("not opened", r["status"])
        self.assertNotIn("columns", r)

    def test_parquet_detected_and_explained(self):
        p = self.path.with_name("t.parquet")
        p.write_bytes(b"PAR1" + b"\x00" * 32 + b"PAR1")
        r, _ = self.inspect(p)
        self.assertEqual(r["detected_format"], "parquet")
        self.assertIn("pyarrow", r["status"])

    def test_damaged_hdf5_is_reported_not_crashed(self):
        p = self.path.with_name("x.h5ad")
        p.write_bytes(b"\x89HDF\r\n\x1a\n" + b"\x00" * 32)
        r, _ = self.inspect(p)
        self.assertEqual(r["detected_format"], "hdf5")
        self.assertIn("could not be opened", r["status"])
        self.assertTrue(any("damaged" in w for w in r["warnings"]))


HG38 = [("chr1", 248_956_422), ("chr2", 242_193_529), ("chrX", 156_040_895),
        ("chrY", 57_227_415), ("chrM", 16_569), ("chr1_KI270706v1_random", 175_055)]


def bigwig_bytes(chroms, endian="<", two_level=False):
    """A minimal bigWig built from the format spec: 64-byte header, then the
    chromosome B+ tree, then stand-in 'signal' bytes that must never be read."""
    key = max(len(c) for c, _ in chroms)

    def leaf(items, first):
        b = struct.pack(endian + "BBH", 1, 0, len(items))
        for i, (c, size) in enumerate(items, first):
            b += c.encode().ljust(key, b"\0") + struct.pack(endian + "II", i, size)
        return b

    tree_at = 64
    if not two_level:
        nodes = leaf(chroms, 0)
    else:                                   # root -> two leaves
        half = len(chroms) // 2
        root_len = 4 + 2 * (key + 8)
        l1 = leaf(chroms[:half], 0)
        at1 = tree_at + 32 + root_len
        at2 = at1 + len(l1)
        root = struct.pack(endian + "BBH", 0, 0, 2)
        root += chroms[0][0].encode().ljust(key, b"\0") + struct.pack(endian + "Q", at1)
        root += chroms[half][0].encode().ljust(key, b"\0") + struct.pack(endian + "Q", at2)
        nodes = root + l1 + leaf(chroms[half:], half)
    tree = struct.pack(endian + "IIIIQQ", 0x78CA8C91, 256, key, 8, len(chroms), 0) + nodes
    data_at = tree_at + len(tree)
    header = struct.pack(endian + "IHHQQQHHQQIQ", 0x888FFC26, 4, 10, tree_at, data_at,
                         0, 0, 0, 0, 0, 32768, 0)
    return header + tree + b"SIGNAL31415.9265" * 8


class TestBigWig(InspectCase):

    def write(self, data, name="ADNC0Astro.bw"):
        p = self.path.with_name(name)
        p.write_bytes(data)
        return p

    def test_header_genome_build_and_chromosomes(self):
        r, md = self.inspect(self.write(bigwig_bytes(HG38)))
        self.assertEqual(r["detected_format"], "bigwig")
        t = r["track"]
        self.assertEqual(t["genome_build"], "GRCh38 / hg38")
        self.assertEqual(t["chromosomes"]["primary"], ["chr1", "chr2", "chrX", "chrY", "chrM"])
        self.assertEqual(t["chromosomes"]["other"], 1)
        self.assertEqual((t["zoom_levels"], t["compressed"]), (10, True))
        self.assertIn("header only", r["status"])

    def test_two_level_tree_and_big_endian(self):
        r, _ = self.inspect(self.write(bigwig_bytes(HG38, endian=">", two_level=True)))
        self.assertEqual(r["track"]["chromosomes"]["n"], len(HG38))
        self.assertEqual(r["track"]["genome_build"], "GRCh38 / hg38")

    def test_hg19_is_told_apart(self):
        r, _ = self.inspect(self.write(bigwig_bytes([("chr1", 249_250_621)])))
        self.assertEqual(r["track"]["genome_build"], "GRCh37 / hg19")

    def test_signal_is_never_read(self):
        r, md = self.inspect(self.write(bigwig_bytes(HG38)))
        self.assertNotIn("31415", md + json.dumps(r))

    def test_damaged_header_is_reported(self):
        r, _ = self.inspect(self.write(bigwig_bytes(HG38)[:70]))
        self.assertIn("could not be read", r["status"])
        self.assertTrue(r["warnings"])


class TestLedgerCheck(InspectCase):

    def ledger(self, sha, path=None):
        led = self.data / "_ledger.csv"
        with open(led, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["downloaded_at", "bucket", "key", "size", "etag",
                        "s3_last_modified", "sha256", "verified", "path"])
            w.writerow(["2026-11-22", SC, "LEC/RNAseq/meta.csv", 1, "e", "", sha,
                        "md5 match", str(self.path) if path is None else path])

    def test_match(self):
        self.write_csv()
        self.ledger(ins.sha256(self.path))
        r, _ = self.inspect()
        self.assertTrue(r["ledger"].startswith("matches ledger"))
        self.assertEqual(r["warnings"], [])

    def test_mismatch_is_a_warning(self):
        self.write_csv()
        self.ledger("0" * 64)
        r, md = self.inspect()
        self.assertIn("DIFFERS", r["ledger"])
        self.assertIn("WARNING", md)

    def test_relative_ledger_row_matches(self):
        self.write_csv()
        self.ledger(ins.sha256(self.path), path=f"{SC}/LEC/RNAseq/meta.csv")
        r, _ = self.inspect()
        self.assertTrue(r["ledger"].startswith("matches ledger"))

    def test_downloaded_files_reads_relative_rows(self):
        self.write_csv()
        self.ledger("0" * 64, path=f"{SC}/LEC/RNAseq/meta.csv")
        self.assertEqual(ins.downloaded_files(self.data), [str(self.path)])

    def test_not_in_ledger(self):
        self.write_csv()
        r, _ = self.inspect()
        self.assertIn("not in the download ledger", r["ledger"])


class TestOutputs(InspectCase):

    def test_saved_report_has_no_absolute_path(self):
        # Reports are committed and shared: no user name or home folder in them.
        self.write_csv()
        r, md = self.inspect()
        saved = (self.out / SC / "LEC/RNAseq/meta.csv.json").read_text()
        self.assertEqual(r["file"], "meta.csv")            # outside the project: name only
        for text in (saved, md):
            self.assertNotIn(str(self.data), text)
            self.assertNotIn(str(self.data.resolve()), text)

    def test_shown_path_is_relative_inside_the_project(self):
        f = ins.ROOT / "data_cache" / SC / "x.csv"
        self.assertEqual(ins.shown_path(f), f"data_cache/{SC}/x.csv")

    def test_json_and_md_written_under_bucket_and_key(self):
        self.write_csv()
        r, _ = self.inspect()
        base = self.out / SC / "LEC/RNAseq"
        self.assertTrue((base / "meta.csv.json").exists())
        self.assertTrue((base / "meta.csv.md").exists())
        self.assertEqual(r["key"], "LEC/RNAseq/meta.csv")

    def test_source_is_recovered_from_the_path(self):
        self.assertEqual(ins.source_of(self.path, self.data), (SC, "LEC/RNAseq/meta.csv"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
