"""
test_images.py — unit tests for scripts/images.py (image headers) and scripts/peek.py.

No network, no real data. TIFFs are built byte by byte from the TIFF 6.0 /
BigTIFF specification (so every tag and pointer is known), plus one written
by Pillow as a check against a real writer. Text that must never reach a
report is planted: Aperio file names and dates, a free-text description,
and stand-in pixel bytes.

    python -m unittest discover -s tests -t tests -v

Why these tests: an image header report must get dimensions, pyramid
levels and pixel type right, must not leak what scanners write into
description text, must survive damaged files, and peek must read only a
few kilobytes of a huge file, stop at its byte limit, and refuse
non-image files.
"""

import csv
import gzip
import io
import json
import pathlib
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import images  # noqa: E402
import index  # noqa: E402
import inspect_data as ins  # noqa: E402
import peek  # noqa: E402
from test_sources import FakeS3  # noqa: E402

NP = "sea-ad-quantitative-neuropathology"
PIXELS = b"SECRETPIXEL" * 64
APERIO = ("Aperio Image Library v12.0.15\r\n4000x3000 [0,0 4000x3000] (256x256) JPEG/RGB Q=70"
          "|AppMag = 20|StripeWidth = 2040|ScanScope ID = SECRETSCANNER|Filename = SECRETFILE"
          "|Date = 01/02/03|Time = 10:11:12|MPP = 0.5040")
OME = ('<?xml version="1.0"?><OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
       '<Image ID="Image:0"><Pixels DimensionOrder="XYCZT" Type="uint16" SizeX="2048" '
       'SizeY="1024" SizeZ="1" SizeC="3" SizeT="1" PhysicalSizeX="0.2125" '
       'PhysicalSizeXUnit="µm" PhysicalSizeY="0.2125" PhysicalSizeYUnit="µm">'
       '<Channel ID="C0" Name="DAPI"/><Channel ID="C1" Name="PolyT"/>'
       '<Channel ID="C2" Name="Aux1"/></Pixels></Image></OME>')


# ---- a TIFF writer that follows the spec ------------------------------------

def page(w, h, bits=8, samples=1, compression=1, tile=None, desc=None, subfile=0):
    t = {256: (4, [w]), 257: (4, [h]), 258: (3, [bits] * samples), 259: (3, [compression]),
         262: (3, [2 if samples == 3 else 1]), 277: (3, [samples]), 254: (4, [subfile])}
    if tile:
        t[322], t[323] = (3, [tile]), (3, [tile])
    if desc is not None:
        t[270] = (2, desc)
    return t


def build_tiff(pages, big=False, endian="<", pixels_first=True, subifds=(), loop=False,
               lead=1):
    """Bytes of a TIFF: header, (pixels x lead), sub-IFDs, then the page IFDs."""
    E = endian
    hdr = 16 if big else 8
    buf = bytearray(hdr)
    if pixels_first:
        buf += PIXELS * lead

    def put_ifd(tags):
        """Append one IFD (+ its overflow values); return (offset, next-pointer position)."""
        at = len(buf)
        n = len(tags)
        cnt_fmt, ent, ptr, inline = ("Q", 20, "Q", 8) if big else ("H", 12, "I", 4)
        size = struct.calcsize(E + cnt_fmt) + n * ent + struct.calcsize(E + ptr)
        overflow = bytearray()
        entries = b""
        for tag in sorted(tags):
            typ, val = tags[tag]
            if typ == 2:
                data = val.encode() + b"\x00"
                count = len(data)
            else:
                code = images.TYPES[typ][0]
                data = struct.pack(E + code * len(val), *val)
                count = len(val)
            if len(data) <= inline:
                field = data.ljust(inline, b"\x00")
            else:
                field = struct.pack(E + ptr, at + size + len(overflow))
                overflow += data
                if len(overflow) % 2:
                    overflow += b"\x00"
            entries += struct.pack(E + ("HHQ" if big else "HHI"), tag, typ, count) + field
        buf.extend(struct.pack(E + cnt_fmt, n) + entries)
        next_at = len(buf)
        buf.extend(struct.pack(E + ptr, 0))
        buf.extend(overflow)
        return at, next_at

    sub_offsets = [put_ifd(t)[0] for t in subifds]
    first_tags = dict(pages[0])
    if sub_offsets:
        first_tags[330] = (16 if big else 13, sub_offsets)
    offsets = []
    prev_next = None
    for k, t in enumerate([first_tags] + list(pages[1:])):
        at, nxt = put_ifd(t)
        offsets.append(at)
        if prev_next is not None:
            struct.pack_into(E + ("Q" if big else "I"), buf, prev_next, at)
        prev_next = nxt
    if loop:                                  # last page points back at the first
        struct.pack_into(E + ("Q" if big else "I"), buf, prev_next, offsets[0])
    if not pixels_first:
        buf += PIXELS * 4096                  # ~2.9 MB of 'pixels' after the header
    bo = b"II" if E == "<" else b"MM"
    if big:
        buf[0:16] = bo + struct.pack(E + "HHHQ", 43, 8, 0, offsets[0])
    else:
        buf[0:8] = bo + struct.pack(E + "HI", 42, offsets[0])
    return bytes(buf)


def svs_bytes(**kw):
    return build_tiff([
        page(4000, 3000, samples=3, compression=7, tile=256, desc=APERIO),
        page(400, 300, samples=3, compression=7, desc="Aperio Image Library v12\r\n400x300 thumb"),
        page(1000, 750, samples=3, compression=7, tile=256, desc="Aperio Image Library v12\r\n1000x750"),
        page(400, 300, samples=3, compression=5, desc="Aperio Image Library v12\r\nlabel 400x300"),
        page(800, 300, samples=3, compression=7, desc="Aperio Image Library v12\r\nmacro 800x300"),
    ], **kw)


class ImageCase(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = pathlib.Path(self.tmp.name)
        self.data = self.base / "data_cache"

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, key, data):
        p = self.data / NP / key
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data) if isinstance(data, bytes) else p.write_text(data)
        return p

    def inspect(self, p):
        return ins.inspect(p, data_dir=self.data, out_dir=self.base / "out")


class TestAperioSlide(ImageCase):

    def setUp(self):
        super().setUp()
        self.r, self.md = self.inspect(self.write("HIP/slide.svs", svs_bytes()))
        self.img = self.r["image"]

    def test_detected_and_counted(self):
        self.assertEqual(self.r["detected_format"], "tiff")
        self.assertEqual((self.img["container"], self.img["n_pages"]), ("TIFF", 5))
        self.assertIn("no pixels", self.r["status"])

    def test_dimensions_and_pixel_type(self):
        p0 = self.img["pages"][0]
        self.assertEqual((p0["width"], p0["height"], p0["channels"], p0["bits"]),
                         (4000, 3000, 3, 8))
        self.assertEqual((p0["compression"], p0["layout"]), ("JPEG", "tiles 256×256"))

    def test_only_magnification_and_pixel_size_taken_from_notes(self):
        self.assertEqual(self.img["aperio"], {"magnification": 20.0, "microns_per_pixel": 0.504})
        blob = self.md + json.dumps(self.r)
        for secret in ("SECRETSCANNER", "SECRETFILE", "01/02/03", "10:11:12", "SECRETPIXEL"):
            self.assertNotIn(secret, blob, secret)

    def test_pages_named_by_role(self):
        roles = [p["role"] for p in self.img["pages"]]
        self.assertEqual(roles, ["full resolution", "thumbnail", "1/4 resolution",
                                 "slide label photo", "whole-slide photo (macro)"])


class TestOmeBigTiff(ImageCase):

    def test_bigtiff_big_endian_ome_with_pyramid(self):
        data = build_tiff([page(2048, 1024, bits=16, tile=512, desc=OME)], big=True,
                          endian=">", subifds=[page(1024, 512, bits=16, tile=512, subfile=1),
                                               page(512, 256, bits=16, tile=512, subfile=1)])
        r, md = self.inspect(self.write("ST/DAPI.ome.tif", data))
        img = r["image"]
        self.assertEqual((img["container"], img["byte_order"]), ("BigTIFF", "big-endian"))
        self.assertEqual([p["width"] for p in img["sub_resolutions"]], [1024, 512])
        s = img["ome"]["series"][0]
        self.assertEqual(s["size"], {"X": 2048, "Y": 1024, "Z": 1, "C": 3, "T": 1})
        self.assertEqual(s["channels"], ["DAPI", "PolyT", "Aux1"])
        self.assertEqual(s["pixel_size"]["X"], "0.2125 µm")
        self.assertIn("16", md)

    def test_one_page_per_channel_named_by_channel(self):
        data = build_tiff([page(2048, 1024, bits=16, tile=512, desc=OME),
                           page(2048, 1024, bits=16, tile=512),
                           page(2048, 1024, bits=16, tile=512)], big=True)
        r, md = self.inspect(self.write("ST/IF.ome.tif", data))
        self.assertEqual([p.get("role") for p in r["image"]["pages"]],
                         ["channel: DAPI", "channel: PolyT", "channel: Aux1"])

    def test_pixel_size_rounded(self):
        data = build_tiff([page(8, 8, desc=OME.replace('PhysicalSizeX="0.2125"',
                                                       'PhysicalSizeX="0.3250006296210897"'))])
        r, _ = self.inspect(self.write("ST/r.ome.tif", data))
        self.assertEqual(r["image"]["ome"]["series"][0]["pixel_size"]["X"], "0.325 µm")


class TestOtherWriters(ImageCase):

    def test_pillow_tiff_free_text_not_shown(self):
        from PIL import Image
        buf = io.BytesIO()
        frames = [Image.new("L", (64, 32)), Image.new("L", (32, 16))]
        frames[0].save(buf, format="TIFF", save_all=True, append_images=frames[1:],
                       description="SECRETDESC lab notes")
        r, md = self.inspect(self.write("x/mask.tif", buf.getvalue()))
        img = r["image"]
        self.assertEqual(img["n_pages"], 2)
        self.assertEqual((img["pages"][0]["width"], img["pages"][0]["height"]), (64, 32))
        self.assertEqual(img["description"], "other free text")
        self.assertNotIn("SECRETDESC", md + json.dumps(r))


class TestCompanionImage(ImageCase):

    def test_same_size_sub_image_is_not_called_reduced(self):
        data = build_tiff([page(30, 70, samples=3, tile=800)],
                          subifds=[page(30, 70, tile=800)])
        _, md = self.inspect(self.write("x/m.tif", data))
        self.assertIn("same-size companion image", md)
        self.assertNotIn("reduced resolution", md.split("### How")[0])


class TestDamaged(ImageCase):

    def test_loop_in_directory_chain(self):
        r, _ = self.inspect(self.write("x/loop.tif", build_tiff(
            [page(10, 10), page(5, 5)], loop=True)))
        self.assertIn("could not be read", r["status"])
        self.assertTrue(any("loops" in w for w in r["warnings"]))

    def test_truncated_file(self):
        r, _ = self.inspect(self.write("x/cut.tif", svs_bytes()[:900]))
        self.assertIn("could not be read", r["status"])


class TestDeepZoom(ImageCase):

    DZI = ('<?xml version="1.0" encoding="UTF-8"?><Image TileSize="256" Overlap="1" '
           'Format="jpeg" xmlns="http://schemas.microsoft.com/deepzoom/2008">'
           '<Size Width="1024" Height="1024"/></Image>')

    def test_geometry_and_expected_tiles(self):
        d = images.dzi(self.DZI)
        # levels 1..256 px: 9 single tiles; 512 px: 2×2; 1024 px: 4×4  -> 29
        self.assertEqual((d["levels"], d["expected_tiles"]), (11, 29))

    def test_tile_folder_compared_with_the_catalog(self):
        files = self.base / "catfiles"
        files.mkdir()
        with gzip.open(files / f"{NP}.csv.gz", "wt", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["bucket", "unit", "kind", "n_objects"])
            w.writerow([NP, "MTG/s1/s1_files/", "collapsed", 20])
        saved = ins.CATALOG_FILES
        ins.CATALOG_FILES = files
        try:
            r, md = self.inspect(self.write("MTG/s1/s1.dzi", self.DZI))
        finally:
            ins.CATALOG_FILES = saved
        self.assertEqual(r["deep_zoom"]["catalog_tiles"], 20)
        self.assertIn("FEWER than expected: 9 missing", md)

    def test_one_extra_object_is_the_libvips_properties_file(self):
        files = self.base / "catfiles"
        files.mkdir()
        with gzip.open(files / f"{NP}.csv.gz", "wt", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["bucket", "unit", "kind", "n_objects"])
            w.writerow([NP, "MTG/s1/s1_files/", "collapsed", 30])
        saved = ins.CATALOG_FILES
        ins.CATALOG_FILES = files
        try:
            _, md = self.inspect(self.write("MTG/s1/s1.dzi", self.DZI))
        finally:
            ins.CATALOG_FILES = saved
        self.assertIn("matches, plus one non-tile object", md)

    def test_other_xml_is_not_a_dzi(self):
        r, _ = self.inspect(self.write("x/odd.dzi", "<Other/>"))
        self.assertIn("not a deep-zoom descriptor", r["status"])


# ---- peek ---------------------------------------------------------------------

class FakeObjectStore(FakeS3):
    """FakeS3 listing plus head/get with byte ranges and If-Match, from {key: bytes}."""

    def __init__(self, objects):
        super().__init__({k: len(v) for k, v in objects.items()})
        self.objects = objects
        self.ranges = []

    def etag(self, key):
        return f"etag-{len(key)}"            # what FakeS3's listing reports

    def head_object(self, Bucket, Key):
        return {"ContentLength": len(self.objects[Key]), "ETag": f'"{self.etag(Key)}"'}

    def get_object(self, Bucket, Key, Range=None, IfMatch=None):
        if IfMatch is not None and IfMatch.strip('"') != self.etag(Key):
            raise RuntimeError("PreconditionFailed")
        a, b = map(int, Range.split("=")[1].split("-"))
        self.ranges.append((a, b))
        return {"Body": io.BytesIO(self.objects[Key][a:b + 1])}


class TestPeek(unittest.TestCase):

    SVS = "MTG/H00.00.000/slide/slide.svs"
    CSV = "MTG/H00.00.000/slide/table.csv"
    FAR = "MTG/H00.00.000/slide/far.svs"

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self.tmp.name)
        self.cat, self.data, self.out = base / "catalog", base / "data_cache", base / "out"
        self.big = svs_bytes(pixels_first=False)     # header first, ~2.9 MB of pixels after
        self.store = FakeObjectStore({self.SVS: self.big, self.CSV: b"a,b\n1,2\n"})
        index.refresh(NP, client=self.store, catalog_dir=self.cat, out=lambda *a: None)
        self.lines = []

    def tearDown(self):
        self.tmp.cleanup()

    def run_peek(self, match="slide.svs", **kw):
        kw.setdefault("yes", True)
        return peek.run(match, client=self.store, catalog_dir=self.cat, data_dir=self.data,
                        out_dir=self.out, out=self.lines.append, **kw)

    def test_reads_a_few_kilobytes_of_a_big_file(self):
        code = self.run_peek()
        self.assertEqual(code, 0)
        rep = json.loads((self.out / NP / (self.SVS + ".json")).read_text())
        self.assertEqual(rep["image"]["n_pages"], 5)
        read = sum(b - a + 1 for a, b in self.store.ranges)
        self.assertLess(read, 100_000)
        self.assertGreater(len(self.big), 2_000_000)
        self.assertFalse(any(self.data.glob(f"{NP}/**/*.svs")))     # nothing saved
        log = list(csv.DictReader((self.data / "_peek_log.csv").read_text().splitlines()))
        self.assertEqual(int(log[0]["bytes_read"]), read)

    def test_byte_limit_stops_the_read(self):
        # header directories placed ~140 KB into the file: beyond a 64 KB limit
        self.store.objects[self.FAR] = svs_bytes(lead=200)
        self.store.keys[self.FAR] = len(self.store.objects[self.FAR])
        index.refresh(NP, client=self.store, catalog_dir=self.cat, out=lambda *a: None)
        code = self.run_peek(match="far.svs", max_kb=64)
        self.assertEqual(code, 7)
        self.assertIn("64 KB limit", "\n".join(self.lines))
        self.assertLessEqual(sum(b - a + 1 for a, b in self.store.ranges), 64 * 1024)

    def test_dry_run_reads_nothing(self):
        self.assertEqual(self.run_peek(dry_run=True), 0)
        self.assertEqual(self.store.ranges, [])

    def test_needs_a_yes(self):
        code = self.run_peek(yes=False, isatty=False)
        self.assertEqual(code, 6)
        self.assertEqual(self.store.ranges, [])

    def test_non_image_refused(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_peek(match="table.csv")
        self.assertIn("image headers only", str(cm.exception))

    def test_changed_object_is_skipped(self):
        self.store.objects[self.SVS] = self.big + b"x"
        code = self.run_peek()
        self.assertEqual(code, 7)
        self.assertIn("changed on S3", "\n".join(self.lines))
        self.assertEqual(self.store.ranges, [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
