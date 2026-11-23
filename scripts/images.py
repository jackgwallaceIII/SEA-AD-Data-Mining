"""
images.py — read the HEADER of an image file, never its pixels.

Not a step. Used by inspect_data.py (files on disk) and peek.py (files on
S3, read in small byte ranges). Covers TIFF and BigTIFF, which includes
Aperio .svs whole-slide images and OME-TIFF microscopy images, and the
small XML descriptors (.dzi) of deep-zoom tile pyramids.

Decisions made here, and why:

1. PURE PYTHON, FROM THE FORMAT SPECIFICATION.
   A TIFF header is a chain of small directories ("IFDs") listing tags:
   width, height, bits per sample, compression, tiling. Reading them needs
   only Python's struct module, so no imaging library is installed and
   nothing ever decodes a pixel. tifffile or SimpleITK would also read
   headers, but both are built to read pixels, and that is exactly what
   this tool must not do.

2. ONE READER INTERFACE, TWO SOURCES.
   The parser asks for bytes by position: read(offset, length). A file on
   disk answers from the file; a file on S3 answers with HTTP range
   requests (peek.py). So the same tested code reads a 20 KB mask on disk
   and the header of a 7 GB image that is never downloaded.

3. ONLY THE TAGS THAT DESCRIBE STRUCTURE ARE READ.
   Pixel locations (strip/tile offset tables, which can be megabytes) are
   never fetched. Tags that can carry personal or lab-internal text
   (date/time, artist, host computer) are ignored.

4. FREE-TEXT DESCRIPTIONS ARE MINED, NOT SHOWN.
   TIFF has an ImageDescription tag. In an Aperio .svs it holds the
   scanner's notes, including file names, dates and scanner IDs; only
   two numbers are kept from it: magnification (AppMag) and microns per
   pixel (MPP). In an OME-TIFF it is OME-XML, from which the image
   dimensions, physical pixel size and channel names are kept (channel
   names, e.g. DAPI, are the image's 'column names'). Any other
   description is reported as present but not shown.

5. SLIDE PAGES ARE NAMED BY ROLE.
   A .svs holds the full-resolution scan, a thumbnail, lower-resolution
   copies, and photographs of the slide label and of the whole glass
   ("macro"). The label photo can show handwriting; the report names
   the page 'label' but, as for every page, reads no pixels.

6. A DEEP-ZOOM DESCRIPTOR PREDICTS ITS TILE COUNT.
   A .dzi gives the full image size and tile size, which fixes how many
   tiles its _files/ folder should hold. The report compares that with
   the number the catalog listed: a shortfall means an incomplete pyramid.

7. LIMITS EVERYWHERE.
   At most 200 pages, loop detection on the directory chain, and at most
   2 MB of description text, so a damaged or hostile file cannot make the
   reader run away.
"""

import math
import re
import struct

MAX_PAGES = 200
MAX_DESCRIPTION = 2 * 1024 * 1024
MAX_CHANNELS = 50

TIFF_MAGIC = (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")

# type id -> (struct code, byte size)
TYPES = {1: ("B", 1), 2: ("s", 1), 3: ("H", 2), 4: ("I", 4), 5: ("II", 8),
         6: ("b", 1), 7: ("B", 1), 8: ("h", 2), 9: ("i", 4), 10: ("ii", 8),
         11: ("f", 4), 12: ("d", 8), 13: ("I", 4), 16: ("Q", 8), 17: ("q", 8),
         18: ("Q", 8)}
WANTED = {254: "subfile", 256: "width", 257: "height", 258: "bits", 259: "compression",
          262: "photometric", 270: "description", 277: "samples", 322: "tile_width",
          323: "tile_height", 330: "subifds", 339: "sample_format"}
COMPRESSION = {1: "none", 5: "LZW", 7: "JPEG", 8: "Deflate", 32946: "Deflate",
               32773: "PackBits", 33003: "JPEG 2000 (Aperio)", 33005: "JPEG 2000 (Aperio)",
               34712: "JPEG 2000", 50000: "Zstandard", 50001: "WebP", 6: "old JPEG"}
SAMPLE_FORMAT = {1: "unsigned int", 2: "signed int", 3: "float"}


class HeaderError(ValueError):
    pass


class LocalReader:
    """read(offset, n) from a file on disk."""

    def __init__(self, path):
        self.path = path
        self.fh = open(path, "rb")
        self.fh.seek(0, 2)
        self.size = self.fh.tell()
        self.bytes_read = 0

    def read(self, offset, n):
        if offset < 0 or offset >= self.size:
            raise HeaderError(f"offset {offset} outside the file")
        self.fh.seek(offset)
        data = self.fh.read(n)
        self.bytes_read += len(data)
        return data

    def close(self):
        self.fh.close()


# ---- TIFF -----------------------------------------------------------------

def _ifd(reader, offset, endian, big):
    """One directory: {tag name: value} for WANTED tags, plus next offset."""
    count_fmt, entry_size, ptr_fmt, inline = ("Q", 20, "Q", 8) if big else ("H", 12, "I", 4)
    head = reader.read(offset, 8 if big else 2)
    (n,) = struct.unpack(endian + count_fmt, head)
    if n > 4096:
        raise HeaderError(f"directory claims {n} entries: damaged?")
    body = reader.read(offset + len(head), n * entry_size + (8 if big else 4))
    if len(body) < n * entry_size + (8 if big else 4):
        raise HeaderError("directory runs past the end of the file")
    tags = {}
    for k in range(n):
        e = body[k * entry_size:(k + 1) * entry_size]
        if big:
            tag, typ, cnt = struct.unpack(endian + "HHQ", e[:12])
            raw = e[12:20]
        else:
            tag, typ, cnt = struct.unpack(endian + "HHI", e[:8])
            raw = e[8:12]
        if tag not in WANTED or typ not in TYPES:
            continue
        code, size = TYPES[typ]
        nbytes = cnt * size
        if tag == 270 and nbytes > MAX_DESCRIPTION:
            tags["description"] = None            # too long: noted, not read
            continue
        if nbytes > 1 << 20:                      # no structural tag is this big
            continue
        if nbytes <= inline:
            data = raw[:nbytes]
        else:
            (ptr,) = struct.unpack(endian + ptr_fmt, raw)
            data = reader.read(ptr, nbytes)
        if typ == 2:
            value = data.split(b"\x00", 1)[0].decode("utf-8", "replace")
        elif code in ("II", "ii"):                # rationals: numerator/denominator
            nums = struct.unpack(endian + code[0] * 2 * cnt, data)
            value = [nums[i] / nums[i + 1] if nums[i + 1] else 0
                     for i in range(0, len(nums), 2)]
        else:
            value = list(struct.unpack(endian + code * cnt, data))
        if isinstance(value, list) and len(value) == 1 and tag != 330:
            value = value[0]
        tags[WANTED[tag]] = value
    (nxt,) = struct.unpack(endian + ptr_fmt, body[n * entry_size:n * entry_size + inline])
    return tags, nxt


def _page_summary(t):
    bits = t.get("bits")
    if isinstance(bits, list):
        bits = bits[0] if len(set(bits)) == 1 else "/".join(map(str, bits))
    tiled = "tile_width" in t
    sub = t.get("subfile") or 0
    return {"width": t.get("width"), "height": t.get("height"),
            "channels": t.get("samples", 1), "bits": bits,
            "sample_format": SAMPLE_FORMAT.get(t.get("sample_format", 1), "other"),
            "compression": COMPRESSION.get(t.get("compression", 1),
                                           f"code {t.get('compression')}"),
            "layout": f"tiles {t['tile_width']}×{t.get('tile_height')}" if tiled else "strips",
            "reduced_resolution": bool(sub & 1), "mask": bool(sub & 4)}


def _aperio(text):
    """AppMag and MPP only; everything else in the text is withheld."""
    out = {}
    for key, name in (("AppMag", "magnification"), ("MPP", "microns_per_pixel")):
        m = re.search(rf"\|\s*{key}\s*=\s*([0-9.]+)", text)
        if m:
            out[name] = float(m.group(1))
    return out


def _ome(text):
    """Dimensions, physical pixel size and channel names from OME-XML."""
    import xml.etree.ElementTree as ET
    try:
        root = ET.fromstring(text)
    except ET.ParseError:
        return {"error": "OME-XML present but not parseable"}

    def local(tag):
        return tag.rsplit("}", 1)[-1]
    images = [e for e in root.iter() if local(e.tag) == "Image"]
    out = {"images": len(images), "series": []}
    for img in images[:20]:
        px = next((e for e in img.iter() if local(e.tag) == "Pixels"), None)
        if px is None:
            continue
        a = px.attrib
        chans = [c.attrib.get("Name", "") for c in px if local(c.tag) == "Channel"]
        out["series"].append({
            "size": {k: int(a[f"Size{k}"]) for k in "XYZCT" if f"Size{k}" in a},
            "type": a.get("Type"), "dimension_order": a.get("DimensionOrder"),
            "pixel_size": {k: f"{float(a[f'PhysicalSize{k}']):.4g} "
                              f"{a.get(f'PhysicalSize{k}Unit', 'µm')}"
                           for k in "XY" if f"PhysicalSize{k}" in a},
            "channels": chans[:MAX_CHANNELS], "n_channels": len(chans)})
    return out


def _description_kind(text):
    if text is None:
        return "too long to read (over 2 MB)"
    if text.startswith("Aperio"):
        return "Aperio scanner notes"
    if "<OME" in text[:2000]:
        return "OME-XML"
    if text.startswith("ImageJ="):
        return "ImageJ"
    return "other free text" if text.strip() else "empty"


def _svs_role(index, page, desc, first):
    """Name a page of an Aperio slide by what it is."""
    low = (desc or "").lower()
    if "label" in low.split("\n", 1)[-1][:80]:
        return "slide label photo"
    if "macro" in low.split("\n", 1)[-1][:80]:
        return "whole-slide photo (macro)"
    if index == 0:
        return "full resolution"
    if index == 1 and page["layout"] == "strips":
        return "thumbnail"
    if first and page["width"] and first["width"]:
        return f"1/{round(first['width'] / page['width'])} resolution"
    return "reduced resolution"


def tiff_header(reader):
    """Structure of a TIFF/BigTIFF from its directories; pixels never read."""
    head = reader.read(0, 16)
    if head[:4] not in TIFF_MAGIC:
        raise HeaderError("not a TIFF (signature missing)")
    endian = "<" if head[:2] == b"II" else ">"
    big = head[2:4] in (b"+\x00", b"\x00+")
    (first,) = struct.unpack(endian + ("Q" if big else "I"), head[8:16] if big else head[4:8])

    pages, descs, seen, offset = [], [], set(), first
    sub_levels = []
    while offset and len(pages) < MAX_PAGES:
        if offset in seen:
            raise HeaderError("directory chain loops back on itself: damaged file")
        seen.add(offset)
        tags, offset = _ifd(reader, offset, endian, big)
        pages.append(_page_summary(tags))
        descs.append(tags.get("description") if "description" in tags else "")
        if len(pages) == 1 and tags.get("subifds"):   # OME pyramid levels live here
            for so in tags["subifds"][:MAX_PAGES]:
                st, _ = _ifd(reader, so, endian, big)
                sub_levels.append(_page_summary(st))
    truncated = bool(offset)

    desc0 = descs[0] if descs else ""
    kind = _description_kind(desc0)
    out = {"container": "BigTIFF" if big else "TIFF",
           "byte_order": "little-endian" if endian == "<" else "big-endian",
           "n_pages": len(pages), "pages_truncated": truncated,
           "description": kind, "pages": pages[:50], "sub_resolutions": sub_levels[:50]}
    if kind == "Aperio scanner notes":
        out["aperio"] = _aperio(desc0)
        for i, p in enumerate(pages):
            p["role"] = _svs_role(i, p, descs[i], pages[0])
    elif kind == "OME-XML":
        out["ome"] = _ome(desc0)
        series = out["ome"].get("series", [])
        if len(series) == 1:
            s = series[0]["size"]
            names = series[0]["channels"]
            order = series[0]["dimension_order"] or ""
            # One page per channel (XYC...), as OME writers store a 1-plane image
            if s.get("Z", 1) == s.get("T", 1) == 1 and order.startswith("XYC") and \
                    s.get("C") == len(pages) == len(names):
                for p, name in zip(pages, names):
                    p["role"] = f"channel: {name or '(unnamed)'}"
    return out


def tiff_markdown(img):
    lines = ["", "## Image header",
             f"- **Container:** {img['container']} ({img['byte_order']})",
             f"- **Pages (images stored in the file):** {img['n_pages']}"
             + (" (listing stopped at the limit)" if img["pages_truncated"] else ""),
             f"- **Description text:** {img['description']}"
             + ("" if img["description"] in ("empty", "OME-XML") else " (not shown)")]
    a = img.get("aperio")
    if a:
        if "magnification" in a:
            lines.append(f"- **Scan magnification:** {a['magnification']:g}×")
        if "microns_per_pixel" in a:
            lines.append(f"- **Pixel size:** {a['microns_per_pixel']:g} µm per pixel")
    for k, s in enumerate((img.get("ome") or {}).get("series", [])[:5]):
        size = " × ".join(f"{k2}={v:,}" for k2, v in s["size"].items())
        lines.append(f"- **OME series {k}:** {size}; {s['type']}; order {s['dimension_order']}")
        if s["pixel_size"]:
            lines.append("  - pixel size: " + ", ".join(f"{k2} {v}" for k2, v in
                                                         s["pixel_size"].items()))
        if s["channels"]:
            more = f" (+{s['n_channels'] - len(s['channels'])} more)" \
                if s["n_channels"] > len(s["channels"]) else ""
            lines.append(f"  - {s['n_channels']} channel(s): " + ", ".join(
                c or "(unnamed)" for c in s["channels"]) + more)
    rows = [("page", img["pages"])] + ([("sub", img["sub_resolutions"])]
                                        if img["sub_resolutions"] else [])
    lines += ["", "| page | role | width × height | channels | bits | type | compression | layout |",
              "|---:|---|---|---:|---|---|---|---|"]
    for label, pages in rows:
        for i, p in enumerate(pages):
            if p.get("role"):
                role = p["role"]
            elif p["mask"]:
                role = "mask"
            elif label == "sub" and (p["width"], p["height"]) == \
                    (img["pages"][0]["width"], img["pages"][0]["height"]):
                role = "same-size companion image"
            elif p["reduced_resolution"] or label == "sub":
                role = "reduced resolution"
            else:
                role = ""
            wh = f"{p['width']:,} × {p['height']:,}" if p["width"] and p["height"] else "?"
            name = f"{i}" if label == "page" else f"0.{i + 1}"
            lines.append(f"| {name} | {role} | {wh} | {p['channels']} | {p['bits']} | "
                         f"{p['sample_format']} | {p['compression']} | {p['layout']} |")
    lines += ["", "### How to read this",
              "- Only the file's directories were read: sizes, bit depth, compression, "
              "tiling. No pixel was decoded or read.",
              "- *bits* and *type* together give the pixel's number range: 8-bit unsigned = "
              "0–255; 16-bit unsigned = 0–65,535 (typical for fluorescence).",
              "- Several pages of shrinking size form a pyramid: viewers show the small "
              "ones when zoomed out. Analyses use the full-resolution page.",
              "- Description text can contain file names, dates or scanner IDs, so only "
              "magnification, pixel size, dimensions and channel names are taken from it."]
    return lines


# ---- deep zoom ------------------------------------------------------------

def dzi(text):
    """Tile pyramid geometry from a .dzi descriptor (a few hundred bytes of XML)."""
    import xml.etree.ElementTree as ET
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise HeaderError(f"not valid XML ({exc})")
    if root.tag.rsplit("}", 1)[-1] != "Image":
        raise HeaderError("not a deep-zoom descriptor (no <Image> root)")
    size = next((e for e in root if e.tag.rsplit("}", 1)[-1] == "Size"), None)
    if size is None:
        raise HeaderError("deep-zoom descriptor without <Size>")
    w, h = int(size.attrib["Width"]), int(size.attrib["Height"])
    tile, overlap = int(root.attrib["TileSize"]), int(root.attrib.get("Overlap", 0))
    levels = math.ceil(math.log2(max(w, h))) + 1 if max(w, h) > 1 else 1
    tiles = 0
    for lv in range(levels):
        scale = 2 ** (levels - 1 - lv)
        lw, lh = max(1, math.ceil(w / scale)), max(1, math.ceil(h / scale))
        tiles += math.ceil(lw / tile) * math.ceil(lh / tile)
    return {"width": w, "height": h, "tile_size": tile, "overlap": overlap,
            "tile_format": root.attrib.get("Format"), "levels": levels,
            "expected_tiles": tiles}


def dzi_markdown(d):
    lines = ["", "## Deep-zoom pyramid",
             f"- **Full image:** {d['width']:,} × {d['height']:,} pixels",
             f"- **Tiles:** {d['tile_size']} px ({d['overlap']} px overlap), "
             f"{d['tile_format']} format, {d['levels']} zoom levels",
             f"- **Tiles the folder should hold:** {d['expected_tiles']:,}"]
    if "catalog_tiles" in d:
        n = d["catalog_tiles"]
        e = d["expected_tiles"]
        if n == e:
            verdict = "matches"
        elif n == e + 1:
            # libvips (the usual tiling tool) also writes vips-properties.xml
            # into the folder: one object that is not a tile.
            verdict = "matches, plus one non-tile object (libvips writes vips-properties.xml)"
        elif n < e:
            verdict = f"FEWER than expected: {e - n:,} missing"
        else:
            verdict = f"{n - e:,} more objects than tiles"
        lines.append(f"- **Tiles the catalog lists in `{d['tile_folder']}`:** {n:,} ({verdict})")
    elif "tile_folder" in d:
        lines.append(f"- **Tile folder `{d['tile_folder']}`:** not in the catalog")
    lines += ["", "### How to read this",
              "- A deep-zoom pyramid stores a huge image as small tiles at every zoom "
              "level, for viewing in a browser. It is a display copy: analyses should use "
              "the original image, not the tiles.",
              "- The expected tile count follows from the image and tile sizes; a folder "
              "with fewer tiles has gaps when zoomed in."]
    return lines
