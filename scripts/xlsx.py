"""
xlsx.py — the structure of an Excel workbook, read with the standard library only.

Not a step. Used by inspect_data.py for .xlsx files. An .xlsx file is a zip
archive of XML files: xl/workbook.xml names the sheets, xl/sharedStrings.xml
holds the text, and xl/worksheets/sheetN.xml holds the cells. Reading those
directly needs no new package (openpyxl is not installed; adding it is your
decision), and every value goes straight into the same column rules as CSV
and h5ad (inspect_data.describe_column), so the withholding rules are the
same in every format.

Decisions made here, and why:

1. THE HEADER IS THE FIRST NEARLY-FULL ROW, AND THE REPORT SAYS WHICH ROW.
   Spreadsheets often have a title line or merged group labels above the
   column names; those rows fill fewer cells. The header is the first row
   filling at least 80% as many cells as the fullest of the first 10 rows.
   The report states the row used and how many rows above it were skipped,
   so a wrong guess is visible rather than silent.

2. TYPES COME FROM THE CELLS, WITH ONE CAVEAT STATED.
   Numbers, text, true/false and Excel errors are told apart by the cell
   type. Excel stores DATES as numbers with a display format; this reader
   does not decode display formats, so a date column shows as float. The
   report says so.

3. LIMITS AGAINST A HOSTILE ARCHIVE.
   At most 200 MB unpacked, 2,000 archive entries and 100 sheets; a file
   beyond that is reported, not read. Formulas are never evaluated: only
   the value Excel saved with the file is read.

4. HIDDEN SHEETS ARE NAMED AS HIDDEN.
   A workbook can hide sheets; their structure is reported like any other
   sheet, marked hidden, because hidden is not the same as empty.
"""

import re
import zipfile
import xml.etree.ElementTree as ET

MAX_UNPACKED = 200 * 1024 * 1024
MAX_ENTRIES = 2000
MAX_SHEETS = 100
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "rel": "http://schemas.openxmlformats.org/package/2006/relationships"}


class WorkbookError(ValueError):
    pass


def _col(ref):
    """'C12' -> 2 (zero-based column)."""
    letters = re.match(r"[A-Z]+", ref or "")
    if not letters:
        return None
    n = 0
    for ch in letters.group(0):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def _text(el):
    """All text inside an <si> or <is> element (rich-text runs joined)."""
    return "".join(t.text or "" for t in el.iter(f"{{{NS['m']}}}t"))


def _check_archive(z):
    infos = z.infolist()
    if len(infos) > MAX_ENTRIES:
        raise WorkbookError(f"{len(infos):,} archive entries (limit {MAX_ENTRIES:,})")
    unpacked = sum(i.file_size for i in infos)
    if unpacked > MAX_UNPACKED:
        raise WorkbookError(f"unpacks to {unpacked / 1e6:,.0f} MB (limit "
                            f"{MAX_UNPACKED // 1024 // 1024} MB)")


def _sheets(z):
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels.findall("rel:Relationship", NS)}
    out = []
    for s in wb.findall("m:sheets/m:sheet", NS)[:MAX_SHEETS]:
        t = target.get(s.get(f"{{{NS['r']}}}id"), "")
        path = t.lstrip("/") if t.startswith("/") else "xl/" + t
        out.append({"name": s.get("name"), "path": path,
                    "hidden": s.get("state") in ("hidden", "veryHidden")})
    return out


def _shared(z):
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    return [_text(si) for si in root.findall("m:si", NS)]


def _cells(z, path, shared):
    """Yield (row number, {column: (kind, value)}) for every row with a cell."""
    with z.open(path) as fh:
        for _, el in ET.iterparse(fh):
            if el.tag != f"{{{NS['m']}}}row":
                continue
            row = {}
            for c in el.findall("m:c", NS):
                j = _col(c.get("r"))
                t = c.get("t", "n")
                v = c.find("m:v", NS)
                if t == "inlineStr":
                    is_ = c.find("m:is", NS)
                    val = ("text", _text(is_)) if is_ is not None else None
                elif v is None or v.text is None:
                    val = None
                elif t == "s":
                    val = ("text", shared[int(v.text)] if int(v.text) < len(shared) else "")
                elif t in ("str",):
                    val = ("text", v.text)
                elif t == "b":
                    val = ("boolean", v.text == "1")
                elif t == "e":
                    val = ("error", v.text)
                else:
                    try:
                        f = float(v.text)
                        val = ("number", int(f) if f.is_integer() else f)
                    except ValueError:
                        val = ("text", v.text)
                if j is not None and val is not None and val[1] != "":
                    row[j] = val
            yield int(el.get("r") or 0), row
            el.clear()


def read_sheet(z, sheet, shared, track_cap):
    """Header, row count and per-column facts for one sheet.

    Each column counts the cells it HAS ('filled'); missing = rows - filled,
    so a column that only starts partway down is still counted correctly.
    """
    header, header_row, cols, n, above = None, None, {}, 0, 0
    rows_iter = _cells(z, sheet["path"], shared)
    # Look at the first rows with content before choosing the header: title or
    # merged group-label rows ("RIPA Buffer Tissue extractions") sit above the
    # real column names and fill fewer cells. The header is the first row that
    # fills at least 80% as many cells as the fullest of the first 10 rows.
    first = []
    for r, cells in rows_iter:
        if cells:
            first.append((r, cells))
        if len(first) == 10:
            break
    if first:
        widest = max(len(c) for _, c in first)
        k = next(i for i, (_, c) in enumerate(first) if len(c) >= 0.8 * widest)
        r, cells = first[k]
        width = max(cells) + 1
        header = [str(cells[j][1]) if j in cells else f"(column {j + 1})" for j in range(width)]
        header_row, above = r, k
        pending = first[k + 1:]
        header = _disambiguate(header, [c for _, c in first[:k]])
    else:
        pending = []

    def body():
        yield from pending
        yield from rows_iter
    for r, cells in body():
        if not cells:                 # a formatted but blank row is not a data row
            continue
        n += 1
        for j, (kind, val) in cells.items():
            c = cols.setdefault(j, {"kinds": set(), "filled": 0, "distinct": set(),
                                    "overflow": False, "first": []})
            c["filled"] += 1
            c["kinds"].add("integer" if kind == "number" and isinstance(val, int)
                           else "float" if kind == "number" else kind)
            if len(c["first"]) < 200 and kind == "text":
                c["first"].append(val)
            if not c["overflow"]:
                c["distinct"].add(val)
                if len(c["distinct"]) > track_cap:
                    c["overflow"], c["distinct"] = True, set()
    header = header or []
    for j in sorted(cols):                       # cells right of the header row
        while j >= len(header):
            header.append(f"(column {len(header) + 1})")
    return header, header_row, cols, n, above


def _disambiguate(header, label_rows):
    """Repeated column names get the group label above them, e.g. two
    'ABeta40 pg/ug' columns under 'RIPA ...' and 'GuHCl ...'. A merged label
    cell sits in its first column, so each label carries rightwards until the
    next one."""
    counts = {}
    for h in header:
        counts[h] = counts.get(h, 0) + 1
    if not label_rows or all(n == 1 for n in counts.values()):
        return header
    labels = [""] * len(header)
    for row in label_rows:
        current = ""
        for j in range(len(header)):
            if j in row:
                current = str(row[j][1])
            if current:
                labels[j] = current
    return [f"{h} [{labels[j]}]" if counts[h] > 1 and labels[j] else h
            for j, h in enumerate(header)]


def _kind(kinds):
    if not kinds:
        return None
    if kinds <= {"integer"}:
        return "integer"
    if kinds <= {"integer", "float"}:
        return "float"
    if kinds == {"boolean"}:
        return "boolean"
    return "text"


def inspect(path, describe_column, looks_like_ids, track_cap=100_000):
    """{'sheets': [...]} with each sheet's columns described by the shared rules."""
    try:
        z = zipfile.ZipFile(path)
    except zipfile.BadZipFile as exc:
        raise WorkbookError(f"not a readable zip archive ({exc})")
    with z:
        _check_archive(z)
        if "xl/workbook.xml" not in z.namelist():
            raise WorkbookError("zip archive without xl/workbook.xml: not an Excel workbook")
        shared = _shared(z)
        sheets = []
        for s in _sheets(z):
            header, header_row, cols, n, above = read_sheet(z, s, shared, track_cap)
            columns = []
            for j, name in enumerate(header):
                c = cols.get(j, {"kinds": set(), "filled": 0, "distinct": set(),
                                 "overflow": False, "first": []})
                columns.append(describe_column(
                    name, _kind(c["kinds"]), n, n - c["filled"],
                    None if c["overflow"] else c["distinct"], looks_like_ids(c["first"])))
            sheets.append({"name": s["name"], "hidden": s["hidden"], "rows": n,
                           "header_row": header_row, "rows_above_header": above,
                           "columns": columns})
    return {"sheets": sheets}
