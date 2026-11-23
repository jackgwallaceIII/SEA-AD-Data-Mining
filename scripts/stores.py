"""
stores.py — describe a Zarr store or a Xenium bundle from its descriptor files alone.

Not a step. Used by inspect_data.py when the file being inspected is a
consolidated Zarr metadata file (zmetadata), a folder of Zarr metadata, an
experiment.xenium or a gene_panel.json — the files `get --descriptors`
fetches. No chunk, parquet part or image is needed or read.

Decisions made here, and why:

1. A SpatialData STORE IS DESCRIBED ELEMENT BY ELEMENT.
   SpatialData keeps five kinds of element in fixed groups: images
   (microscope pictures), labels (segmentation masks: which pixel belongs
   to which cell), points (transcripts), shapes (cell outlines) and tables
   (a cell-by-gene AnnData). The report lists each element with its size,
   number type and storage, so the next person knows what is there before
   downloading gigabytes.

2. COORDINATE SYSTEMS AND SCALES ARE THE POINT.
   Each element states how its pixels map to a named coordinate system
   (e.g. 'global', in micrometres). Pixel size, and whether elements share
   a coordinate system, decide whether an image, its masks and its
   transcripts line up. The report shows the scale factors and the type
   of each transformation (identity, scale, affine); affine matrices are
   named, not printed.

3. TABLES: NAMES, TYPES AND SIZES ONLY.
   The table's metadata gives its shape, its obs column names and types,
   and for category columns the NUMBER of categories (from the shape of
   the categories array). The categories themselves live in chunks that
   are never fetched, so no level, donor ID or barcode can appear.

4. XENIUM RUN FILES: A WHITELIST.
   experiment.xenium mixes structure (pixel size, panel, software
   versions, cell count) with run identifiers (instrument serial number,
   slide ID, run name, UUIDs, times). Only whitelisted fields are shown;
   the names of the others are listed so the reader knows they exist.
   gene_panel.json: panel name, organism, number of targets by type.
   Gene names are counted, not listed.
"""

import json
import pathlib

ELEMENT_GROUPS = ("images", "labels", "points", "shapes", "tables")
XENIUM_SHOWN = ("major_version", "minor_version", "patch_version", "analysis_sw_version",
                "instrument_sw_version", "chemistry_version", "preservation_method",
                "panel_name", "panel_organism", "panel_type", "panel_tissue_type",
                "panel_num_targets_predesigned", "panel_num_targets_custom",
                "segmentation_stain", "pixel_size", "z_step_size", "num_cells")
MAX_JSON = 50_000_000


def load_json(path):
    path = pathlib.Path(path)
    if path.stat().st_size > MAX_JSON:
        raise ValueError(f"{path.name} is over {MAX_JSON // 1_000_000} MB: not a descriptor")
    return json.loads(path.read_text(encoding="utf-8"))


# ---- Zarr ---------------------------------------------------------------------

def metadata_from_consolidated(path):
    d = load_json(path)
    if "metadata" not in d:
        raise ValueError("no 'metadata' key: not consolidated Zarr metadata")
    return d["metadata"]


def metadata_from_folder(folder):
    """{relative path: json} for every metadata file under a folder."""
    folder = pathlib.Path(folder)
    out = {}
    for p in folder.rglob("*"):
        if p.is_file() and p.name in (".zgroup", ".zarray", ".zattrs", "zarr.json"):
            out[p.relative_to(folder).as_posix()] = load_json(p)
    return out


def _dtype(z):
    dt = z.get("dtype", z.get("data_type", "?"))
    if isinstance(dt, str) and dt in ("|O", "object") or \
            any(f.get("id") == "vlen-utf8" for f in (z.get("filters") or []) if isinstance(f, dict)):
        return "text"
    names = {"u1": "uint8", "u2": "uint16", "u4": "uint32", "u8": "uint64", "i1": "int8",
             "i2": "int16", "i4": "int32", "i8": "int64", "f2": "float16", "f4": "float32",
             "f8": "float64", "b1": "bool"}
    return names.get(str(dt).lstrip("<>|="), str(dt))


def _compressor(z):
    c = z.get("compressor")
    if isinstance(c, dict):
        return c.get("cname") or c.get("id") or "?"
    return "none" if c is None else str(c)


def _transforms(attrs):
    """[(coordinate system, kind, scale or None)] from SpatialData / NGFF attrs."""
    found = []
    ms = (attrs.get("multiscales") or [{}])[0] if attrs.get("multiscales") else None
    lists = [attrs.get("coordinateTransformations") or []]
    if ms:
        lists.append(ms.get("coordinateTransformations") or [])
    for lst in lists:
        for t in lst:
            out = t.get("output")
            system = out.get("name") if isinstance(out, dict) else out
            kind = t.get("type", "?")
            scale = t.get("scale") if kind == "scale" else None
            if kind == "sequence":
                kinds = [x.get("type") for x in t.get("transformations", [])]
                scales = [x.get("scale") for x in t.get("transformations", []) if x.get("scale")]
                kind = "sequence: " + " → ".join(kinds)
                scale = scales[0] if scales else None
            found.append((system or "?", kind, scale))
    return found


def _raster(meta, base, attrs):
    """An image or labels element: pyramid levels, axes, channels, transforms."""
    ms = (attrs.get("multiscales") or [{}])[0]
    axes = [a.get("name", "?") + (f" ({a['unit']})" if a.get("unit") else "")
            for a in ms.get("axes", [])]
    levels = []
    for ds in ms.get("datasets", []):
        z = meta.get(f"{base}/{ds.get('path')}/.zarray") or {}
        scale = next((t.get("scale") for t in ds.get("coordinateTransformations", [])
                      if t.get("type") == "scale"), None)
        levels.append({"path": ds.get("path"), "shape": z.get("shape"), "dtype": _dtype(z),
                       "chunks": z.get("chunks"), "compressor": _compressor(z),
                       "scale": scale})
    channels = [c.get("label", "") for c in (attrs.get("omero") or {}).get("channels", [])]
    return {"axes": axes, "levels": levels, "channels": channels,
            "transforms": _transforms(attrs), "ngff_version": ms.get("version")}


def _column(meta, base, name):
    attrs = meta.get(f"{base}/{name}/.zattrs") or {}
    enc = attrs.get("encoding-type", "")
    if enc == "categorical":
        cats = meta.get(f"{base}/{name}/categories/.zarray") or {}
        n = (cats.get("shape") or [None])[0]
        return {"name": name, "type": f"category ({n} levels)" if n is not None else "category"}
    if enc.startswith("nullable"):
        z = meta.get(f"{base}/{name}/values/.zarray") or {}
        return {"name": name, "type": f"{_dtype(z)} (nullable)"}
    z = meta.get(f"{base}/{name}/.zarray") or {}
    return {"name": name, "type": _dtype(z) if z else (enc or "?")}


def _frame(meta, base):
    attrs = meta.get(f"{base}/.zattrs") or {}
    index = attrs.get("_index", "_index")
    idx = meta.get(f"{base}/{index}/.zarray") or meta.get(f"{base}/{index}/values/.zarray") or {}
    cols = attrs.get("column-order") or []
    if isinstance(cols, str):
        cols = [cols]
    return {"rows": (idx.get("shape") or [None])[0], "index": index,
            "columns": [_column(meta, base, c) for c in cols]}


def _matrix(meta, base):
    attrs = meta.get(f"{base}/.zattrs") or {}
    enc = attrs.get("encoding-type")
    if enc in ("csr_matrix", "csc_matrix"):
        data = meta.get(f"{base}/data/.zarray") or {}
        return {"storage": f"{enc[:3]} sparse", "shape": attrs.get("shape"),
                "dtype": _dtype(data), "stored": (data.get("shape") or [None])[0]}
    z = meta.get(f"{base}/.zarray")
    if z:
        return {"storage": "dense", "shape": z.get("shape"), "dtype": _dtype(z), "stored": None}
    return {"storage": enc or "?", "shape": None, "dtype": "", "stored": None}


def _children(meta, base):
    """Immediate child names of a group, from the metadata paths."""
    pre = f"{base}/" if base else ""
    names = set()
    for k in meta:
        if k.startswith(pre):
            rest = k[len(pre):]
            if "/" in rest:
                names.add(rest.split("/", 1)[0])
    return sorted(names)


def _table(meta, base):
    x = _matrix(meta, f"{base}/X")
    obs, var = _frame(meta, f"{base}/obs"), _frame(meta, f"{base}/var")
    return {"X": x, "n_obs": obs["rows"] or (x["shape"] or [None])[0],
            "n_vars": var["rows"] or (x["shape"] or [None, None])[1],
            "layers": {n: _matrix(meta, f"{base}/layers/{n}")
                       for n in _children(meta, f"{base}/layers")},
            "obs": obs, "var": var,
            "obsm": {n: (_matrix(meta, f"{base}/obsm/{n}")["shape"])
                     for n in _children(meta, f"{base}/obsm")},
            "uns": _children(meta, f"{base}/uns")}


def describe_zarr(meta, source):
    root = meta.get(".zattrs") or {}
    sd = root.get("spatialdata_attrs") or {}
    out = {"source": source, "zarr_format": 3 if any(k.endswith("zarr.json") for k in meta)
           else 2, "spatialdata_version": sd.get("version"), "elements": {},
           "metadata_entries": len(meta)}
    systems = set()
    for group in ELEMENT_GROUPS:
        for name in _children(meta, group):
            base = f"{group}/{name}"
            attrs = meta.get(f"{base}/.zattrs") or {}
            if group in ("images", "labels"):
                el = _raster(meta, base, attrs)
            elif group == "tables":
                el = _table(meta, base)
            else:
                el = {"transforms": _transforms(attrs), "storage": "parquet",
                      "attrs": sorted(k for k in (attrs.get("spatialdata_attrs") or {}))}
                if group == "shapes":
                    el["encoding"] = attrs.get("encoding-type") or \
                        (attrs.get("spatialdata_attrs") or {}).get("geos", {}).get("name")
            for system, _, _ in el.get("transforms", []):
                systems.add(system)
            out["elements"].setdefault(group, {})[name] = el
    out["coordinate_systems"] = sorted(systems)
    if out["zarr_format"] == 3:
        # The v3 layout keeps metadata in zarr.json files with different keys;
        # the element details above are read from the v2 layout only.
        out["note"] = "Zarr v3 layout: element details are not interpreted yet"
    return out


def _shape(s):
    return " × ".join(f"{x:,}" for x in s) if s else "?"


def _scale(sc):
    return "[" + ", ".join(f"{x:.4g}" for x in sc) + "]" if sc else ""


def _unit_note(z):
    """If images map to a system unchanged, that system is in image pixels; a
    points/shapes scale s into it means 1 stored unit = s image pixels."""
    pixel_systems = {sys_ for el in z["elements"].get("images", {}).values()
                     for sys_, kind, _ in el["transforms"] if kind == "identity"}
    notes = []
    for group in ("points", "shapes"):
        for name, el in z["elements"].get(group, {}).items():
            for sys_, kind, sc in el["transforms"]:
                if sys_ in pixel_systems and kind == "scale" and sc and sc[0] > 1 \
                        and abs(sc[0] - sc[1]) < 1e-9:
                    notes.append((group, name, sys_, sc[0]))
    if not notes:
        return []
    s = notes[0][3]
    same = all(abs(n[3] - s) < 1e-9 for n in notes)
    who = ", ".join(f"{g}/{n}" for g, n, _, _ in notes)
    return [f"- **Units:** `{notes[0][2]}` is measured in image pixels (the image maps to it "
            f"unchanged). {who} {'are' if len(notes) > 1 else 'is'} stored in a unit "
            f"{s:.4g}× larger" + (": 1 ÷ " + f"{s:.4g} = {1 / s:.4g}. If the image pixel is "
            f"{1 / s:.4g} µm (Xenium's standard is 0.2125 µm), these are stored in "
            f"micrometres." if same else ".")]


def zarr_markdown(z):
    lines = ["", "## Zarr store",
             f"- **Described from:** {z['source']} ({z['metadata_entries']:,} metadata entries; "
             f"no data chunks read)",
             f"- **Zarr format:** v{z['zarr_format']}"
             + (f"; SpatialData {z['spatialdata_version']}" if z["spatialdata_version"] else ""),
             f"- **Coordinate systems:** {', '.join(z['coordinate_systems']) or 'none declared'}",
             ] + ([f"- **Note:** {z['note']}"] if z.get("note") else []) + [
             "- **Elements:** " + "; ".join(f"{g} {len(e)} ({', '.join(e)})"
                                         for g, e in z["elements"].items())] + _unit_note(z)
    if not any("(" in a for el in z["elements"].get("images", {}).values() for a in el["axes"]):
        lines.append("- **Axis units:** not declared in the metadata")
    for group in ("images", "labels"):
        for name, el in z["elements"].get(group, {}).items():
            lines += ["", f"### {group}/{name}",
                      f"- axes: {', '.join(el['axes']) or '?'}"
                      + (f"; channels: {', '.join(el['channels'])}" if el["channels"] else "")]
            for system, kind, scale in el["transforms"]:
                lines.append(f"- to `{system}`: {kind}" + (f" {_scale(scale)}" if scale else ""))
            lines += ["", "| level | shape | dtype | chunks | compression | scale |",
                      "|---|---|---|---|---|---|"]
            for lv in el["levels"]:
                lines.append(f"| {lv['path']} | {_shape(lv['shape'])} | {lv['dtype']} | "
                             f"{_shape(lv['chunks'])} | {lv['compressor']} | "
                             f"{_scale(lv['scale'])} |")
    for group in ("points", "shapes"):
        for name, el in z["elements"].get(group, {}).items():
            tr = "; ".join(f"to `{s}`: {k}" + (f" {_scale(sc)}" if sc else "")
                           for s, k, sc in el["transforms"]) or "no transform declared"
            extra = f"; geometry: {el['encoding']}" if el.get("encoding") else ""
            lines += ["", f"### {group}/{name}",
                      f"- stored as {el['storage']} (row count needs the parquet footer: "
                      f"pyarrow){extra}", f"- {tr}"]
    for name, t in z["elements"].get("tables", {}).items():
        x = t["X"]
        lines += ["", f"### tables/{name} (AnnData)",
                  f"- **{t['n_obs']:,} cells × {t['n_vars']:,} genes**"
                  if t["n_obs"] and t["n_vars"] else "- size: ?",
                  f"- X: {x['storage']}, {x['dtype']}"
                  + (f", {x['stored']:,} stored values" if x.get("stored") else ""),
                  "- layers: " + (", ".join(f"{n} ({m['storage']}, {m['dtype']})"
                                            for n, m in t["layers"].items()) or "none"),
                  "- obsm: " + (", ".join(f"{n} {_shape(s)}" for n, s in t["obsm"].items())
                                or "none"),
                  "- uns keys: " + (", ".join(t["uns"]) or "none"),
                  "", f"| obs column ({len(t['obs']['columns'])}) | type |", "|---|---|"]
        for c in t["obs"]["columns"]:
            lines.append(f"| {c['name']} | {c['type']} |")
        lines.append(f"\nvar: {len(t['var']['columns'])} column(s): "
                     + ", ".join(f"{c['name']} ({c['type']})" for c in t["var"]["columns"]))
    lines += ["", "### How to read this",
              "- Only metadata was read: names, shapes, number types, chunking and "
              "coordinate transformations. No chunk (pixel, transcript or cell value) was "
              "fetched, so no category level or identifier can appear here.",
              "- Each element states how its own coordinates map into a named coordinate "
              "system: *identity* = unchanged; *scale s* = multiply by s; *affine* = a "
              "registration (shift, rotation, scaling). Elements mapped into the same "
              "system line up. In the level tables, *scale* is relative to level 0: each "
              "pyramid level halves the resolution, so the factor doubles.",
              "- Coordinates in different units (image pixels vs micrometres) only line up "
              "through these factors; ignoring them misplaces transcripts by the factor.",
              "- Table values (X vs layers/counts) were not read, so whether X is normalised "
              "cannot be judged from metadata alone.",
              "- *labels* are segmentation masks: each pixel holds the ID of the cell (or "
              "nucleus) it belongs to, 0 for background.",
              "- *category (N levels)* gives the number of categories from the metadata; "
              "the categories themselves were not read."]
    return lines


# ---- Xenium -------------------------------------------------------------------

def describe_xenium(path):
    d = load_json(path)
    shown = {k: d[k] for k in XENIUM_SHOWN if k in d and not isinstance(d[k], (dict, list))}
    ranger = d.get("xenium_ranger")
    if isinstance(ranger, dict) and isinstance(ranger.get("version"), str):
        # Re-processed with Xenium Ranger (e.g. new cell segmentation). Only the
        # version is shown: its command line can contain paths and sample names.
        shown["xenium_ranger_version"] = ranger["version"]
    nested = {k: sorted(v) for k, v in d.items() if isinstance(v, dict)}
    withheld = sorted(k for k in d if k not in shown and k not in nested)
    return {"fields": shown, "file_groups": nested, "withheld_fields": withheld}


def describe_panel(path):
    """Named fields only. The file also holds the customer's and the panel
    designer's names and e-mail addresses (payload.customer, payload.designer);
    those keys are never read into the report."""
    d = load_json(path)
    payload = d.get("payload", d)
    panel = payload.get("panel") or d.get("panel") or {}
    ident = panel.get("identity") or {}
    targets = payload.get("targets") or []
    kinds = {}
    for t in targets:
        typ = t.get("type") or {}
        kind = typ.get("descriptor") or typ.get("type") or "?"
        kinds[kind] = kinds.get(kind, 0) + 1
    return {"panel_name": ident.get("name") or panel.get("name"),
            "panel_version": ident.get("version"),
            "organism": panel.get("species") or panel.get("organism"),
            "tissue": panel.get("tissue"),
            "panel_type": (panel.get("type") or {}).get("descriptor"),
            "declared_gene_targets": panel.get("num_gene_targets"),
            "chemistry_version": (payload.get("chemistry") or {}).get("version"),
            "spec_version": payload.get("spec_version"),
            "n_targets": len(targets), "targets_by_type": kinds}


def xenium_markdown(x):
    lines = ["", "## Xenium run settings (experiment.xenium)"]
    for k, v in x["fields"].items():
        lines.append(f"- **{k}:** {v}")
    for k, v in x["file_groups"].items():
        lines.append(f"- **{k}:** {', '.join(v)}")
    if x["withheld_fields"]:
        lines.append(f"- *not shown (run identifiers, times, results):* "
                     f"{', '.join(x['withheld_fields'])}")
    lines += ["", "### How to read this",
              "- *pixel_size* is micrometres per pixel of the Xenium morphology images; "
              "transcript and cell coordinates in the bundle are in micrometres.",
              "- Fields that identify the run, slide or instrument are deliberately not "
              "shown."]
    return lines


def panel_markdown(p):
    kinds = ", ".join(f"{n:,} {k}" for k, n in sorted(p["targets_by_type"].items()))
    extra = ", ".join(f"{label} {p[k]}" for k, label in
                      (("panel_version", "version"), ("tissue", "tissue"),
                       ("panel_type", "type"), ("chemistry_version", "chemistry"))
                      if p.get(k))
    genes = p["targets_by_type"].get("gene")
    check = ""
    if p.get("declared_gene_targets") is not None and genes is not None:
        check = (" (matches the declared gene count)" if genes == p["declared_gene_targets"]
                 else f" (the panel declares {p['declared_gene_targets']:,} gene targets)")
    return ["", "## Xenium gene panel (gene_panel.json)",
            f"- **Panel:** {p['panel_name'] or '?'} ({p['organism'] or 'organism ?'})"
            + (f"; {extra}" if extra else ""),
            f"- **Targets:** {p['n_targets']:,}" + (f" — {kinds}" if kinds else "") + check,
            "- *not read:* customer and designer names and e-mail addresses",
            "", "### How to read this",
            "- A Xenium panel measures a fixed list of genes (hundreds to ~5,000), not the "
            "whole transcriptome; genes outside the panel are simply absent, not zero.",
            "- Control targets (negative probes, codewords) estimate background; they are "
            "not genes."]
