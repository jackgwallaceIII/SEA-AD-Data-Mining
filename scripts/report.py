"""
report.py — one readable data-structure report per brain region, for the next person.

Not a step. Run through index.py:

    python scripts/index.py report LEC          # results/reports/LEC.md
    python scripts/index.py report --all        # every region + results/reports/README.md

Built ONLY from what is already on disk: the catalog (catalog/files/),
the per-file structure reports (results/inspect/**.json) and the known
upstream issues (catalog/known_issues.csv). It never contacts S3 and
never opens a data file, so it can be regenerated at any time and gives
the same output for the same inputs.

Decisions made here, and why:

1. EVERY STATEMENT SAYS WHERE IT COMES FROM.
   Three sources, kept visibly apart:
   - the catalog: what EXISTS (names, sizes, counts), from listing the
     buckets on the date shown;
   - inspected files: what is INSIDE, checked on specific files (named,
     with the date and how they were read);
   - known issues: problems found upstream, with their evidence.
   A reader can always tell a measured fact from a catalog label.

2. "CHECKED ON 1 OF 86" IS SAID OUT LOUD.
   Structure is usually verified on one or two example files per kind
   (one per-donor object stands for 86). The report gives that ratio, so
   nobody mistakes an example for a census. Kinds of file with no
   inspected example are listed as NOT YET INSPECTED rather than left out.

3. NOTHING NEW IS REVEALED.
   The report only summarises the inspect reports, which already applied
   the withholding rules (identifiers, few-donor files, scanner notes, run
   IDs). It adds no category levels at all: shapes, counts, verdicts and
   links to the per-file reports, where the details live.

4. REGION-ONLY AND SHARED ARE SEPARATE, AS IN `show`.
   Units covering several regions (the multi-region release, region-pair
   neuropathology folders) are reported in their own section, so a size
   total for "LEC" means LEC only.

5. DETERMINISTIC.
   Sorted throughout; the only dates in the output are the catalog's
   listing dates, not the time of running. Re-running without new inputs
   changes nothing, so a diff of the report shows real changes only.
"""

import collections
import json
import pathlib
import re

import browse
import index
import issues as issues_mod
import plan
import regions as regions_mod

ROOT = pathlib.Path(__file__).resolve().parent.parent
REPORT_DIR = ROOT / "results" / "reports"
INSPECT_DIR = ROOT / "results" / "inspect"


# ---- inputs -------------------------------------------------------------------

def load_inspections(inspect_dir=INSPECT_DIR):
    """Every per-file report JSON under results/inspect/."""
    out = []
    for p in sorted(pathlib.Path(inspect_dir).rglob("*.json")):
        try:
            r = json.loads(p.read_text())
        except (ValueError, OSError):
            continue
        if r.get("bucket") and r.get("key"):
            r["_json"] = p
            out.append(r)
    return out


def unit_for(rows_by_key, prefixes, bucket, key):
    """The catalog unit holding a file: exact key, or the folder unit around it."""
    row = rows_by_key.get((bucket, key))
    if row:
        return row
    for pre, r in prefixes.get(bucket, []):
        if key.startswith(pre):
            return r
    return None


def attach(rows, inspections):
    """{(bucket, unit): [inspection, ...]} for inspections that map to a unit."""
    by_key = {(r["bucket"], r["unit"]): r for r in rows}
    prefixes = collections.defaultdict(list)
    for r in rows:
        if r["kind"] == "collapsed":
            prefixes[r["bucket"]].append((r["unit"], r))
    for b in prefixes:
        prefixes[b].sort(key=lambda x: -len(x[0]))         # deepest folder first
    found = collections.defaultdict(list)
    for rep in inspections:
        row = unit_for(by_key, prefixes, rep["bucket"], rep["key"])
        if row:
            found[(row["bucket"], row["unit"])].append(rep)
    return found


# ---- one line per inspected file ------------------------------------------------

def _n(x):
    return f"{x:,}" if isinstance(x, int) else str(x)


VERDICTS = (("only 0/1", "binarised (0/1)"), ("negative", "scaled/centred"),
            ("integer type", "counts (integers)"), ("whole numbers only", "raw counts"),
            ("non-whole", "normalised/transformed"))


def verdict(text):
    """The inspect verdict in two or three words."""
    return next((short for key, short in VERDICTS if text.startswith(key)), text or "?")


def summarise(rep):
    """A one-line structural summary of an inspect report."""
    if "anndata" in rep:
        a = rep["anndata"]
        feat = (a.get("var", {}).get("index") or {}).get("looks_like", "features")
        feat = "peaks" if "peaks" in feat else "genes" if "gene" in feat else "features"
        mats = "; ".join(f"{m['name']}: {verdict(m.get('values', ''))}"
                         for m in a.get("matrices", []))
        return f"{_n(a.get('n_obs'))} cells × {_n(a.get('n_vars'))} {feat}; {mats}"
    if "columns" in rep:
        cols = rep["columns"]
        empty = sum(1 for c in cols if c.get("empty"))
        held = sum(1 for c in cols if c.get("levels_withheld"))
        return (f"table, {_n(rep.get('rows'))} rows × {len(cols)} columns; {empty} empty, "
                f"{held} withheld")
    if "image" in rep:
        im = rep["image"]
        p0 = (im.get("pages") or [{}])[0]
        bits = f"{p0.get('bits')}-bit {p0.get('sample_format', '')}".strip()
        comp = "uncompressed" if p0.get("compression") == "none" else p0.get("compression")
        s = (f"{im['container']}, {im['n_pages']} page(s); full {_n(p0.get('width'))} × "
             f"{_n(p0.get('height'))}, {p0.get('channels')} ch, {bits}, {comp}")
        ap, ome = im.get("aperio") or {}, (im.get("ome") or {}).get("series") or []
        if ap:
            s += f"; {ap.get('magnification', '?'):g}×, {ap.get('microns_per_pixel', '?')} µm/px"
        if ome:
            ch = ome[0].get("channels") or []
            px = (ome[0].get("pixel_size") or {}).get("X")
            if px:                           # older reports stored unrounded numbers
                num, _, unit = px.partition(" ")
                try:
                    px = f"{float(num):.4g} {unit}".strip()
                except ValueError:
                    pass
            s += f"; channels {', '.join(ch)}" + (f"; {px}/px" if px else "")
        return s
    if "deep_zoom" in rep:
        d = rep["deep_zoom"]
        ok = ""
        if "catalog_tiles" in d:
            ok = ("; tiles complete" if d["catalog_tiles"] in (d["expected_tiles"],
                                                                d["expected_tiles"] + 1)
                  else f"; {d['catalog_tiles']:,} of {d['expected_tiles']:,} tiles")
        return f"deep-zoom {d['width']:,} × {d['height']:,}, {d['levels']} levels{ok}"
    if "track" in rep:
        t = rep["track"]
        return (f"bigWig, {t['genome_build']}, {t['chromosomes']['n']} chromosomes/contigs")
    if "zarr" in rep:
        z = rep["zarr"]
        parts = [f"{g} {len(e)}" for g, e in z["elements"].items()]
        tab = next(iter(z["elements"].get("tables", {}).values()), None)
        t = f"; table {_n(tab['n_obs'])} cells × {_n(tab['n_vars'])} genes" if tab else ""
        return (f"SpatialData: {', '.join(parts)}{t}; coordinate systems "
                f"{', '.join(z['coordinate_systems']) or 'none'}")
    if "xenium" in rep:
        f = rep["xenium"]["fields"]
        return (f"Xenium run: panel {f.get('panel_name', '?')}, pixel {f.get('pixel_size', '?')} "
                f"µm, {_n(f.get('num_cells', '?'))} cells")
    if "panel" in rep:
        p = rep["panel"]
        kinds = ", ".join(f"{n} {k}" for k, n in sorted(p["targets_by_type"].items()))
        return f"gene panel {p.get('panel_name') or '?'}: {kinds}"
    if "text_lines" in rep:
        return f"plain text, {rep['text_lines']:,} lines"
    return rep.get("status", "?")


def how_read(rep):
    s = rep.get("status", "")
    m = re.search(r"(\d{4}-\d{2}-\d{2})", rep.get("ledger", ""))
    when = f", {m.group(1)}" if m else ""
    if "peek" in s:
        return "header via peek"
    if "descriptor" in s:
        return f"descriptor file{when}"
    return f"downloaded{when}"


# ---- the report -------------------------------------------------------------------

def _md_table(head, rows):
    out = ["| " + " | ".join(head) + " |",
           "|" + "|".join("---:" if h in ("units", "size", "inspected", "region-only units",
                                           "shared units", "inspected (region-only / shared)",
                                           "known issues") else "---"
                          for h in head) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def _block_table(block, inspected):
    rows = []
    for (b, m, p), rs in browse.breakdown(block):
        fmts = ", ".join(sorted({r["format"] or "?" for r in rs}))
        n_in = sum(1 for r in rs if (r["bucket"], r["unit"]) in inspected)
        rows.append((browse._short_bucket(b), m or "-", p, len(rs), browse._size(rs), fmts,
                     f"{n_in} of {len(rs)}" if n_in else "—"))
    return _md_table(["bucket", "data type", "product", "units", "size", "formats",
                      "inspected"], rows)


def build_region(code, rows, hidden, inspections, issues, dates, vocab, report_dir):
    reg = vocab.regions[code]
    only = [r for r in rows if r["region_list"] == [code]]
    shared = [r for r in rows if code in r["region_list"] and len(r["region_list"]) > 1]
    mine = only + shared
    hidden_here = [r for r in hidden if code in r["region_list"]]
    inspected = attach(mine, inspections)

    L = [f"# {code} — {reg.name}", "",
         f"Allen Human Brain Atlas structure: {reg.dhba_id}"
         + (f"; also called: {', '.join(reg.aliases)}" if reg.aliases else ""), "",
         "Generated by `python scripts/index.py report " + code + "` from the catalog, the "
         "per-file inspect reports and `catalog/known_issues.csv`. Nothing here was read from "
         "S3 or from a data file while writing it.", "",
         "**Sources**", ""]
    for b in sorted({r["bucket"] for r in mine}):
        L.append(f"- catalog of `{b}`: listed {dates.get(b, 'unknown')}")
    L += [f"- inspected files: {sum(len(v) for v in inspected.values())} (each linked below, "
          f"with how it was read)", ""]

    total = browse._size(only)
    L += [f"## At a glance", "",
          f"**{code} only:** {len(only):,} units, {total}. **Shared with other regions:** "
          f"{len(shared):,} units, {browse._size(shared)}."
          + (f" **Hidden (deprecated, `previous_objects/`):** {len(hidden_here):,} units, "
             f"{browse._size(hidden_here)}." if hidden_here else ""), ""]
    if only:
        L += [f"### {code} only", ""] + _block_table(only, inspected) + [""]
    if shared:
        combos = collections.Counter(browse.region_label(r["region_list"]) for r in shared)
        L += ["### Shared with other regions", "",
              "These units cover several regions at once; their sizes are not part of the "
              f"{code}-only total. Combinations: "
              + "; ".join(f"{c} ({n})" for c, n in sorted(combos.items())), ""]
        L += _block_table(shared, inspected) + [""]

    # what was checked inside
    L += ["## What was checked inside the files", ""]
    if inspected:
        rows_out = []
        for (b, u), reps in sorted(inspected.items()):
            row = next(r for r in mine if (r["bucket"], r["unit"]) == (b, u))
            for rep in reps:
                md = rep["_json"].with_suffix(".md")
                try:
                    link = md.relative_to(report_dir.parent).as_posix()
                    link = f"[report](../{link})"
                except ValueError:
                    link = md.name
                name = rep["key"][len(u):] if rep["key"] != u else u.rsplit("/", 1)[-1]
                name = name or u.rstrip("/").rsplit("/", 1)[-1]
                rows_out.append((row["product"], f"`{name}`", how_read(rep), summarise(rep),
                                 link))
        L += _md_table(["product", "file", "how read", "structure", "details"], rows_out)
        same = collections.Counter(r["product"] for r in mine)
        checked = collections.Counter(next(r for r in mine if (r["bucket"], r["unit"]) == k)
                                      ["product"] for k in inspected)
        L += ["", f"Coverage ({code}-only and shared units together): "
              + "; ".join(f"{p}: {checked[p]} of {same[p]:,}" for p in sorted(checked))
              + ". An example stands for its kind only as far as the files were made the same "
              "way; it is not a census.", ""]
    else:
        L += ["No file of this region has been inspected yet.", ""]
    missing = sorted({r["product"] for r in mine} -
                     {next(r for r in mine if (r["bucket"], r["unit"]) == k)["product"]
                      for k in inspected})
    if missing:
        L += ["**Not yet inspected** (no example of these kinds has been looked inside): "
              + ", ".join(missing), ""]

    # issues
    known = issues_mod.matching(issues, mine)
    L += ["## Known upstream issues", ""]
    if known:
        for i, n in known:
            L += [f"- **{i['severity'].upper()} {i['id']}** — {i['summary']} "
                  f"({n:,} unit(s)). {i['detail']} *Evidence: {i['evidence']} "
                  f"(found {i['found']}).*"]
    else:
        L.append("None recorded for this region's units.")
    L += [""]

    # how to get it
    # Example download: a file someone already looked inside (its structure is
    # known), else the region's main processed dataset, else the smallest unit.
    seen = [r for r in mine if (r["bucket"], r["unit"]) in inspected and r["kind"] == "file"]
    main = [r for r in only if r["product"] in ("dataset", "cell-by-gene table")
            and r["kind"] == "file"]
    pick = sorted(seen, key=lambda r: r["size"]) or sorted(main, key=lambda r: r["size"]) \
        or sorted(only or shared, key=lambda r: r["size"])
    L += ["## How to get it", "", "```bash", "conda activate seaad",
          f"python scripts/index.py show {code}",
          f"python scripts/index.py show {code} --files"]
    if pick:
        L.append(f"python scripts/index.py get --match \"{pick[0]['unit']}\" --dry-run")
    L += ["```", "",
          "`get` always previews first and asks before downloading; `--dry-run` only "
          "previews. Folder units (Zarr stores, Xenium bundles) can be described without "
          "downloading them: `get --descriptors --match <folder>`. Image headers can be read "
          "without downloading: `peek --match <file>`.", ""]

    L += ["## How to read this report", "",
          "- *units* are files, or folders treated as one item (a Zarr store, a Xenium "
          "bundle, a tile pyramid). Sizes come from the catalog listing.",
          "- *inspected* counts units with at least one file whose structure was checked. "
          "'—' means none yet.",
          "- *structure* summarises the linked per-file report: shapes, number types and, for "
          "expression matrices, whether values look like raw counts. Category levels and "
          "identifiers are not repeated here; the per-file report shows what is safe to show.",
          "- Known issues are listed when their pattern matches one of this region's units; "
          "each says how it was established."]
    only_keys = {(r["bucket"], r["unit"]) for r in only}
    insp_only = sum(len(v) for k, v in inspected.items() if k in only_keys)
    return "\n".join(L) + "\n", {"code": code, "only": len(only), "shared": len(shared),
                                 "size_only": sum(r["size"] for r in only),
                                 "size_shared": sum(r["size"] for r in shared),
                                 "inspected": sum(len(v) for v in inspected.values()),
                                 "inspected_only": insp_only,
                                 "inspected_shared": sum(len(v) for v in inspected.values())
                                 - insp_only,
                                 "issues": len(known)}


def external_section(sources_path=None, status_path=None):
    """Sources outside the buckets, from catalog/external_sources.csv."""
    import external
    rows, status = external.load(sources_path), external.load_status(status_path)
    if not rows:
        return []
    out = ["", "## Other sources (not in the S3 buckets)", "",
           "Open files download with `python scripts/index.py get --external <id>`; "
           "controlled-access studies are linked only and never downloaded by this tool.", ""]
    table = []
    for r in rows:
        st = status.get(r["id"], {})
        size = index.human(int(st["size"])) if st.get("size") else "—"
        where = (f"[request access]({r['url']})" if r["access"] == "controlled"
                 else "`get --external " + r["id"] + "`")
        table.append((f"`{r['id']}`", r["title"], r["access"], r["modality"],
                      r["regions"] or "cohort-wide", size, where))
    out += _md_table(["id", "what", "access", "data type", "regions", "size", "how"], table)
    checked = sorted({s["checked"] for s in status.values()})
    out += ["", f"Sizes from a header-only check on {checked[-1]}." if checked else
            "Sizes not yet checked (`python scripts/index.py sources --check`)."]
    return out


def run(region=None, all_regions=False, catalog_dir=None, inspect_dir=None, report_dir=None,
        issues_path=None, vocab=None, out=print):
    catalog_dir = catalog_dir or index.CATALOG_DIR
    report_dir = pathlib.Path(report_dir or REPORT_DIR)
    vocab = vocab or regions_mod.default()
    rows, hidden = browse.filter_units(browse.load_units(catalog_dir), None, None, None, False)
    inspections = load_inspections(inspect_dir or INSPECT_DIR)
    issues = issues_mod.load(issues_path or pathlib.Path(catalog_dir) / "known_issues.csv")
    dates = plan.catalog_dates(rows, catalog_dir)
    if all_regions:
        codes = sorted(c for c in vocab.regions if any(c in r["region_list"] for r in rows))
    elif region:
        try:
            codes = [vocab.resolve(region)]
        except regions_mod.RegionError as exc:
            raise SystemExit(f"\n{exc}")
    else:
        raise SystemExit("\nGive a region (e.g. LEC) or --all.")
    report_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for code in codes:
        text, summary = build_region(code, rows, hidden, inspections, issues, dates, vocab,
                                     report_dir)
        path = report_dir / f"{code}.md"
        path.write_text(text)
        summaries.append(summary)
        out(f"  wrote {path}  ({summary['only']:,} region-only units, "
            f"{summary['inspected']} inspected file(s), {summary['issues']} known issue(s))")
    if all_regions:
        lines = [(f"[{x['code']}]({x['code']}.md)", vocab.regions[x["code"]].name,
                  f"{x['only']:,}", index.human(x["size_only"]),
                  f"{x['shared']:,}" if x["shared"] else "—",
                  f"{x['inspected_only']} / {x['inspected_shared']}", x["issues"])
                 for x in summaries]
        # Totals that add up to the catalog: every current unit is region-only
        # (exactly one region), shared (several) or has no single region.
        single = [r for r in rows if len(r["region_list"]) == 1]
        multi = [r for r in rows if len(r["region_list"]) > 1]
        none = [r for r in rows if not r["region_list"]]
        attached = attach(rows, inspections)

        def n_insp(block):
            keys = {(r["bucket"], r["unit"]) for r in block}
            return sum(len(v) for k, v in attached.items() if k in keys)
        lines += [("**all regions, region-only**", "", f"{len(single):,}",
                   index.human(sum(r["size"] for r in single)), "", str(n_insp(single)), ""),
                  ("**shared by several regions**", "each counted once", f"{len(multi):,}",
                   index.human(sum(r["size"] for r in multi)), "", str(n_insp(multi)), ""),
                  ("**no single region**", "project-wide, documentation, unassigned",
                   f"{len(none):,}", index.human(sum(r["size"] for r in none)), "",
                   str(n_insp(none)), ""),
                  ("**total current data**", "matches the catalog, deprecated excluded",
                   f"{len(rows):,}", index.human(sum(r["size"] for r in rows)), "",
                   str(n_insp(rows)), "")]
        idx = ["# SEA-AD data-structure reports", "",
               "One page per region, built from the catalog, the per-file inspect reports and "
               "the known-issues file (`python scripts/index.py report --all`).", ""]
        idx += _md_table(["region", "name", "region-only units", "size", "shared units",
                          "inspected (region-only / shared)", "known issues"], lines)
        idx += ["",
                "- *region-only* units belong to exactly one region; their sizes add up to the "
                "'all regions, region-only' row.",
                "- *shared units* cover several regions (the multi-region release, region-pair "
                "neuropathology folders such as HIP+MEC or V1C+ESOC). Each region's page lists "
                "them, so the per-region column counts a shared unit under every region it "
                "covers; the 'shared by several regions' row counts each once.",
                "- *inspected* = files whose structure was checked. A file in a shared unit is "
                "shown after the slash for every region it covers, so the per-region numbers "
                "do not add up; the bold rows count each file once.",
                "- *known issues* = issues whose pattern matches any of the region's units, "
                "region-only or shared.",
                f"- Deprecated units (`previous_objects/`) are left out everywhere: "
                f"{len(hidden):,} units, {index.human(sum(r['size'] for r in hidden))}."]
        idx += external_section()
        (report_dir / "README.md").write_text("\n".join(idx) + "\n")
        out(f"  wrote {report_dir / 'README.md'}")
    return summaries
