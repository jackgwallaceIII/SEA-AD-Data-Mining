"""
issues.py — known problems in the upstream data, recorded once, shown wherever they matter.

Not a step. Reads catalog/known_issues.csv (written by hand, reviewed,
committed). index.py show / search / get / inspect print the issues that
apply to the units on screen; `index.py issues` lists them all.

    python scripts/index.py issues
    python scripts/index.py issues --all        # include resolved ones

Decisions made here, and why:

1. ISSUES ARE DATA, NOT CODE.
   Each problem found while inspecting the buckets (a broken link file,
   an undocumented release, a naming change between releases) is one row
   of a CSV. Adding one needs no code change, the file is easy to review
   in a diff, and it lives next to the catalog it annotates.

2. AN ISSUE SAYS WHERE IT APPLIES WITH A PATH PATTERN.
   applies_to is a glob (or several, separated by ';') matched against
   unit paths, ignoring case, the same way `index.py search` matches ('*'
   also crosses '/'). Patterns cover only what was actually checked: an
   issue seen in LEC is not claimed for every region. Whatever
   command puts a matching unit on screen also prints the issue, so the
   warning reaches the person at the moment they are about to use the file.

3. EVERY ISSUE CARRIES ITS EVIDENCE AND DATE.
   'evidence' says how it was established (which files were compared, what
   was listed) and 'found' says when. An issue may be fixed upstream
   later; the date tells the reader how old the claim is, and status
   'resolved' hides it without deleting the history.

4. THREE SEVERITIES, DEFINED BY WHAT THE READER SHOULD DO.
   problem: something is broken or wrong; do not use it as is.
   caution: usable, but it changes how results can be interpreted.
   note:    worth knowing; nothing to do.

5. THE FILE IS VALIDATED WHEN READ.
   A typo in a severity, a duplicate id, a malformed date or an unknown
   bucket stops with a message naming the row, rather than silently
   hiding an issue from everyone. 'web' (files from external_sources.csv,
   matched by host/path) is accepted alongside the three buckets.
"""

import csv
import datetime
import fnmatch
import pathlib

import pathrules

ROOT = pathlib.Path(__file__).resolve().parent.parent
ISSUES_FILE = ROOT / "catalog" / "known_issues.csv"

FIELDS = ["id", "bucket", "applies_to", "severity", "summary", "detail",
          "evidence", "found", "status"]
SEVERITIES = ("problem", "caution", "note")
STATUSES = ("open", "resolved")


class IssueError(ValueError):
    pass


def load(path=None, include_resolved=False):
    path = pathlib.Path(path or ISSUES_FILE)
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != FIELDS:
            raise IssueError(f"{path.name}: columns must be {', '.join(FIELDS)}; "
                             f"found {', '.join(reader.fieldnames or [])}")
        rows = list(reader)
    seen = set()
    for n, r in enumerate(rows, start=2):           # line 1 is the header
        where = f"{path.name} line {n} ({r['id'] or 'no id'})"
        if not r["id"] or r["id"] in seen:
            raise IssueError(f"{where}: id missing or repeated")
        seen.add(r["id"])
        if r["bucket"] not in pathrules.BUCKETS and r["bucket"] != "web":
            raise IssueError(f"{where}: unknown bucket {r['bucket']!r}")
        if r["severity"] not in SEVERITIES:
            raise IssueError(f"{where}: severity must be one of {', '.join(SEVERITIES)}")
        if r["status"] not in STATUSES:
            raise IssueError(f"{where}: status must be one of {', '.join(STATUSES)}")
        if not r["applies_to"] or not r["summary"]:
            raise IssueError(f"{where}: applies_to and summary are required")
        try:
            datetime.date.fromisoformat(r["found"])
        except ValueError:
            raise IssueError(f"{where}: found must be a date like 2026-11-22")
    return [r for r in rows if include_resolved or r["status"] == "open"]


def applies(issue, bucket, unit):
    """True if the unit matches any of the issue's ';'-separated patterns."""
    return issue["bucket"] == bucket and any(
        fnmatch.fnmatch(unit.lower(), pat.strip().lower())
        for pat in issue["applies_to"].split(";") if pat.strip())


def matching(issues, rows):
    """[(issue, number of rows it applies to)] for rows with bucket/unit keys."""
    out = []
    for i in issues:
        n = sum(1 for r in rows if applies(i, r["bucket"], r["unit"]))
        if n:
            out.append((i, n))
    order = {s: k for k, s in enumerate(SEVERITIES)}
    return sorted(out, key=lambda m: (order[m[0]["severity"]], m[0]["id"]))


def lines(matches, detail=False):
    """Printable lines, most severe first."""
    out = []
    for i, n in matches:
        out.append(f"{i['severity'].upper():8s} {i['id']}  {i['summary']}  "
                   f"[{n:,} unit(s); found {i['found']}]")
        if detail:
            out.append(f"         {i['detail']}")
            out.append(f"         evidence: {i['evidence']}")
    return out


def report(include_resolved=False, path=None, out=print):
    issues = load(path, include_resolved)
    out(f"\n{'=' * 78}\nKNOWN UPSTREAM ISSUES ({len(issues)})\n{'=' * 78}")
    for i in issues:
        status = "" if i["status"] == "open" else "  (resolved)"
        out(f"\n  {i['severity'].upper()} {i['id']}{status}  found {i['found']}")
        out(f"    {i['summary']}")
        for pat in i["applies_to"].split(";"):
            out(f"    where:    s3://{i['bucket']}/{pat.strip()}")
        out(f"    detail:   {i['detail']}")
        out(f"    evidence: {i['evidence']}")
    if not issues:
        out("  none recorded.")
    out(f"\n  severities: problem = do not use as is; caution = changes interpretation; "
        f"note = worth knowing")
    return issues
