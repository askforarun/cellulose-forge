"""Check that the repository is in a releasable state.

Run before tagging a release::

    uv run python validation/release_gate.py --release --tag v1.3.0 \
        --junit core-tests.xml --junit openff-tests.xml

Without ``--release`` the gate runs in development mode: it reports every
finding but fails only on structural problems (undeclared or missing result
files, inconsistent versions, failed tests). ``--release`` additionally fails on
a dirty tree, stale or unprovenanced evidence, unreleased changelog entries, a
tag that does not match, and skipped tests that were not explicitly allowed.

Each finding is ``error`` (fails this mode), ``warning`` (reported only), or
``ok``. ``--json`` writes the findings for CI or an agent to read.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tomllib
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "results"
MANIFEST = ROOT / "validation" / "release_manifest.json"
ROLES = frozenset({"evidence", "unprovenanced-evidence", "reference"})
# Paths whose changes invalidate every stored result: the package (code and
# crystallographic data) and the shared provenance helper.
SOURCE_PATHS = ("src/cellulose_forge", "validation/_provenance.py")
CITING_PATHS = ("README.md", "CHANGELOG.md", "docs", "validation")
_RESULT_REFERENCE = re.compile(r"results/([A-Za-z0-9_]+\.json)")
_COMMIT = re.compile(r"[0-9a-f]{40}")


@dataclass(frozen=True)
class Finding:
    check: str
    status: str
    message: str


class Gate:
    def __init__(self, release: bool) -> None:
        self.release = release
        self.findings: list[Finding] = []

    def ok(self, check: str, message: str) -> None:
        self.findings.append(Finding(check, "ok", message))

    def warn(self, check: str, message: str) -> None:
        self.findings.append(Finding(check, "warning", message))

    def error(self, check: str, message: str) -> None:
        self.findings.append(Finding(check, "error", message))

    def release_error(self, check: str, message: str) -> None:
        """An error at release time and a warning during development."""
        (self.error if self.release else self.warn)(check, message)

    @property
    def passed(self) -> bool:
        return all(finding.status != "error" for finding in self.findings)


def _git(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", *arguments], cwd=ROOT, text=True, capture_output=True, check=False
    )


def check_tree(gate: Gate) -> None:
    status = _git("status", "--porcelain", "--untracked-files=no").stdout.strip()
    if status:
        changed = len(status.splitlines())
        gate.release_error("tree", f"{changed} tracked file(s) have uncommitted changes")
    else:
        gate.ok("tree", "tracked tree is clean")


def _manifest(gate: Gate) -> dict[str, dict[str, str | None]]:
    try:
        results = json.loads(MANIFEST.read_text(encoding="utf-8"))["results"]
    except (OSError, ValueError, KeyError) as error:
        gate.error("manifest", f"cannot read {MANIFEST.relative_to(ROOT)}: {error}")
        return {}
    for name, entry in results.items():
        if entry.get("role") not in ROLES:
            gate.error("manifest", f"{name}: unknown role {entry.get('role')!r}")
        validator = entry.get("validator")
        if validator is not None and not (ROOT / validator).is_file():
            gate.error("manifest", f"{name}: validator {validator} does not exist")
    return dict(results)


def check_declared_results(gate: Gate, manifest: dict[str, dict[str, str | None]]) -> None:
    tracked = set(_git("ls-files", "validation/results").stdout.split())
    on_disk = {path.name for path in RESULTS.glob("*.json")}
    for name in sorted(on_disk - manifest.keys()):
        gate.error("manifest", f"{name} is in validation/results but not declared in the manifest")
    for name in sorted(manifest):
        if name not in on_disk:
            gate.error("results", f"{name} is declared but missing")
        elif f"validation/results/{name}" not in tracked:
            gate.error("results", f"{name} exists but is not tracked by git")
    cited: dict[str, set[str]] = {}
    listing = _git("ls-files", "--", *CITING_PATHS).stdout.split()
    for relative in listing:
        if relative.startswith("validation/results/") or not relative.endswith((".md", ".py")):
            continue
        text = (ROOT / relative).read_text(encoding="utf-8", errors="replace")
        for match in _RESULT_REFERENCE.findall(text):
            cited.setdefault(match, set()).add(relative)
    for name, sources in sorted(cited.items()):
        if name not in manifest and "{" not in name:
            gate.error(
                "citations", f"{name} is cited in {', '.join(sorted(sources))} but undeclared"
            )
    gate.ok("citations", f"{len(cited)} cited result file(s) checked against the manifest")


def _provenance(document: object) -> dict[str, object] | None:
    """Return the first mapping that records a CelluloseForge commit."""
    if isinstance(document, dict):
        if "cellulose_forge_commit" in document:
            return document
        for value in document.values():
            found = _provenance(value)
            if found is not None:
                return found
    elif isinstance(document, list):
        for value in document:
            found = _provenance(value)
            if found is not None:
                return found
    return None


def check_evidence(gate: Gate, manifest: dict[str, dict[str, str | None]]) -> None:
    fresh = 0
    if _git("rev-parse", "--is-shallow-repository").stdout.strip() == "true":
        gate.release_error(
            "freshness", "shallow clone: commit history unavailable, provenance not checked"
        )
        return
    for name, entry in sorted(manifest.items()):
        role = entry.get("role")
        path = RESULTS / name
        if role == "reference" or not path.is_file():
            continue
        if role == "unprovenanced-evidence":
            gate.release_error(
                "provenance", f"{name}: validator does not record provenance yet; regenerate it"
            )
            continue
        provenance = _provenance(json.loads(path.read_text(encoding="utf-8")))
        if provenance is None:
            gate.error("provenance", f"{name}: declared evidence has no provenance block")
            continue
        commit = str(provenance.get("cellulose_forge_commit", ""))
        if not _COMMIT.fullmatch(commit):
            gate.error("provenance", f"{name}: recorded commit {commit!r} is not a full SHA")
            continue
        if _git("merge-base", "--is-ancestor", commit, "HEAD").returncode != 0:
            gate.error("provenance", f"{name}: commit {commit[:10]} is not an ancestor of HEAD")
            continue
        if provenance.get("tracked_tree_clean") is not True:
            gate.release_error("provenance", f"{name}: generated from a tree not recorded as clean")
            continue
        validator = entry.get("validator")
        watched = [*SOURCE_PATHS, *([validator] if validator else [])]
        changed = _git("diff", "--name-only", commit, "HEAD", "--", *watched).stdout.split()
        if changed:
            gate.release_error(
                "freshness",
                f"{name}: {len(changed)} source file(s) changed since {commit[:10]} "
                f"(e.g. {changed[0]})",
            )
            continue
        fresh += 1
    gate.ok("freshness", f"{fresh} evidence file(s) are current with HEAD")


def _changelog_versions() -> tuple[list[tuple[str, str]], list[str]]:
    lines = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8").splitlines()
    released: list[tuple[str, str]] = []
    unreleased: list[str] = []
    section: str | None = None
    for line in lines:
        heading = re.match(r"^## (.+)$", line)
        if heading:
            section = heading.group(1).strip()
            match = re.match(r"(\d+\.\d+\.\d+)\s+-\s+(\d{4}-\d{2}-\d{2})$", section)
            if match:
                released.append((match.group(1), match.group(2)))
            continue
        if section == "Unreleased" and line.strip().startswith(("-", "*")):
            unreleased.append(line.strip())
    return released, unreleased


def _citation_field(text: str, field: str) -> str | None:
    match = re.search(rf"^{field}:\s*\"?([^\"\n]+)\"?\s*$", text, re.MULTILINE)
    return match.group(1).strip() if match else None


def check_versions(gate: Gate, tag: str | None) -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    version = str(project["project"]["version"])
    init = (ROOT / "src" / "cellulose_forge" / "__init__.py").read_text(encoding="utf-8")
    fallback = re.search(r'__version__ = "([^"]+)"', init)
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    sources = {
        "src/cellulose_forge/__init__.py fallback": fallback.group(1) if fallback else None,
        "CITATION.cff version": _citation_field(citation, "version"),
    }
    released, unreleased = _changelog_versions()
    sources["CHANGELOG.md latest release"] = released[0][0] if released else None
    mismatched = {where: found for where, found in sources.items() if found != version}
    for where, found in mismatched.items():
        gate.error("version", f"{where} is {found!r}, pyproject.toml is {version!r}")
    if not mismatched:
        gate.ok("version", f"version {version} is consistent across {len(sources) + 1} files")
    date = _citation_field(citation, "date-released")
    if released and date != released[0][1]:
        gate.release_error(
            "version", f"CITATION.cff date-released {date} differs from CHANGELOG {released[0][1]}"
        )
    if unreleased:
        gate.release_error(
            "changelog", f"{len(unreleased)} entry(ies) remain under '## Unreleased'"
        )
    doi = re.search(r"type:\s*doi\s*\n\s*value:\s*\"?([^\"\s]+)", citation)
    if doi is None:
        gate.warn("doi", "CITATION.cff records no DOI")
    else:
        gate.warn(
            "doi",
            f"CITATION.cff DOI is {doi.group(1)}; confirm it is the intended concept "
            "or version DOI after Zenodo archives the release",
        )
    if tag is not None:
        if tag != f"v{version}":
            gate.error("tag", f"tag {tag} does not match version v{version}")
        tagged = _git("rev-parse", "--verify", "--quiet", f"refs/tags/{tag}^{{commit}}").stdout
        head = _git("rev-parse", "HEAD").stdout.strip()
        if tagged.strip() and tagged.strip() != head:
            gate.error("tag", f"existing tag {tag} points at {tagged.strip()[:10]}, not HEAD")
        elif tag == f"v{version}":
            gate.ok("tag", f"tag {tag} matches the version and HEAD")


def check_tests(gate: Gate, reports: list[Path], allowed_skips: list[str]) -> None:
    if not reports:
        gate.release_error("tests", "no --junit report supplied; test results not checked")
        return
    patterns = [re.compile(pattern) for pattern in allowed_skips]
    totals = {"tests": 0, "failed": 0, "skipped": 0, "allowed": 0}
    for report in reports:
        try:
            root = ET.parse(report).getroot()
        except (OSError, ET.ParseError) as error:
            gate.error("tests", f"cannot read {report}: {error}")
            continue
        for case in root.iter("testcase"):
            totals["tests"] += 1
            name = f"{case.get('classname')}::{case.get('name')}"
            if case.find("failure") is not None or case.find("error") is not None:
                totals["failed"] += 1
                gate.error("tests", f"{report.name}: {name} failed")
            skipped = case.find("skipped")
            if skipped is not None:
                reason = skipped.get("message", "")
                if any(pattern.search(f"{name} {reason}") for pattern in patterns):
                    totals["allowed"] += 1
                else:
                    totals["skipped"] += 1
                    gate.release_error("tests", f"{report.name}: {name} skipped ({reason})")
    if totals["tests"] == 0:
        gate.error("tests", "the supplied reports contain no test cases")
    else:
        gate.ok(
            "tests",
            f"{totals['tests']} test(s) read: {totals['failed']} failed, "
            f"{totals['skipped']} unexpected skip(s), {totals['allowed']} allowed skip(s)",
        )


def run(arguments: argparse.Namespace) -> Gate:
    gate = Gate(release=arguments.release)
    check_tree(gate)
    manifest = _manifest(gate)
    if manifest:
        check_declared_results(gate, manifest)
        check_evidence(gate, manifest)
    check_versions(gate, arguments.tag)
    check_tests(gate, arguments.junit, arguments.allow_skip)
    return gate


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--release", action="store_true", help="apply release-time strictness")
    parser.add_argument("--tag", help="release tag to check, for example v1.3.0")
    parser.add_argument(
        "--junit", type=Path, action="append", default=[], help="pytest JUnit XML report"
    )
    parser.add_argument(
        "--allow-skip",
        action="append",
        default=[],
        help="regular expression matched against 'test-id skip-reason'",
    )
    parser.add_argument("--json", type=Path, help="write findings as JSON to this path")
    arguments = parser.parse_args(argv)
    gate = run(arguments)
    for finding in gate.findings:
        marker = {"ok": "  ok ", "warning": "warn ", "error": " FAIL"}[finding.status]
        print(f"{marker} [{finding.check}] {finding.message}")
    mode = "release" if gate.release else "development"
    print(f"\nrelease gate ({mode} mode): {'PASSED' if gate.passed else 'FAILED'}")
    if arguments.json is not None:
        arguments.json.write_text(
            json.dumps(
                {
                    "mode": mode,
                    "passed": gate.passed,
                    "findings": [asdict(finding) for finding in gate.findings],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    return 0 if gate.passed else 1


if __name__ == "__main__":
    sys.exit(main())
