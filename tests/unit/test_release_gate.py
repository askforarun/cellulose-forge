"""Tests for the release gate's individual checks."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

_PATH = Path(__file__).resolve().parents[2] / "validation" / "release_gate.py"
_SPEC = importlib.util.spec_from_file_location("release_gate", _PATH)
assert _SPEC is not None and _SPEC.loader is not None
release_gate = importlib.util.module_from_spec(_SPEC)
sys.modules["release_gate"] = release_gate
_SPEC.loader.exec_module(release_gate)


def _statuses(gate: object, check: str) -> list[str]:
    return [finding.status for finding in gate.findings if finding.check == check]  # type: ignore[attr-defined]


def _junit(path: Path, body: str) -> Path:
    path.write_text(f'<testsuites><testsuite name="s">{body}</testsuite></testsuites>')
    return path


def test_development_mode_passes_on_the_repository() -> None:
    gate = release_gate.run(
        release_gate.argparse.Namespace(release=False, tag=None, junit=[], allow_skip=[])
    )
    errors = [finding for finding in gate.findings if finding.status == "error"]
    assert errors == []


def test_release_mode_requires_test_reports() -> None:
    gate = release_gate.Gate(release=True)
    release_gate.check_tests(gate, [], [])
    assert _statuses(gate, "tests") == ["error"]


def test_failures_fail_in_every_mode(tmp_path: Path) -> None:
    report = _junit(tmp_path / "r.xml", '<testcase classname="t" name="a"><failure/></testcase>')
    gate = release_gate.Gate(release=False)
    release_gate.check_tests(gate, [report], [])
    assert "error" in _statuses(gate, "tests")


@pytest.mark.parametrize(("release", "expected"), [(False, "warning"), (True, "error")])
def test_unexpected_skips_block_only_a_release(
    tmp_path: Path, release: bool, expected: str
) -> None:
    report = _junit(
        tmp_path / "r.xml",
        '<testcase classname="t" name="a"><skipped message="optional OpenFF environment"/>'
        "</testcase>",
    )
    gate = release_gate.Gate(release=release)
    release_gate.check_tests(gate, [report], [])
    assert _statuses(gate, "tests")[0] == expected


def test_allowed_skips_are_counted_not_reported(tmp_path: Path) -> None:
    report = _junit(
        tmp_path / "r.xml",
        '<testcase classname="t" name="a"><skipped message="needs NAMD"/></testcase>',
    )
    gate = release_gate.Gate(release=True)
    release_gate.check_tests(gate, [report], ["NAMD"])
    assert _statuses(gate, "tests") == ["ok"]


def test_mismatched_tag_is_an_error() -> None:
    gate = release_gate.Gate(release=True)
    release_gate.check_versions(gate, "v0.0.1")
    assert "error" in _statuses(gate, "tag")


def test_missing_release_tag_is_an_error(monkeypatch: pytest.MonkeyPatch) -> None:
    real_git = release_gate._git

    def fake_git(*arguments: str) -> subprocess.CompletedProcess[str]:
        if arguments[:3] == ("rev-parse", "--verify", "--quiet"):
            return subprocess.CompletedProcess(["git", *arguments], 1, "", "")
        return real_git(*arguments)

    monkeypatch.setattr(release_gate, "_git", fake_git)
    gate = release_gate.Gate(release=True)
    release_gate.check_versions(gate, "v1.3.0")
    assert any(
        finding.status == "error" and "does not exist" in finding.message
        for finding in gate.findings
        if finding.check == "tag"
    )


def test_provenance_finder_searches_nested_documents() -> None:
    commit = "a" * 40
    document = {"cases": [{"meta": {"cellulose_forge_commit": commit}}]}
    found = release_gate._provenance(document)
    assert found == {"cellulose_forge_commit": commit}
    assert release_gate._provenance({"cases": []}) is None
