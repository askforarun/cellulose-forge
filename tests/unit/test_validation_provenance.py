"""Stored validation results must record where and how they were produced."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

VALIDATION = Path(__file__).resolve().parents[2] / "validation"
RESULTS = VALIDATION / "results"
MANIFEST = json.loads((VALIDATION / "release_manifest.json").read_text(encoding="utf-8"))
PROVENANCED_RESULTS = tuple(
    sorted(name for name, entry in MANIFEST["results"].items() if entry["role"] == "evidence")
)


@pytest.mark.parametrize("name", PROVENANCED_RESULTS)
def test_result_records_commit_tree_state_and_versions(name: str) -> None:
    result = json.loads((RESULTS / name).read_text(encoding="utf-8"))
    provenance = result["provenance"]
    assert len(provenance["cellulose_forge_commit"]) == 40
    assert isinstance(provenance["tracked_tree_clean"], bool)
    assert provenance["software_versions"]
    assert "not installed" not in provenance["software_versions"].values()


def test_manifest_declares_every_stored_result() -> None:
    stored = {path.name for path in RESULTS.glob("*.json")}
    assert stored == set(MANIFEST["results"])
