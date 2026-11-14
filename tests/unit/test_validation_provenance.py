"""Stored validation results must record where and how they were produced."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

RESULTS = Path(__file__).resolve().parents[2] / "validation" / "results"
PROVENANCED_RESULTS = (
    "openff_fibril_nagl_matrix.json",
    "openff_finite_slab_matrix.json",
    "openff_finite_slab_nagl_matrix.json",
    "openff_periodic_slab_nagl_matrix.json",
    "openff_periodic_bulk_nagl_matrix.json",
    "openff_shipped_outputs.json",
    "openff_mixed_precision_accumulation.json",
    "cellulose_ii_hydroxyl_network.json",
    "opls_cellulose_matrix.json",
    "openmm_tilted_box.json",
)


@pytest.mark.parametrize("name", PROVENANCED_RESULTS)
def test_result_records_commit_tree_state_and_versions(name: str) -> None:
    result = json.loads((RESULTS / name).read_text(encoding="utf-8"))
    provenance = result["provenance"]
    assert len(provenance["cellulose_forge_commit"]) == 40
    assert isinstance(provenance["tracked_tree_clean"], bool)
    assert provenance["software_versions"]
    assert "not installed" not in provenance["software_versions"].values()
