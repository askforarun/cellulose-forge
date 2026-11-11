"""Byte-level regression guard for coordinate-only cellulose builds.

The recorded SHA-256 hashes were produced before polymer generalization began.
Any change to cellulose coordinates, ordering, or annotations fails this test.
Regenerate deliberately with ``CELLULOSE_FORGE_UPDATE_HASHES=1`` and review the
diff of ``cellulose_hashes.json`` before committing.
"""

from __future__ import annotations

import hashlib
import json
import os
import warnings
from pathlib import Path

import pytest

from cellulose_forge.cli import main

HASH_FILE = Path(__file__).with_name("cellulose_hashes.json")
CASES: dict[str, list[str]] = {
    **{
        f"unit_cell_{key}": ["--allomorph", name, "--shape", "unit-cell"]
        for key, name in (
            ("i_alpha", "I-alpha"),
            ("i_beta", "I-beta"),
            ("ii", "II"),
            ("iii_i", "III-I"),
        )
    },
    **{
        f"bulk_{key}": ["--allomorph", name, "--shape", "bulk", "--cells", "2", "2", "2"]
        for key, name in (
            ("i_alpha", "I-alpha"),
            ("i_beta", "I-beta"),
            ("ii", "II"),
            ("iii_i", "III-I"),
        )
    },
    **{
        f"slab_{key}": [
            "--allomorph", name, "--shape", "slab", "--face", *face.split(),
            "--layers", "2", "--chains-per-layer", "4",
            "--degree-polymerization", "4", "--chain-boundary", "periodic",
            "--vacuum", "20",
        ]
        for key, name, face in (
            ("i_alpha_010", "I-alpha", "0 1 0"),
            ("i_beta_110", "I-beta", "1 1 0"),
            ("ii_110", "II", "1 1 0"),
            ("iii_i_110", "III-I", "1 1 0"),
        )
    },
    "fibril_i_beta_36": [
        "--allomorph", "I-beta", "--shape", "fibril", "--fibril-model", "36-chain",
        "--degree-polymerization", "4", "--chain-boundary", "periodic", "--vacuum", "20",
    ],
    "fibril_i_beta_18": [
        "--allomorph", "I-beta", "--shape", "fibril", "--fibril-model", "18-chain",
        "--degree-polymerization", "4", "--chain-boundary", "periodic", "--vacuum", "20",
    ],
}


def _build_hash(arguments: list[str], output: Path) -> str:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert main(["build", *arguments, "--output", str(output)]) == 0
    return hashlib.sha256(output.read_bytes()).hexdigest()


def test_cellulose_coordinate_outputs_are_unchanged(tmp_path: Path) -> None:
    actual = {
        name: _build_hash(arguments, tmp_path / f"{name}.pdb")
        for name, arguments in CASES.items()
    }
    if os.environ.get("CELLULOSE_FORGE_UPDATE_HASHES") == "1":
        HASH_FILE.write_text(json.dumps(actual, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
        pytest.skip("recorded new cellulose hashes")
    expected = json.loads(HASH_FILE.read_text(encoding="utf-8"))
    assert actual == expected
