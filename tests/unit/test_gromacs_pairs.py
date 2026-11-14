"""Tests for preserving explicit GROMACS pair parameters across engines."""

from __future__ import annotations

from pathlib import Path

import pytest

from cellulose_forge.gromacs_pairs import prepare_explicit_pairs


def test_explicit_pairs_are_rewritten_and_expanded(tmp_path: Path) -> None:
    source = tmp_path / "source.top"
    source.write_text(
        """[ defaults ]
1 2 yes 1 1
[ moleculetype ]
CHAIN 3
[ atoms ]
1 A 1 RES A 1 0.3 1
2 B 1 RES B 2 -0.6 16
[ pairs ]
1 2 2 0.8333333333 0.3 -0.6 0.2 0.12
[ system ]
test
[ molecules ]
CHAIN 2
""",
        encoding="utf-8",
    )
    rewritten, pairs = prepare_explicit_pairs(source, tmp_path / "parser.top")
    assert "1       2 1" in rewritten.read_text(encoding="utf-8")
    assert [(pair.first, pair.second) for pair in pairs] == [(0, 1), (2, 3)]
    assert pairs[0].charge_product_e2 == pytest.approx(-0.15)
    assert pairs[0].sigma_nm == pytest.approx(0.2)
    assert pairs[0].epsilon_kj_mol == pytest.approx(0.12)


def test_topology_without_explicit_pairs_is_not_copied(tmp_path: Path) -> None:
    source = tmp_path / "source.top"
    source.write_text("[ pairs ]\n1 2 1\n", encoding="utf-8")
    rewritten, pairs = prepare_explicit_pairs(source, tmp_path / "unused.top")
    assert rewritten == source
    assert pairs == ()
    assert not (tmp_path / "unused.top").exists()
