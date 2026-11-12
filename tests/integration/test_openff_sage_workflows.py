"""Optional integration tests for experimental production OpenFF output."""

from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path

import pytest

from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import build_bulk_crystal, build_chain_parallel_slab
from cellulose_forge.openff_sage import build_sage_system, write_sage_outputs

pytestmark = pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")


@pytest.mark.parametrize(("boundary", "padding"), (("finite", 12.0), ("periodic", 0.0)))
def test_openff_slab_outputs(tmp_path: Path, boundary: str, padding: float) -> None:
    crystal = build_chain_parallel_slab(
        load_allomorph("i_beta"),
        MillerIndex(1, 1, 0),
        layers=2,
        chains_per_layer=4,
        degree_polymerization=6,
        vacuum=30.0,
        padding_z=padding,
        chain_boundary=boundary,
    )
    build = build_sage_system(crystal, boundary)
    report = write_sage_outputs(build, tmp_path, ("gromacs", "openmm"), dispersion_correction=False)
    assert report["dispersion_correction"] is False
    assert report["gromacs_coordinates"] == "structure.gro"
    assert (tmp_path / "system.top").is_file()
    assert "DispCorr = no" in (tmp_path / "minim.mdp").read_text()


def test_openff_periodic_bulk_outputs(tmp_path: Path) -> None:
    crystal = build_bulk_crystal(load_allomorph("i_beta"), (3, 3, 2))
    build = build_sage_system(crystal, "periodic")
    report = write_sage_outputs(build, tmp_path, ("gromacs", "openmm"), dispersion_correction=True)
    assert report["dispersion_correction"] is True
    assert "DispCorr = EnerPres" in (tmp_path / "minim.mdp").read_text()
    assert (tmp_path / "system.xml").is_file()
    assert (tmp_path / "topology.cif").is_file()
