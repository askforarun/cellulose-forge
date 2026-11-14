"""One nonbonded policy and run-input format for every force field."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from cellulose_forge.run_inputs import (
    dispersion_correction,
    write_amber_run_inputs,
    write_gromacs_run_inputs,
)


def _setting(text: str, key: str) -> str:
    match = re.search(rf"^{key}\s*=\s*(\S+)", text, flags=re.MULTILINE)
    assert match is not None, key
    return match.group(1)


@pytest.mark.parametrize("force_field", ["glycam", "openff"])
def test_fitted_force_fields_use_dispersion_correction_only_without_vacuum(
    force_field: str,
) -> None:
    for shape in ("bulk", "crystal", "unit-cell"):
        assert dispersion_correction(force_field, shape) is True
    for shape in ("slab", "fibril"):
        assert dispersion_correction(force_field, shape) is False


@pytest.mark.parametrize("shape", ["bulk", "crystal", "unit-cell", "slab", "fibril"])
def test_charmm36_never_uses_dispersion_correction(shape: str) -> None:
    assert dispersion_correction("charmm36", shape) is False


def test_policy_rejects_unknown_inputs() -> None:
    with pytest.raises(ValueError, match="force field"):
        dispersion_correction("opls", "bulk")
    with pytest.raises(ValueError, match="shape"):
        dispersion_correction("glycam", "sphere")


@pytest.mark.parametrize(
    ("force_field", "cutoff", "modifier"),
    [("glycam", "1.0", None), ("charmm36", "1.0", None), ("openff", "0.9", "Potential-switch")],
)
@pytest.mark.parametrize("periodic", [False, True])
def test_gromacs_inputs_share_one_layout(
    tmp_path: Path, force_field: str, cutoff: str, modifier: str | None, periodic: bool
) -> None:
    dispersion = dispersion_correction(force_field, "bulk")
    files = write_gromacs_run_inputs(
        tmp_path, force_field=force_field, periodic_chains=periodic, dispersion=dispersion
    )
    assert files == {"gromacs_minimization_input": "minim.mdp", "gromacs_nvt_input": "nvt.mdp"}
    for name in ("minim.mdp", "nvt.mdp"):
        text = (tmp_path / name).read_text()
        assert _setting(text, "rvdw") == cutoff
        assert _setting(text, "rcoulomb") == cutoff
        assert _setting(text, "DispCorr") == ("EnerPres" if dispersion else "no")
        assert _setting(text, "periodic-molecules") == ("yes" if periodic else "no")
        if modifier is None:
            assert re.search(r"^vdw-modifier", text, flags=re.MULTILINE) is None
        else:
            assert _setting(text, "vdw-modifier") == modifier
    nvt = (tmp_path / "nvt.mdp").read_text()
    assert _setting(nvt, "nsteps") == "1000"
    assert _setting(nvt, "dt") == "0.001"


def test_charmm36_inputs_refuse_dispersion_and_name_force_switch(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="CHARMM36"):
        write_gromacs_run_inputs(
            tmp_path, force_field="charmm36", periodic_chains=True, dispersion=True
        )
    write_gromacs_run_inputs(
        tmp_path, force_field="charmm36", periodic_chains=True, dispersion=False
    )
    assert "force-switch" in (tmp_path / "minim.mdp").read_text()


def test_restraint_define_is_offered_but_not_enabled(tmp_path: Path) -> None:
    write_gromacs_run_inputs(
        tmp_path,
        force_field="glycam",
        periodic_chains=False,
        dispersion=False,
        restraint_define="POSRES_CELLULOSE",
    )
    text = (tmp_path / "nvt.mdp").read_text()
    assert "; define                  = -DPOSRES_CELLULOSE" in text
    assert re.search(r"^define", text, flags=re.MULTILINE) is None


def test_amber_inputs_disable_dispersion_in_the_ewald_namelist(
    tmp_path: Path,
) -> None:
    files = write_amber_run_inputs(tmp_path, dispersion=False)
    assert files == {"amber_minimization_input": "min.in", "amber_nvt_input": "nvt.in"}
    for name in ("min.in", "nvt.in"):
        text = (tmp_path / name).read_text()
        assert "&ewald\n  vdwmeth=0,\n/" in text
        assert text.index("vdwmeth=0") > text.index("&ewald")
        assert "cut=10.0" in text
    assert "nstlim=1000, dt=0.001" in (tmp_path / "nvt.in").read_text()


def test_amber_inputs_reject_dispersion_for_finite_chain_outputs(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="cannot enable dispersion correction"):
        write_amber_run_inputs(tmp_path, dispersion=True)
