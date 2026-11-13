"""Force-field-consistent nonbonded policy and starting run inputs.

Every force-field build applies the same rules, whichever force field or MD
engine is selected:

* Long-range dispersion correction is used only for homogeneous periodic
  systems (bulk crystals and unit cells) and only for force fields fitted with
  it (GLYCAM06 and OpenFF Sage).  CHARMM36 is parameterized without it, and
  slabs and fibrils contain vacuum, which the isotropic correction does not
  describe.
* Each engine receives a short minimization and a 1 ps NVT smoke-test input in
  its native format.  These are starting protocols, not production settings.
"""

from __future__ import annotations

from pathlib import Path

_HOMOGENEOUS_SHAPES = frozenset({"bulk", "crystal", "unit-cell"})
_VACUUM_SHAPES = frozenset({"slab", "fibril"})
_FORCE_FIELDS = frozenset({"glycam", "charmm36", "openff"})

# Longest nonbonded interaction range of any generated engine input, in
# angstrom (CHARMM36 reaches 12 A in the LAMMPS 10--12 A force switch).  Used
# only to warn when slab or fibril vacuum is thinner than this range.
INTERACTION_RANGE_ANGSTROM = {"glycam": 10.0, "charmm36": 12.0, "openff": 9.0}
# Real-space cutoff of the generated GROMACS and OpenMM inputs, in angstrom.
# GROMACS and OpenMM evaluate pairs by the minimum-image convention, so every
# perpendicular box height must exceed twice this cutoff.  LAMMPS builds ghost
# atoms for every periodic image within its cutoff and has no such limit, but
# its files are produced through an OpenMM system, so the same box rule
# applies to every build.
GROMACS_CUTOFF_ANGSTROM = {"glycam": 10.0, "charmm36": 10.0, "openff": 9.0}


def _require_force_field(force_field: str) -> str:
    name = force_field.lower()
    if name not in _FORCE_FIELDS:
        raise ValueError(f"no run-input policy is defined for force field {force_field!r}")
    return name


def dispersion_correction(force_field: str, shape: str) -> bool:
    """Return whether long-range dispersion correction applies to a build."""
    name = _require_force_field(force_field)
    if shape not in _HOMOGENEOUS_SHAPES | _VACUUM_SHAPES:
        raise ValueError(f"dispersion policy is undefined for shape {shape!r}")
    if name == "charmm36":
        return False
    return shape in _HOMOGENEOUS_SHAPES


def is_vacuum_shape(shape: str) -> bool:
    """Return whether a shape is built with transverse vacuum."""
    return shape in _VACUUM_SHAPES


def _nonbonded_block(force_field: str, *, dispersion: bool) -> str:
    dispcorr = "EnerPres" if dispersion else "no"
    if force_field == "openff":
        return f"""cutoff-scheme           = Verlet
nstlist                 = 20
coulombtype             = PME
rcoulomb                = 0.9
fourierspacing          = 0.12
pme-order               = 4
vdwtype                 = Cut-off
vdw-modifier            = Potential-switch
rvdw-switch             = 0.8
rvdw                    = 0.9
DispCorr                = {dispcorr}
"""
    block = f"""cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
fourierspacing          = 0.12
pme-order               = 4
vdwtype                 = Cut-off
rvdw                    = 1.0
DispCorr                = {dispcorr}
"""
    if force_field == "charmm36":
        block = (
            "; CHARMM36 CROSS-ENGINE SMOKE-TEST SETTINGS, NOT THE CHARMM36 PRODUCTION\n"
            "; PROTOCOL.  The plain 1.0 nm cut-off matches the validated GROMACS/OpenMM\n"
            "; comparison (OpenMM cannot reproduce the CHARMM force switch).  For\n"
            "; production use vdw-modifier = force-switch, rvdw-switch = 1.0 and\n"
            "; rvdw = rcoulomb = 1.2, which requires every box height above 2.4 nm.\n"
            + block
        )
    return block


def write_gromacs_run_inputs(
    directory: Path,
    *,
    force_field: str,
    periodic_chains: bool,
    dispersion: bool,
    restraint_define: str | None = None,
) -> dict[str, str]:
    """Write ``minim.mdp`` and ``nvt.mdp`` smoke-test inputs."""
    name = _require_force_field(force_field)
    if name == "charmm36" and dispersion:
        raise ValueError("CHARMM36 run inputs must not enable dispersion correction")
    nonbonded = _nonbonded_block(name, dispersion=dispersion)
    periodic = "yes" if periodic_chains else "no"
    boundary = "covalently periodic" if periodic_chains else "chemically capped finite"
    restraint = (
        f"; Generated position restraints: uncomment to enable them.\n"
        f"; define                  = -D{restraint_define}\n"
        if restraint_define
        else ""
    )
    minim = f"""; Steepest-descent minimization for {boundary} chains.
integrator              = steep
nsteps                  = 10000
emtol                   = 1000.0
emstep                  = 0.01
{restraint}{nonbonded}pbc                     = xyz
periodic-molecules      = {periodic}
constraints             = none
"""
    nvt = f"""; One-ps NVT smoke test for {boundary} chains; not a production protocol.
integrator              = md
dt                      = 0.001
nsteps                  = 1000
continuation            = no
gen-vel                 = yes
gen-temp                = 300
gen-seed                = -1
{restraint}{nonbonded}pbc                     = xyz
periodic-molecules      = {periodic}
constraints             = h-bonds
constraint-algorithm    = lincs
lincs-order             = 4
lincs-iter              = 1
tcoupl                  = v-rescale
tc-grps                 = System
tau-t                   = 0.5
ref-t                   = 300
pcoupl                  = no
nstxout-compressed      = 100
nstenergy               = 100
nstlog                  = 100
nstcalcenergy           = 100
comm-mode               = Linear
nstcomm                 = 100
"""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "minim.mdp").write_text(minim, encoding="utf-8")
    (directory / "nvt.mdp").write_text(nvt, encoding="utf-8")
    return {"gromacs_minimization_input": "minim.mdp", "gromacs_nvt_input": "nvt.mdp"}


def write_amber_run_inputs(directory: Path, *, dispersion: bool) -> dict[str, str]:
    """Write GLYCAM ``min.in`` and ``nvt.in`` smoke-test inputs for Amber.

    ``vdwmeth`` selects Amber's continuum long-range dispersion correction, so
    it follows the same geometry policy as the other engines.
    """
    vdwmeth = 1 if dispersion else 0
    minim = f"""Minimization (GLYCAM06j-1, finite chains)
&cntrl
  imin=1, maxcyc=1000, ncyc=200,
  ntb=1, cut=10.0, vdwmeth={vdwmeth},
/
"""
    nvt = f"""One-ps NVT smoke test (GLYCAM06j-1, finite chains); not a production protocol
&cntrl
  imin=0, irest=0, ntx=1,
  nstlim=1000, dt=0.001,
  ntc=2, ntf=2,
  ntb=1, cut=10.0, vdwmeth={vdwmeth},
  ntt=3, gamma_ln=1.0, temp0=300.0, tempi=300.0, ig=-1,
  ntpr=100, ntwx=100, ntwe=100,
/
"""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "min.in").write_text(minim, encoding="utf-8")
    (directory / "nvt.in").write_text(nvt, encoding="utf-8")
    return {"amber_minimization_input": "min.in", "amber_nvt_input": "nvt.in"}
