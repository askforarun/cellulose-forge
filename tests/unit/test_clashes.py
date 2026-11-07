"""Tests for periodic-image safety validation."""

import numpy as np
import pytest

from cellulose_forge.clashes import (
    closest_longitudinal_periodic_image_contact,
    validate_finite_periodic_image_contact,
    validate_periodic_glycosidic_bonds,
    validate_periodic_gromacs_topology_bonds,
)
from cellulose_forge.glycam_mapping import GlycamMappedSystem


def _system(distance: float) -> GlycamMappedSystem:
    return GlycamMappedSystem(
        atom_names=("O4", "O1"),
        residue_names=("0GB", "4GB"),
        residue_numbers=(1, 2),
        chain_numbers=(1, 1),
        cartesian_coordinates=np.asarray([[0.0, 0.0, 0.0], [0.0, 10.0 - distance, 0.0]]),
        basis_vectors=np.diag([20.0, 10.0, 20.0]),
    )


def test_reports_closest_longitudinal_periodic_image_contact() -> None:
    contact = closest_longitudinal_periodic_image_contact(_system(4.0))
    assert contact.distance_angstrom == pytest.approx(4.0)
    assert {contact.atom_name, contact.image_atom_name} == {"O1", "O4"}


def test_rejects_heavy_atom_periodic_image_clash() -> None:
    with pytest.raises(ValueError, match=r"periodic-image clash.*0\.140 angstrom"):
        validate_finite_periodic_image_contact(_system(0.14))


def test_accepts_separated_finite_chain_ends() -> None:
    contact = validate_finite_periodic_image_contact(_system(4.0))
    assert contact.distance_angstrom == pytest.approx(4.0)


def test_rejects_long_periodic_glycosidic_bond() -> None:
    system = GlycamMappedSystem(
        atom_names=("O4", "C1", "O4", "C1"),
        residue_names=("4GB",) * 4,
        residue_numbers=(1, 1, 2, 2),
        chain_numbers=(1,) * 4,
        cartesian_coordinates=np.asarray(
            [[0.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 8.0, 0.0], [0.0, 12.0, 0.0]]
        ),
        basis_vectors=np.diag([20.0, 16.0, 20.0]),
    )
    with pytest.raises(ValueError, match=r"glycosidic bond is too long"):
        validate_periodic_glycosidic_bonds(system)


def test_validates_bonds_read_from_written_gromacs_topology(tmp_path) -> None:
    system = GlycamMappedSystem(
        atom_names=("O4", "C1", "O4", "C1"),
        residue_names=("4GB",) * 4,
        residue_numbers=(1, 1, 2, 2),
        chain_numbers=(1,) * 4,
        cartesian_coordinates=np.asarray(
            [[0.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 8.0, 0.0], [0.0, 1.0, 0.0]]
        ),
        basis_vectors=np.diag([20.0, 10.0, 20.0]),
    )
    topology = tmp_path / "cellulose.itp"
    topology.write_text(
        """[ moleculetype ]
CELLULOSE 3
[ atoms ]
1 Os 1 4GB O4 1 0 16
2 Cg 1 4GB C1 2 0 12
3 Os 2 4GB O4 3 0 16
4 Cg 2 4GB C1 4 0 12
[ bonds ]
1 4 1
3 2 1
[ system ]
Cellulose
[ molecules ]
CELLULOSE 1
""",
        encoding="utf-8",
    )
    result = validate_periodic_gromacs_topology_bonds(system, topology)
    assert result == {
        "molecule_type_count": 1,
        "molecule_count": 1,
        "bond_count": 2,
        "maximum_distance_angstrom": pytest.approx(1.0),
    }


def test_written_topology_validation_checks_every_molecule_copy(tmp_path) -> None:
    first_chain = np.asarray(
        [[0.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 8.0, 0.0], [0.0, 1.0, 0.0]]
    )
    second_chain = np.asarray(
        [[5.0, 0.0, 0.0], [5.0, 5.0, 0.0], [5.0, 8.0, 0.0], [5.0, 4.0, 0.0]]
    )
    system = GlycamMappedSystem(
        atom_names=("O4", "C1", "O4", "C1") * 2,
        residue_names=("4GB",) * 8,
        residue_numbers=(1, 1, 2, 2, 3, 3, 4, 4),
        chain_numbers=(1,) * 4 + (2,) * 4,
        cartesian_coordinates=np.vstack([first_chain, second_chain]),
        basis_vectors=np.diag([20.0, 10.0, 20.0]),
    )
    topology = tmp_path / "cellulose.top"
    topology.write_text(
        """[ moleculetype ]
CELLULOSE 3
[ atoms ]
1 Os 1 4GB O4 1 0 16
2 Cg 1 4GB C1 2 0 12
3 Os 2 4GB O4 3 0 16
4 Cg 2 4GB C1 4 0 12
[ bonds ]
1 4 1
3 2 1
[ system ]
Cellulose
[ molecules ]
CELLULOSE 2
""",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="written GROMACS glycosidic bond is too long"):
        validate_periodic_gromacs_topology_bonds(system, topology)
