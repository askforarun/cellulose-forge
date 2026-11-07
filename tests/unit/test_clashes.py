"""Tests for periodic-image safety validation."""

import numpy as np
import pytest

from cellulose_forge.clashes import (
    closest_longitudinal_periodic_image_contact,
    validate_finite_periodic_image_contact,
    validate_periodic_glycosidic_bonds,
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
