import numpy as np
import pytest

from src.hamiltonian import build_blockade_hamiltonian, build_ideal_hamiltonian, single_atom_hamiltonian


def test_single_atom_hamiltonian_matrix_elements():
    omega = 2.0
    expected = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    np.testing.assert_allclose(single_atom_hamiltonian(omega).full(), expected)


def test_ideal_atom1_hamiltonian_projects_out_rr():
    omega = 2.0
    expected = np.zeros((3, 3), dtype=complex)
    expected[0, 2] = expected[2, 0] = 1.0
    np.testing.assert_allclose(build_ideal_hamiltonian(omega, atom_index=1).full(), expected)


def test_ideal_atom2_hamiltonian_projects_out_rr():
    omega = 2.0
    expected = np.zeros((3, 3), dtype=complex)
    expected[0, 1] = expected[1, 0] = 1.0
    np.testing.assert_allclose(build_ideal_hamiltonian(omega, atom_index=2).full(), expected)


def test_blockade_hamiltonian_matrix_elements_atom2():
    omega = 2.0
    blockade = 10.0
    expected = np.zeros((4, 4), dtype=complex)
    expected[0, 1] = expected[1, 0] = 1.0
    expected[2, 3] = expected[3, 2] = 1.0
    expected[3, 3] = blockade
    np.testing.assert_allclose(build_blockade_hamiltonian(omega, blockade, atom_index=2).full(), expected)


def test_invalid_atom_index_is_rejected():
    with pytest.raises(ValueError, match="atom_index"):
        build_ideal_hamiltonian(1.0, atom_index=0)
