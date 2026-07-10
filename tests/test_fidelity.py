import numpy as np

from src.fidelity import (
    CZ_TARGET,
    LOCAL_Z_PRODUCT,
    apply_local_z_correction,
    entangling_phase_from_diagonal,
    pedersen_fidelity,
    pedersen_fidelity_kraus,
)


def test_cz_vs_cz_is_perfect():
    assert pedersen_fidelity(CZ_TARGET, CZ_TARGET) == 1.0


def test_cz_vs_identity_known_value():
    # M = CZ, Tr(M†M)=4 and |Tr(M)|^2=4, so F=(4+4)/(4*5)=0.4.
    assert np.isclose(pedersen_fidelity(np.eye(4), CZ_TARGET), 0.4)


def test_local_z_corrects_raw_blockade_phase_gate():
    raw_gate = CZ_TARGET @ LOCAL_Z_PRODUCT  # diag(1, -1, -1, -1)
    assert np.isclose(pedersen_fidelity(raw_gate, CZ_TARGET), 0.2)
    corrected = apply_local_z_correction(raw_gate)
    np.testing.assert_allclose(corrected, CZ_TARGET, atol=1e-12)
    assert pedersen_fidelity(raw_gate, CZ_TARGET, correct_local_z=True) == 1.0


def test_kraus_formula_matches_unitary_case():
    assert pedersen_fidelity_kraus([CZ_TARGET], CZ_TARGET) == 1.0


def test_entangling_phase_for_raw_blockade_gate_is_pi():
    raw_gate = np.diag([1.0, -1.0, -1.0, -1.0])
    assert np.isclose(abs(entangling_phase_from_diagonal(raw_gate)), np.pi)


def test_random_unitary_fidelity_is_bounded():
    rng = np.random.default_rng(1234)
    z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    q, r = np.linalg.qr(z)
    phases = np.diag(r) / np.abs(np.diag(r))
    unitary = q @ np.diag(phases.conj())
    fidelity = pedersen_fidelity(unitary, CZ_TARGET)
    assert 0.0 <= fidelity <= 1.0
