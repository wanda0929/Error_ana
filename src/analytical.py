"""Closed-form error estimates for the Rydberg CZ error budget."""

from __future__ import annotations

import math

import numpy as np


def _check_positive_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite, got {value!r}")
    return value


def _check_nonnegative_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be non-negative and finite, got {value!r}")
    return value


def rydberg_decay_exposure_time(omega: float) -> float:
    """Return the effective Rydberg occupation time ``7*pi/(4*omega)``.

    ``omega`` is an angular Rabi frequency.  The returned time uses the inverse
    unit of ``omega``; for the project default, ``omega`` is in rad/us and the
    time is in microseconds.
    """

    omega = _check_positive_finite("omega", omega)
    return 7.0 * math.pi / (4.0 * omega)


def epsilon_decay_from_gamma(omega: float, gamma: float) -> float:
    """Leading Rydberg-decay infidelity ``gamma * 7*pi/(4*omega)``.

    ``gamma`` must be expressed in the inverse time unit corresponding to
    ``omega``.  In this repository that means ``gamma`` in us^-1 when ``omega``
    is in rad/us.
    """

    gamma = _check_nonnegative_finite("gamma", gamma)
    return gamma * rydberg_decay_exposure_time(omega)


def epsilon_decay(omega: float, tau: float) -> float:
    """Leading Rydberg-decay infidelity for lifetime ``tau``.

    ``tau`` must be expressed in the time unit corresponding to ``omega``.  In
    this repository that means ``tau`` in us when ``omega`` is in rad/us.
    """

    tau = _check_positive_finite("tau", tau)
    return epsilon_decay_from_gamma(omega, 1.0 / tau)


def epsilon_blockade(omega: float, blockade_shift: float) -> float:
    """Return the perturbative finite-blockade infidelity ``Ω²/(8U²)``.

    Parameters use the same angular-frequency units.  The formula only depends
    on the ratio, so ``omega`` and ``blockade_shift`` may be rad/us, rad/s, or
    any consistent unit.
    """

    omega = _check_positive_finite("omega", omega)
    blockade_shift = _check_positive_finite("blockade_shift", blockade_shift)
    return omega**2 / (8.0 * blockade_shift**2)


def epsilon_scattering(
    gamma_e: float,
    delta_p: float,
    omega1_over_omega2: float = 1.0,
) -> float:
    """Leading intermediate-state-scattering infidelity.

    ``gamma_e`` is the intermediate-state decay rate and ``delta_p`` is the
    one-photon detuning.  They must use the same angular-frequency units; in this
    repository both are normally rad/us.  ``omega1_over_omega2`` is the lower-leg
    to upper-leg single-photon Rabi-frequency ratio ``q``.

    The perturbative CZ estimate is

    ``(7*pi/8) * (gamma_e/delta_p) * (q + 1/q)/2``.
    """

    gamma_e = _check_nonnegative_finite("gamma_e", gamma_e)
    delta_p = _check_positive_finite("delta_p", delta_p)
    omega1_over_omega2 = _check_positive_finite("omega1_over_omega2", omega1_over_omega2)
    beam_imbalance = 0.5 * (omega1_over_omega2 + 1.0 / omega1_over_omega2)
    return (7.0 * math.pi / 8.0) * (gamma_e / delta_p) * beam_imbalance
