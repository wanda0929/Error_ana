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
    unit of ``omega``; for the project default, ``omega`` is rad/us and the time
    is in microseconds.
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
