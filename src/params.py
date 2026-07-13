"""Physical parameters for the Rb-87 |70S_1/2> Rydberg CZ model.

ARC is the source of atomic physics numbers.  Experimental choices such as the
Rabi frequency and atom spacing are defaults, not atomic constants.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import logging
import math

import numpy as np
from arc import PairStateInteractions, Rubidium87

logger = logging.getLogger(__name__)

# Project quantum numbers: Rb-87 |70S_1/2>.
DEFAULT_N = 70
DEFAULT_L = 0
DEFAULT_J = 0.5
DEFAULT_MJ = 0.5

# Evered-like operating point, expressed in the time units used by the solver.
DEFAULT_OMEGA_MHZ = 4.0  # Omega / 2pi, cycles per microsecond.
DEFAULT_OMEGA_RAD_PER_US = 2.0 * math.pi * DEFAULT_OMEGA_MHZ

# 3 um gives U/Omega ~= 300 for n=70; 6 um would not be strong blockade.
DEFAULT_DISTANCE_UM = 3.0

# Thermal/Doppler constants for the 780 nm + 480 nm two-photon Rb excitation.
BOLTZMANN_J_PER_K = 1.380649e-23
ATOMIC_MASS_UNIT_KG = 1.66053906660e-27
RB87_MASS_KG = 87.0 * ATOMIC_MASS_UNIT_KG
LAMBDA_LOWER_NM = 780.0
LAMBDA_UPPER_NM = 480.0
LAMBDA_LOWER_UM = LAMBDA_LOWER_NM * 1e-3
LAMBDA_UPPER_UM = LAMBDA_UPPER_NM * 1e-3
K_EFF_CYCLES_PER_UM = abs(1.0 / LAMBDA_LOWER_UM - 1.0 / LAMBDA_UPPER_UM)
K_EFF_RAD_PER_UM = 2.0 * math.pi * K_EFF_CYCLES_PER_UM
DEFAULT_TEMPERATURE_K = 10e-6

# ARC perturbative C6 calculation controls.
DEFAULT_C6_THETA = 0.0
DEFAULT_C6_PHI = 0.0
DEFAULT_C6_N_RANGE = 5
DEFAULT_C6_ENERGY_DELTA_HZ = 25e9
DEFAULT_INCLUDE_LEVELS_UP_TO = 100


@dataclass(frozen=True)
class RydbergParams:
    """Baseline parameters for the Rb-Rb Rydberg CZ simulation.

    ``c6_ghz_um6`` keeps ARC's sign.  ``blockade_shift_rad_per_us`` uses the
    magnitude because the ideal/blockade error budget cares about the detuning
    scale |U| unless a later study deliberately explores the sign.
    """

    n: int
    l: int
    j: float
    mj: float
    omega_mhz: float
    omega_rad_per_us: float
    distance_um: float
    rydberg_lifetime_s: float
    rydberg_decay_rate_s: float
    rydberg_decay_rate_per_us: float
    c6_ghz_um6: float
    c6_abs_ghz_um6: float
    blockade_shift_rad_per_us: float
    blockade_shift_mhz: float
    blockade_to_rabi: float
    intermediate_lifetime_s: float
    intermediate_decay_rate_s: float
    intermediate_decay_rate_per_us: float
    intermediate_linewidth_mhz: float

    @property
    def rydberg_lifetime_us(self) -> float:
        return self.rydberg_lifetime_s * 1e6

    def summary(self) -> str:
        return (
            f"Rb87 |{self.n}S1/2>, tau0={self.rydberg_lifetime_us:.1f} us, "
            f"C6={self.c6_ghz_um6:.1f} GHz um^6, R={self.distance_um:g} um, "
            f"U/Ω={self.blockade_to_rabi:.1f}"
        )


def thermal_velocity_rms_um_per_us(
    temperature_K: float,
    mass_kg: float = RB87_MASS_KG,
) -> float:
    """Return the 1D thermal velocity RMS in um/us.

    Numerically, ``1 m/s == 1 um/us``, so the Maxwell-Boltzmann value
    ``sqrt(k_B T / m)`` needs no scale factor before it is multiplied by a
    wavevector in rad/um to produce a detuning in rad/us.
    """

    temperature_K = float(temperature_K)
    mass_kg = float(mass_kg)
    if not np.isfinite(temperature_K) or temperature_K < 0.0:
        raise ValueError(f"temperature_K must be non-negative and finite, got {temperature_K!r}")
    if not np.isfinite(mass_kg) or mass_kg <= 0.0:
        raise ValueError(f"mass_kg must be positive and finite, got {mass_kg!r}")
    return math.sqrt(BOLTZMANN_J_PER_K * temperature_K / mass_kg)


def _validate_quantum_numbers(n: int, l: int, j: float, mj: float) -> None:
    if not isinstance(n, int) or n <= 0:
        raise ValueError(f"n must be a positive integer, got {n!r}")
    if not isinstance(l, int) or l < 0 or l >= n:
        raise ValueError(f"l must satisfy 0 <= l < n, got l={l!r}, n={n!r}")
    if j <= 0 or not np.isfinite(j):
        raise ValueError(f"j must be positive and finite, got {j!r}")
    if abs(mj) > j:
        raise ValueError(f"|mj| must be <= j, got mj={mj!r}, j={j!r}")


def _c6_from_arc(atom: Rubidium87, n: int, l: int, j: float, mj: float) -> float:
    pair = PairStateInteractions(atom, n, l, j, n, l, j, mj, mj)
    return float(
        pair.getC6perturbatively(
            DEFAULT_C6_THETA,
            DEFAULT_C6_PHI,
            DEFAULT_C6_N_RANGE,
            DEFAULT_C6_ENERGY_DELTA_HZ,
        )
    )


@lru_cache(maxsize=8)
def get_rydberg_params(
    n: int = DEFAULT_N,
    l: int = DEFAULT_L,
    j: float = DEFAULT_J,
    mj: float = DEFAULT_MJ,
    omega_mhz: float = DEFAULT_OMEGA_MHZ,
    distance_um: float = DEFAULT_DISTANCE_UM,
) -> RydbergParams:
    """Query ARC and return the default parameter bundle.

    The Rydberg lifetime is radiative-only (T=0 K).  Room-temperature blackbody
    redistribution is a different physical channel and should not be silently
    folded into the simple Lindblad decay model used later.
    """

    _validate_quantum_numbers(n, l, j, mj)
    if omega_mhz <= 0 or not np.isfinite(omega_mhz):
        raise ValueError(f"omega_mhz must be positive and finite, got {omega_mhz!r}")
    if distance_um <= 0 or not np.isfinite(distance_um):
        raise ValueError(f"distance_um must be positive and finite, got {distance_um!r}")

    atom = Rubidium87()
    lifetime_s = float(
        atom.getStateLifetime(
            n,
            l,
            j,
            temperature=0,
            includeLevelsUpTo=DEFAULT_INCLUDE_LEVELS_UP_TO,
        )
    )
    c6_ghz_um6 = _c6_from_arc(atom, n, l, j, mj)

    omega_rad_per_us = 2.0 * math.pi * omega_mhz
    c6_abs_ghz_um6 = abs(c6_ghz_um6)
    blockade_shift_mhz = c6_abs_ghz_um6 * 1000.0 / (distance_um**6)
    blockade_shift_rad_per_us = 2.0 * math.pi * blockade_shift_mhz

    rydberg_decay_rate_s = 1.0 / lifetime_s
    rydberg_decay_rate_per_us = rydberg_decay_rate_s * 1e-6

    intermediate_lifetime_s = float(
        atom.getStateLifetime(5, 1, 1.5, temperature=0, includeLevelsUpTo=DEFAULT_INCLUDE_LEVELS_UP_TO)
    )
    intermediate_decay_rate_s = 1.0 / intermediate_lifetime_s
    intermediate_decay_rate_per_us = intermediate_decay_rate_s * 1e-6
    intermediate_linewidth_mhz = intermediate_decay_rate_s / (2.0 * math.pi * 1e6)

    params = RydbergParams(
        n=n,
        l=l,
        j=j,
        mj=mj,
        omega_mhz=omega_mhz,
        omega_rad_per_us=omega_rad_per_us,
        distance_um=distance_um,
        rydberg_lifetime_s=lifetime_s,
        rydberg_decay_rate_s=rydberg_decay_rate_s,
        rydberg_decay_rate_per_us=rydberg_decay_rate_per_us,
        c6_ghz_um6=c6_ghz_um6,
        c6_abs_ghz_um6=c6_abs_ghz_um6,
        blockade_shift_rad_per_us=blockade_shift_rad_per_us,
        blockade_shift_mhz=blockade_shift_mhz,
        blockade_to_rabi=blockade_shift_rad_per_us / omega_rad_per_us,
        intermediate_lifetime_s=intermediate_lifetime_s,
        intermediate_decay_rate_s=intermediate_decay_rate_s,
        intermediate_decay_rate_per_us=intermediate_decay_rate_per_us,
        intermediate_linewidth_mhz=intermediate_linewidth_mhz,
    )
    logger.info("ARC parameter lookup: %s", params.summary())
    return params
