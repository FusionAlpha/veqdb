"""Standalone reconstruction roots for the VEQDB compact codec."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class FitResult:
    """Raw reconstruction roots restored from a compact record.

    The arrays are stored in the representation used by the compact schema;
    no solver or full equilibrium object is required to decode or re-encode them.
    """

    Nr: int
    Nt: int
    K_max: int | None
    R0: float
    Z0: float
    a: float
    span: float
    P0: float
    B0: float
    kappa_lcfs: float
    c_lcfs: np.ndarray
    s_lcfs: np.ndarray
    h_coeffs: np.ndarray
    v_coeffs: np.ndarray
    kappa_coeffs: np.ndarray
    c_coeffs: np.ndarray
    s_coeffs: np.ndarray
    psi_r: np.ndarray
    psi_rr: np.ndarray
    P_psi: np.ndarray
    FF_psi: np.ndarray
    psin: np.ndarray
    detail: str


def _harmonic_order(
    c_lcfs: np.ndarray, s_lcfs: np.ndarray, c_coeffs: np.ndarray, s_coeffs: np.ndarray
) -> int:
    return max(0, c_lcfs.size - 1, s_lcfs.size, c_coeffs.shape[0] - 1, s_coeffs.shape[0])


def _pad(values: np.ndarray, count: int) -> np.ndarray:
    if values.shape[0] == count:
        return values
    padded = np.zeros((count, *values.shape[1:]), dtype=np.float64)
    padded[: values.shape[0]] = values
    return padded
