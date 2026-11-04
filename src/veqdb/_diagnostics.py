"""Module: fitting.diagnostics

Role:
- Evaluate the display-only diagnostic block of the compact schema from a
  ``FitResult``'s raw roots, reproducing the fusionprime ``Equilibrium``
  property chain operation-for-operation.

Public API:
- ``Diagnostics`` -- frozen scalars (``R_axis`` ... ``psi_lcfs``) plus the
  intermediate arrays (``q``, ``psin``, ``psi``, ``F``, ``P``, ``V_r``) the
  compact writer and the fidelity layer consume.
- ``compute`` -- evaluate ``Diagnostics`` for one accepted fit.

Notes:
- Sign and gauge discipline: every sign- or gauge-critical quantity is derived
  from the raw roots through the lifted kernels in ``veqdb._geometry``;
  sign-regularized or zero-edge-gauge solver profiles are never read here.
- Bitwise contract: the scalar expressions mirror the fusionprime properties
  exactly -- ``_update_matvec``-style sequential sums instead of BLAS,
  ``np.tensordot(weights, values)`` for ``integrate_radial``, ``psi[-1]``
  for ``psi_lcfs``, analytic series evaluation for the axis coordinates --
  so the same fit yields identical bytes on both implementations.
- ``psin`` here is the fusionprime property (``psi / integrate_radial(psi_r)``,
  unpinned), not the pinned array stored on ``FitResult``; ``q95``
  interpolates on this array, which is what the reference writer does.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._geometry import (
    MU0,
    _evaluate_profile_at_axis,
    _integrate_G1,
    _materialize_basic_profile,
    _uniform_theta,
    _update_betat,
    _update_chebyshev_fields,
    _update_covariant_metric,
    _update_edge_anchored_integral,
    _update_eta,
    _update_eta_gradient,
    _update_eta_hessian,
    _update_fourier_fields,
    _update_G_fields,
    _update_gs_radial_fields,
    _update_gs_surface_fields,
    _update_jacobian,
    _update_matvec,
    _update_poloidal_integral,
    _update_q,
    _update_R,
    _update_R_gradient,
    _update_R_hessian,
    _update_scaled_product,
    _update_scaled_profile,
    _update_shaping_fields,
    _update_signed_sqrt,
    _update_trigonometry,
    _update_Z,
    _update_Z_gradient,
    _update_Z_hessian,
)
from ._spectral import _lobatto_quadrature, _spectral_calculus
from ._result import FitResult, _harmonic_order, _pad

__all__ = ["Diagnostics", "compute"]


@dataclass(frozen=True, slots=True)
class Diagnostics:
    """Compact-schema diagnostics derived from one fit's raw roots."""

    R_axis: float
    Z_axis: float
    B_axis: float
    P_axis: float
    qmin: float
    q0: float
    q95: float
    Ip: float
    betat: float
    Wth: float
    psi_lcfs: float
    q: np.ndarray
    psin: np.ndarray
    psi: np.ndarray
    F: np.ndarray
    P: np.ndarray
    V_r: np.ndarray


def compute(result: FitResult) -> Diagnostics:
    """Evaluate the compact diagnostics for one accepted fit."""

    Nr, Nt, K_max = result.Nr, result.Nt, result.K_max
    R0, Z0, a = result.R0, result.Z0, result.a
    B0, P0 = result.B0, result.P0

    r, weights = _lobatto_quadrature(Nr)
    accumulator, _ = _spectral_calculus(Nr)
    theta = _uniform_theta(Nt)

    M_max = _harmonic_order(result.c_lcfs, result.s_lcfs, result.c_coeffs, result.s_coeffs)
    aligned_c_lcfs = _pad(result.c_lcfs, M_max + 1)
    aligned_s_lcfs = _pad(result.s_lcfs, M_max)
    aligned_c_coeffs = _pad(result.c_coeffs, M_max + 1)
    aligned_s_coeffs = _pad(result.s_coeffs, M_max)

    L_max = int(result.h_coeffs.size) - 1
    chebyshev = np.empty((3, L_max + 1, Nr), dtype=np.float64)
    _update_chebyshev_fields(chebyshev, r)
    T, T_r, T_rr = chebyshev[0], chebyshev[1], chebyshev[2]

    K_values = np.arange(M_max + 1, dtype=np.int64)
    if K_max is not None:
        K_values[1:] = np.minimum(K_values[1:], K_max)

    h_fields = _materialize_basic_profile(result.h_coeffs, T, T_r, T_rr, r)
    v_fields = _materialize_basic_profile(result.v_coeffs, T, T_r, T_rr, r)
    kappa_fields = _materialize_basic_profile(
        result.kappa_coeffs, T, T_r, T_rr, r, lcfs=result.kappa_lcfs
    )
    c_fields = np.empty((M_max + 1, 3, Nr), dtype=np.float64)
    s_fields = np.empty((M_max + 1, 3, Nr), dtype=np.float64)
    _update_shaping_fields(c_fields, aligned_c_coeffs, aligned_c_lcfs, T, T_r, T_rr, r, K_values, 0)
    _update_shaping_fields(s_fields, aligned_s_coeffs, aligned_s_lcfs, T, T_r, T_rr, r, K_values, 1)

    phase = np.empty((6, M_max + 1, Nt), dtype=np.float64)
    _update_fourier_fields(phase, theta)
    cos_mtheta, sin_mtheta = phase[0], phase[1]
    m_cos_mtheta, m_sin_mtheta = phase[2], phase[3]
    m2_cos_mtheta, m2_sin_mtheta = phase[4], phase[5]

    eta = np.empty((Nr, Nt), dtype=np.float64)
    _update_eta(eta, theta, c_fields, s_fields, cos_mtheta, sin_mtheta)
    eta_gradient = np.empty((2, Nr, Nt), dtype=np.float64)
    _update_eta_gradient(
        eta_gradient, c_fields, s_fields, cos_mtheta, sin_mtheta, m_cos_mtheta, m_sin_mtheta
    )
    eta_hessian = np.empty((3, Nr, Nt), dtype=np.float64)
    _update_eta_hessian(
        eta_hessian,
        c_fields,
        s_fields,
        cos_mtheta,
        sin_mtheta,
        m_cos_mtheta,
        m_sin_mtheta,
        m2_cos_mtheta,
        m2_sin_mtheta,
    )
    eta_trigonometry = np.empty((2, Nr, Nt), dtype=np.float64)
    _update_trigonometry(eta_trigonometry, eta)
    cosine, sine = eta_trigonometry[0], eta_trigonometry[1]
    theta_trigonometry = np.empty((2, Nt), dtype=np.float64)
    _update_trigonometry(theta_trigonometry, theta)
    theta_sine, theta_cosine = theta_trigonometry[1][None, :], theta_trigonometry[0][None, :]

    R = np.empty((Nr, Nt), dtype=np.float64)
    _update_R(R, R0, a, r, h_fields[0], cosine)
    if not np.all(np.isfinite(R)) or np.any(R <= 0.0):
        raise ValueError("Equilibrium produces a non-positive or non-finite major radius")
    R_gradient = np.empty((2, Nr, Nt), dtype=np.float64)
    _update_R_gradient(
        R_gradient, a, r, h_fields[1], cosine, sine, eta_gradient[0], eta_gradient[1]
    )
    R_hessian = np.empty((3, Nr, Nt), dtype=np.float64)
    _update_R_hessian(
        R_hessian,
        a,
        r,
        h_fields[2],
        cosine,
        sine,
        eta_gradient[0],
        eta_gradient[1],
        eta_hessian[0],
        eta_hessian[1],
        eta_hessian[2],
    )
    Z = np.empty((Nr, Nt), dtype=np.float64)
    _update_Z(Z, Z0, a, r, v_fields[0], kappa_fields[0], theta_sine[0])
    if not np.all(np.isfinite(Z)):
        raise ValueError("Equilibrium produces a non-finite vertical coordinate")
    Z_gradient = np.empty((2, Nr, Nt), dtype=np.float64)
    _update_Z_gradient(Z_gradient, a, r, v_fields[1], kappa_fields, theta_sine[0], theta_cosine[0])
    Z_hessian = np.empty((3, Nr, Nt), dtype=np.float64)
    _update_Z_hessian(Z_hessian, a, r, v_fields[2], kappa_fields, theta_sine[0], theta_cosine[0])

    jacobian = np.empty((Nr, Nt), dtype=np.float64)
    _update_jacobian(jacobian, R_gradient[0], R_gradient[1], Z_gradient[0], Z_gradient[1])
    if not np.all(np.isfinite(jacobian)):
        raise ValueError("Equilibrium Jacobian contains non-finite values")
    first_physical_row = 1 if abs(float(r[0])) <= 1.0e-14 else 0
    if np.any(jacobian[first_physical_row:] <= 0.0):
        raise ValueError("Equilibrium is singular or reverses orientation away from the axis")
    metric = np.empty((3, Nr, Nt), dtype=np.float64)
    _update_covariant_metric(metric, R_gradient[0], R_gradient[1], Z_gradient[0], Z_gradient[1])
    if not np.all(np.isfinite(metric)):
        raise ValueError("Equilibrium metric contains non-finite values")

    gs_surface = np.empty((4, Nr, Nt), dtype=np.float64)
    invalid = _update_gs_surface_fields(
        gs_surface,
        r,
        R,
        R_gradient[0],
        R_gradient[1],
        R_hessian[0],
        R_hessian[1],
        R_hessian[2],
        Z_gradient[0],
        Z_gradient[1],
        Z_hessian[0],
        Z_hessian[1],
        Z_hessian[2],
    )
    if invalid:
        raise ValueError(
            "Equilibrium GS metric is singular away from the magnetic axis at radial index "
            f"{invalid - 1}"
        )
    JdivR, gttdivJR, gttdivJR_r, grtdivJR_t = (
        gs_surface[0],
        gs_surface[1],
        gs_surface[2],
        gs_surface[3],
    )
    gs_radial = np.empty((3, Nr), dtype=np.float64)
    _update_gs_radial_fields(gs_radial, gttdivJR, gttdivJR_r, JdivR)
    Ln_r = gs_radial[2]

    psi_r = result.psi_r
    psi = np.empty(Nr, dtype=np.float64)
    _update_matvec(psi, accumulator, psi_r)
    # integrate_radial(psi_r): fusionprime's _axis_to_lcfs_flux.
    span = float(np.tensordot(weights, psi_r, axes=(0, 0)))
    psin = np.empty_like(psi)
    _update_scaled_profile(psin, psi, 1.0 / span)

    FF_r = np.empty(Nr, dtype=np.float64)
    _update_scaled_product(FF_r, result.FF_psi, psi_r, 1.0)
    F2 = np.empty(Nr, dtype=np.float64)
    _update_edge_anchored_integral(F2, FF_r, accumulator, weights, (R0 * B0) ** 2, 2.0)
    if np.any(F2 <= 0.0) or not np.all(np.isfinite(F2)):
        raise ValueError("integrated F2 must remain finite and positive")
    # Signed F anchored by F_edge = R0*B0 (plan Step 4): never read
    # Equilibrium.F, whose square root is sign-regularized to positive.
    F = np.empty(Nr, dtype=np.float64)
    _update_signed_sqrt(F, F2, R0 * B0)

    P_r = np.empty(Nr, dtype=np.float64)
    _update_scaled_product(P_r, result.P_psi, psi_r, 1.0)
    P = np.empty(Nr, dtype=np.float64)
    _update_edge_anchored_integral(P, P_r, accumulator, weights, P0, 1.0)
    P_axis = float(P0 - np.tensordot(weights, P_r, axes=(0, 0)))

    V_r = np.empty(Nr, dtype=np.float64)
    _update_poloidal_integral(V_r, jacobian, R, 4.0 * np.pi * np.pi / Nt)

    G_fields = np.empty((3, Nr, Nt), dtype=np.float64)
    invalid = _update_G_fields(
        G_fields,
        R,
        JdivR,
        gttdivJR,
        gttdivJR_r,
        grtdivJR_t,
        result.FF_psi,
        result.P_psi,
        psi_r,
        result.psi_rr,
    )
    if invalid:
        raise ValueError(f"GS fields contain a non-finite value at radial index {invalid - 1}")
    Ip = float(_integrate_G1(G_fields[0], weights) / -MU0)

    q = np.empty(Nr, dtype=np.float64)
    invalid = _update_q(q, F, Ln_r, psi_r, r)
    if invalid:
        raise ValueError(f"q is non-finite at radial index {invalid - 1}")

    betat = float(_update_betat(P, V_r, weights, B0))

    R_axis = float(R0 + a * _evaluate_profile_at_axis(result.h_coeffs))
    Z_axis = float(Z0 + a * _evaluate_profile_at_axis(result.v_coeffs))
    axis_weights = np.zeros(Nr, dtype=np.float64)
    axis_weights[0] = 1.0
    F_axis = float(axis_weights @ F)
    B_axis = F_axis / R_axis
    if not np.isfinite(B_axis) or B_axis == 0.0:
        raise ValueError("B_axis must be finite and nonzero")

    return Diagnostics(
        R_axis=R_axis,
        Z_axis=Z_axis,
        B_axis=B_axis,
        P_axis=P_axis,
        qmin=float(np.min(q)),
        q0=float(q[0]),
        q95=float(np.interp(0.95, psin, q)),
        Ip=Ip,
        betat=betat,
        Wth=float(1.5 * np.tensordot(weights, P * V_r, axes=(0, 0))),
        psi_lcfs=float(psi[-1]),
        q=q,
        psin=psin,
        psi=psi,
        F=F,
        P=P,
        V_r=V_r,
    )
