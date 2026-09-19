"""Joint LCFS, interior geometry and flux fitting; private GEQDSK implementation.

This is a least-squares projection of a file, not a Grad--Shafranov solve.
Only a converged, full-rank, geometrically and physically accepted candidate
is returned as a frozen State. No public solve-map derivative is provided.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, TypeVar

import numpy as np

from ....state import Equilibrium
from ....state.geometry import _lobatto_quadrature, _spectral_calculus
from ..payload import Geqdsk
from ..vacuum import inside_polygon
from .flux import flux_coefficients, flux_derivative, setup
from .geometry import local_profiles, profile_basis, projection_diagnostics, radial_derivatives
from .initial import _sample_contour, extract_surfaces, fit_radial, fit_surface
from .kernel import fit, joint_system

_Equilibrium = TypeVar("_Equilibrium", bound=Equilibrium)


@dataclass(frozen=True, slots=True)
class Settings:
    Nr: int
    Nt: int
    c_order: int
    s_order: int
    K_max: int | None
    cosine_count: int
    sine_count: int
    nc: int
    ns: int
    counts: np.ndarray
    powers: np.ndarray
    boundary_indices: np.ndarray


def settings(
    topology: Any, Nr: int, Nt: int, c_order: int, s_order: int, K_max: int | None
) -> Settings:
    for name, value, minimum in (
        ("Nr", Nr, 4),
        ("Nt", Nt, 4),
        ("c_order", c_order, 0),
        ("s_order", s_order, 0),
    ):
        if type(value) is not int:
            raise TypeError(f"GEQDSK configuration: {name} must be int")
        if value < minimum:
            raise ValueError(f"GEQDSK configuration: {name} must be at least {minimum}")
    if K_max is not None and (type(K_max) is not int or K_max < 2):
        raise ValueError("GEQDSK configuration: K_max must be None or an int >= 2")
    if topology is None:
        topology = dict(h_count=4, v_count=4, kappa_count=4, c_counts=(3, 2, 2), s_counts=(3, 2, 2))
    if not isinstance(topology, Mapping):
        raise TypeError("GEQDSK configuration: topology must be a mapping")
    if set(topology) != {"h_count", "v_count", "kappa_count", "c_counts", "s_counts"}:
        raise ValueError(
            "GEQDSK configuration: topology requires exactly "
            "h_count, v_count, kappa_count, c_counts, s_counts"
        )
    axes = []
    for name in ("c_counts", "s_counts"):
        values = topology[name]
        if not isinstance(values, (list, tuple)):
            raise TypeError(f"GEQDSK configuration: {name} must be a list or tuple of int")
        axes.append(tuple(values))
    c, s = axes
    basic = (topology["h_count"], topology["v_count"], topology["kappa_count"])
    for value in (*basic, *c, *s):
        if type(value) is not int:
            raise TypeError("GEQDSK configuration: coefficient counts must be int")
        if value < 0:
            raise ValueError("GEQDSK configuration: coefficient counts must be nonnegative")
    if basic[0] == 0 or basic[1] == 0:
        raise ValueError(
            "GEQDSK configuration: joint fitting requires positive h_count and v_count"
        )
    width = max(*basic, *c, *s)
    # s_beta has degree 2*width+2: its derivative must fit the native polynomial grid.
    if Nr < 2 * width + 2:
        raise ValueError(f"GEQDSK configuration: Nr must be >= 2*max(counts)+2={2 * width + 2}")
    nc, ns = max(c_order + 1, len(c)), max(s_order, len(s))
    if Nt < 2 * max(nc - 1, ns) + 1:
        raise ValueError(
            "GEQDSK configuration: Nt must resolve all boundary and interior harmonics"
        )
    counts = np.array(
        (*basic, *c, *((0,) * (nc - len(c))), *s, *((0,) * (ns - len(s)))), dtype=np.int64
    )
    powers = np.r_[0, 0, 0, np.arange(nc), np.arange(1, ns + 1)].astype(np.int64)
    if K_max is not None:
        powers = np.minimum(powers, K_max)
    indices = np.r_[np.arange(4 + c_order + 1), np.arange(4 + nc, 4 + nc + s_order)].astype(
        np.int64
    )
    return Settings(
        Nr, Nt, c_order, s_order, K_max, len(c), len(s), nc, ns, counts, powers, indices
    )


def validate_payload(payload: Geqdsk) -> None:
    if not isinstance(payload, Geqdsk):
        raise TypeError("GEQDSK input: gfile must be Geqdsk")
    if min(payload.NR, payload.NZ) < 4 or payload.boundary.shape[0] < 4:
        raise ValueError("GEQDSK input: at least a 4 x 4 grid and four LCFS points are required")
    if payload.psi_axis == payload.psi_bound or payload.F[-1] == 0.0:
        raise ValueError("GEQDSK input: nonzero flux span and edge F are required")
    if not (
        0 < payload.Rmin < payload.Raxis < payload.Rmax
        and payload.Zmin < payload.Zaxis < payload.Zmax
    ):
        raise ValueError("GEQDSK input: positive-R rectangle must contain the magnetic axis")
    b = payload.boundary
    if (
        np.any(b[:, 0] < payload.Rmin)
        or np.any(b[:, 0] > payload.Rmax)
        or np.any(b[:, 1] < payload.Zmin)
        or np.any(b[:, 1] > payload.Zmax)
    ):
        raise ValueError("GEQDSK input: LCFS must lie within the rectangular grid")

    contained, on_edge = inside_polygon(np.array([payload.Raxis]), np.array([payload.Zaxis]), b)
    if not contained[0] or on_edge[0]:
        raise ValueError("GEQDSK input: LCFS must enclose the magnetic axis")


def initialize(payload: Geqdsk, config: Settings):
    levels, contours = extract_surfaces(payload, config)
    contours = [_sample_contour(c, max(64, 2 * config.Nt)) for c in contours]
    raw = np.array([fit_surface(c, config.nc - 1, config.ns) for c in contours])
    boundary = raw[-1].copy()
    # The common numerical harmonic axes do not introduce new boundary freedom.
    inactive = np.ones(boundary.size, dtype=bool)
    inactive[config.boundary_indices] = False
    boundary[inactive] = 0.0
    r = raw[:, 2] / boundary[2]
    if np.any(np.diff(r) <= 0) or r[0] <= 0:
        raise ValueError("non-nested geometric radii in contour initial guess")
    values = np.column_stack(
        (
            (raw[:, 0] - boundary[0]) / boundary[2],
            (raw[:, 1] - boundary[1]) / boundary[2],
            raw[:, 3:],
        )
    )
    edge = np.r_[0.0, 0.0, boundary[3:]]
    axis = np.array([payload.Raxis, payload.Zaxis])
    coefficients = fit_radial(
        r, values, edge, config.counts, config.powers, (axis - boundary[:2]) / boundary[2]
    )
    beta = flux_coefficients(r, levels, int(max(config.counts)))
    return boundary, coefficients, beta, contours[-1], axis


def normalized_jacobian(boundary, coefficients, r, theta, config):
    values = profile_basis(r, config.counts, config.powers)
    raw = local_profiles(boundary, coefficients, r, values, config.powers, config.nc, config.ns)
    dr = radial_derivatives(boundary, coefficients, r, config.powers, config.nc, config.ns)
    phase = np.array(
        [np.cos(m * theta) for m in range(config.nc)]
        + [np.sin(m * theta) for m in range(1, config.ns + 1)]
    )
    phase_t = np.array(
        [-m * np.sin(m * theta) for m in range(config.nc)]
        + [m * np.cos(m * theta) for m in range(1, config.ns + 1)]
    )
    eta, eta_t, eta_r = theta + raw[:, 4:] @ phase, 1 + raw[:, 4:] @ phase_t, dr[:, 4:] @ phase
    a, ar, kappa = boundary[2], boundary[2] * r[:, None], raw[:, 3, None]
    Rt = -ar * np.sin(eta) * eta_t
    Zt = -ar * kappa * np.cos(theta)
    Rr = dr[:, 0, None] + a * np.cos(eta) - ar * np.sin(eta) * eta_r
    Zr = dr[:, 1, None] - (a * kappa + ar * dr[:, 3, None]) * np.sin(theta)
    R = raw[:, 0, None] + ar * np.cos(eta)
    return (Rt * Zr - Rr * Zt) / (a * a * r[:, None]), R, kappa


def accept_geometry(boundary, coefficients, beta, r, config):
    check_r = np.unique(np.r_[r[1:], np.linspace(0, 1, max(129, 4 * config.Nr + 1))[1:]])
    theta = np.unique(
        np.r_[
            np.linspace(0, 2 * np.pi, max(256, 8 * config.Nt), endpoint=False),
            np.linspace(0, 2 * np.pi, config.Nt, endpoint=False),
        ]
    )
    J, R, kappa = normalized_jacobian(boundary, coefficients, check_r, theta, config)
    derivative = flux_derivative(beta, check_r)
    if not np.all(np.isfinite(J)) or np.any(J <= 0):
        raise RuntimeError(
            f"GEQDSK geometry acceptance: nonpositive or nonfinite J; min_J={np.min(J):.9g}"
        )
    if not np.all(np.isfinite(R)) or np.any(R <= 0) or np.any(kappa <= 0):
        raise RuntimeError("GEQDSK geometry acceptance: invalid radius or elongation")
    if not np.all(np.isfinite(derivative)) or np.any(derivative <= 0):
        raise RuntimeError("GEQDSK flux acceptance: nonmonotonic normalized flux")


def materialize(payload, boundary, coefficients, beta, config, state_type):
    r, weights = _lobatto_quadrature(config.Nr)
    accumulator, _ = _spectral_calculus(config.Nr)
    span = payload.psi_bound - payload.psi_axis
    psi_r = span * flux_derivative(beta, r)
    integral = float(weights @ psi_r)
    if not np.isfinite(integral) or integral == 0:
        raise RuntimeError("invalid poloidal flux integral")
    psi_r *= span / integral
    psin = accumulator @ psi_r / span
    psin[0], psin[-1] = 0.0, 1.0
    if np.any(np.diff(psin) <= 0):
        raise RuntimeError("nonmonotonic native poloidal flux")
    nodes = np.linspace(0, 1, payload.NR)
    candidate = state_type(
        Nr=config.Nr,
        Nt=config.Nt,
        K_max=config.K_max,
        R0=boundary[0],
        Z0=boundary[1],
        a=boundary[2],
        kappa_lcfs=boundary[3],
        c_lcfs=boundary[4 : 5 + config.c_order],
        s_lcfs=boundary[4 + config.nc : 4 + config.nc + config.s_order],
        h_coeffs=coefficients[0],
        v_coeffs=coefficients[1],
        kappa_coeffs=coefficients[2],
        c_coeffs=coefficients[3 : 3 + config.cosine_count],
        s_coeffs=coefficients[3 + config.nc : 3 + config.nc + config.sine_count],
        B0=float(payload.F[-1] / boundary[0]),
        P0=float(payload.P[-1]),
        psi_r=psi_r,
        P_psi=np.interp(psin, nodes, payload.P_psi),
        FF_psi=np.interp(psin, nodes, payload.FF_psi),
    )
    for name in (
        "J",
        "rho",
        "F",
        "P",
        "q",
        "Ip",
        "R",
        "Z",
        "gm1",
        "gm2",
        "gm9",
        "ftrap",
        "Itor",
        "jtor",
        "jpara",
        "jtotal",
        "S",
        "V",
        "A",
        "Phi",
        "Phi_rho",
    ):
        if not np.all(np.isfinite(getattr(candidate, name))):
            raise RuntimeError(f"nonfinite physical field {name}")
    return candidate.freeze()


def from_geqdsk(
    payload: Geqdsk,
    *,
    state_type: type[_Equilibrium],
    topology: Mapping[str, int | Sequence[int]] | None = None,
    Nr: int = 32,
    Nt: int = 32,
    c_order: int = 5,
    s_order: int = 5,
    K_max: int | None = None,
) -> _Equilibrium:
    config = settings(topology, Nr, Nt, c_order, s_order, K_max)
    validate_payload(payload)
    try:
        boundary, coefficients, beta, points, axis = initialize(payload, config)
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        raise RuntimeError(f"GEQDSK initialization failed: {error}") from error
    r = _lobatto_quadrature(Nr)[0]
    theta = np.linspace(0, 2 * np.pi, max(64, 2 * Nt), endpoint=False)
    grid = (payload.psi - payload.psi_axis) / (payload.psi_bound - payload.psi_axis)
    head = (grid, payload.Rmin, payload.Zmin, payload.dR, payload.dZ)
    try:
        b, c, beta, cost, iterations, status, angles = fit(
            *head,
            boundary,
            coefficients,
            beta,
            r,
            theta,
            config.counts,
            config.powers,
            config.nc,
            config.ns,
            points,
            axis,
            config.boundary_indices,
        )
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        raise RuntimeError(f"GEQDSK joint solve failed: {error}") from error
    detail = f"status={status}, iteration={iterations}, residual={np.sqrt(cost):.9g}"
    if status not in ("gradient", "relative_cost") or not np.isfinite(cost):
        raise RuntimeError(f"GEQDSK joint solve rejected: {detail}")
    step, failures = projection_diagnostics(points, b, angles, config.nc, config.ns)
    if failures:
        raise RuntimeError(
            f"GEQDSK LCFS projection rejected: failures={failures}, step={step:.9g}; {detail}"
        )
    values, design, flux, family, _ = setup(r, config.counts, config.powers)
    A, residual, _, _ = joint_system(
        *head,
        b,
        c,
        beta,
        r,
        theta,
        values,
        design,
        flux,
        family,
        config.powers,
        config.nc,
        config.ns,
        points,
        angles,
        axis,
        2 / boundary[2],
        1.0,
        True,
    )
    columns = np.r_[config.boundary_indices, np.arange(len(b), A.shape[1])]
    A = A[:, columns]
    scales = np.linalg.norm(A, axis=0)
    if (
        not np.all(np.isfinite(A))
        or not np.all(np.isfinite(residual))
        or not np.all(np.isfinite(scales))
        or np.any(scales == 0)
    ):
        raise RuntimeError(f"GEQDSK rank acceptance: zero or nonfinite Jacobian column; {detail}")
    try:
        singular = np.linalg.svd(A / scales, compute_uv=False)
    except np.linalg.LinAlgError as error:
        raise RuntimeError(f"GEQDSK rank acceptance failed: {error}; {detail}") from error
    rank = int(np.sum(singular > singular[0] * max(A.shape) * np.finfo(float).eps))
    if rank != A.shape[1]:
        raise RuntimeError(f"GEQDSK rank acceptance: rank={rank}/{A.shape[1]}; {detail}")
    accept_geometry(b, c, beta, r, config)
    try:
        return materialize(payload, b, c, beta, config, state_type)
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        raise RuntimeError(f"GEQDSK physical reconstruction rejected: {error}; {detail}") from error
