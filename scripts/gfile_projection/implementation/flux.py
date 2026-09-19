"""Bilinear flux-map residual and anchored radial degrees of freedom."""

import numpy as np
from numba import njit

from .geometry import local_profiles, profile_basis, surface_point


@njit(cache=True, nogil=True)
def interpolate(grid, Rmin, Zmin, dR, dZ, R, Z):
    x = (R - Rmin) / dR
    y = (Z - Zmin) / dZ
    if (
        not np.isfinite(x)
        or not np.isfinite(y)
        or x < 0
        or x > grid.shape[0] - 1
        or y < 0
        or y > grid.shape[1] - 1
    ):
        return np.nan, np.nan, np.nan
    i = min(int(x), grid.shape[0] - 2)
    j = min(int(y), grid.shape[1] - 2)
    x -= i
    y -= j
    v00, v01, v10, v11 = grid[i, j], grid[i, j + 1], grid[i + 1, j], grid[i + 1, j + 1]
    value = (1 - x) * ((1 - y) * v00 + y * v01) + x * ((1 - y) * v10 + y * v11)
    dx = ((1 - y) * (v10 - v00) + y * (v11 - v01)) / dR
    dy = ((1 - x) * (v01 - v00) + x * (v11 - v10)) / dZ
    return value, dx, dy


@njit(cache=True, nogil=True)
def setup(r, counts, powers):
    values = profile_basis(r, counts, powers)
    geo_count = np.sum(counts) - int(counts[0] > 0) - int(counts[1] > 0)
    family = np.empty(geo_count, dtype=np.int64)
    orders = np.empty(geo_count, dtype=np.int64)
    idx = 0
    for k in range(len(counts)):
        for degree in range(int(k < 2), counts[k]):
            family[idx] = k
            orders[idx] = degree
            idx += 1
    design = np.empty((len(r), geo_count))
    flux = np.empty((len(r), max(counts)))
    for i in range(len(r)):
        x = 2 * r[i] ** 2 - 1
        T = np.empty(max(counts))
        T[0] = 1.0
        if len(T) > 1:
            T[1] = x
        for degree in range(2, len(T)):
            T[degree] = 2 * x * T[degree - 1] - T[degree - 2]
        for degree in range(len(T)):
            flux[i, degree] = r[i] ** 2 * (1 - r[i] ** 2) * T[degree]
        for j in range(geo_count):
            k, degree = family[j], orders[j]
            design[i, j] = values[i, k, degree]
            if k < 2:
                design[i, j] -= (-1.0) ** degree * values[i, k, 0]
    return values, design, flux, family, orders


@njit(cache=True, nogil=True)
def system(
    grid,
    Rmin,
    Zmin,
    dR,
    dZ,
    boundary,
    coefs,
    beta,
    r,
    theta,
    values,
    design,
    flux,
    family,
    powers,
    nc,
    ns,
    derivative,
):
    raw = local_profiles(boundary, coefs, r, values, powers, nc, ns)
    count = len(family) + len(beta)
    matrix = np.empty((len(r) * len(theta), count)) if derivative else np.empty((0, 0))
    residual = np.empty(len(r) * len(theta))
    cost = 0.0
    for i in range(len(r)):
        source = r[i] ** 2 + flux[i] @ beta
        for t in range(len(theta)):
            idx = i * len(theta) + t
            R, Z, _, _, _, _, _, se = surface_point(raw[i], theta[t], nc, ns)
            s, sR, sZ = interpolate(grid, Rmin, Zmin, dR, dZ, R, Z)
            value = s - source
            residual[idx] = value
            cost += value * value
            if derivative:
                ct, st = np.cos(theta[t]), np.sin(theta[t])
                local = np.empty(len(powers))
                local[0] = boundary[2] * sR
                local[1] = boundary[2] * sZ
                local[2] = -boundary[2] * r[i] * st * sZ
                local[3] = -boundary[2] * r[i] * se * sR
                before_c, value_c, before_s, value_s = 1.0, ct, 0.0, st
                for m in range(1, max(nc - 1, ns) + 1):
                    if m < nc:
                        local[3 + m] = local[3] * value_c
                    if m <= ns:
                        local[3 + nc + m - 1] = local[3] * value_s
                    before_c, value_c = value_c, 2 * ct * value_c - before_c
                    before_s, value_s = value_s, 2 * ct * value_s - before_s
                for j in range(len(family)):
                    matrix[idx, j] = local[family[j]] * design[i, j]
                for j in range(len(beta)):
                    matrix[idx, len(family) + j] = -flux[i, j]
    return matrix, residual, cost


def flux_coefficients(r, levels, count):
    T = np.polynomial.chebyshev.chebvander(2 * r * r - 1, count - 1)
    solution, _, rank, _ = np.linalg.lstsq(
        (r * r * (1 - r * r))[:, None] * T, levels - r * r, rcond=None
    )
    if rank != count:
        raise ValueError("GEQDSK initialization: rank-deficient flux basis")
    return solution


def flux_derivative(beta, r):
    T = np.polynomial.chebyshev.chebval(2 * r * r - 1, beta)
    D = np.polynomial.chebyshev.chebval(2 * r * r - 1, np.polynomial.chebyshev.chebder(beta))
    return 2 * r * (1 + (1 - 2 * r * r) * T + 2 * r * r * (1 - r * r) * D)
