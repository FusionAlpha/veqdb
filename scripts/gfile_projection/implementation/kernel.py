"""One Numba solve for LCFS, interior geometry and flux roots.

Loss = mean((s_grid(X)-s(r))**2) + weight*mean((2*d_LCFS/a_initial)**2).
The fixed initial halfwidth normalizes units; it is not an optimized weight.
Both h and v must be represented so the measured magnetic axis stays fixed.
"""

import numpy as np
from numba import njit

from .flux import interpolate, setup, system
from .geometry import (
    initial_angles,
    local_profiles,
    normal_system,
    project,
    surface_point,
)


@njit(cache=True, nogil=True)
def anchor(boundary, coefficients, axis):
    for k in range(2):
        coefficients[k, 0] = (axis[k] - boundary[k]) / boundary[2]
        for j in range(1, coefficients.shape[1]):
            coefficients[k, 0] -= (-1.0) ** j * coefficients[k, j]


@njit(cache=True, nogil=True)
def joint_system(
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
    points,
    angles,
    axis,
    scale,
    weight,
    derivative,
):
    A, residual, cost = system(
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
    )
    projected, edge_cost = project(points, boundary, angles, nc, ns)
    edge_A, edge_residual = normal_system(points, boundary, projected, nc, ns)
    nf, ne, nb = len(residual), len(points), len(boundary)
    wf, we = 1 / np.sqrt(nf), scale * np.sqrt(weight / ne)
    target = np.concatenate((residual * wf, edge_residual * we))
    matrix = np.zeros((nf + ne, nb + len(family) + len(beta))) if derivative else np.empty((0, 0))
    if derivative:
        matrix[:nf, nb:] = A * wf
        matrix[nf:, :nb] = edge_A * we
        raw = local_profiles(boundary, coefs, r, values, powers, nc, ns)
        for i in range(len(r)):
            for t in range(len(theta)):
                idx = i * len(theta) + t
                R, Z, _, _, _, _, _, se = surface_point(raw[i], theta[t], nc, ns)
                _, sR, sZ = interpolate(grid, Rmin, Zmin, dR, dZ, R, Z)
                st = np.sin(theta[t])
                matrix[idx, 0] = sR * r[i] ** 2 * wf
                matrix[idx, 1] = sZ * r[i] ** 2 * wf
                # Remove the axis-anchor's 1/a dependence from the scale derivative.
                dRa = (R - boundary[0] - (1 - r[i] ** 2) * (axis[0] - boundary[0])) / boundary[2]
                dZa = (Z - boundary[1] - (1 - r[i] ** 2) * (axis[1] - boundary[1])) / boundary[2]
                matrix[idx, 2] = (sR * dRa + sZ * dZa) * wf
                matrix[idx, 3] = -sZ * boundary[2] * r[i] * st * wf
                phase = -sR * boundary[2] * r[i] * se * wf
                matrix[idx, 4] = phase * r[i] ** powers[3]
                for m in range(1, nc):
                    matrix[idx, 4 + m] = phase * np.cos(m * theta[t]) * r[i] ** powers[3 + m]
                for m in range(1, ns + 1):
                    matrix[idx, 3 + nc + m] = (
                        phase * np.sin(m * theta[t]) * r[i] ** powers[2 + nc + m]
                    )
    return matrix, target, cost / nf + weight * scale**2 * edge_cost / ne, projected


@njit(cache=True, nogil=True)
def fit(
    grid,
    Rmin,
    Zmin,
    dR,
    dZ,
    boundary_initial,
    initial,
    beta_initial,
    r,
    theta,
    counts,
    powers,
    nc,
    ns,
    points,
    axis,
    boundary_indices,
    weight=1.0,
    relative_tolerance=1e-6,
    max_iterations=80,
):
    if counts[0] == 0 or counts[1] == 0:
        raise ValueError("joint LCFS kernel requires nonzero h and v counts")
    if not np.isfinite(weight) or weight <= 0:
        raise ValueError("LCFS weight must be positive and finite")
    values, design, flux, family, orders = setup(r, counts, powers)
    boundary, coefs, beta = boundary_initial.copy(), initial.copy(), beta_initial.copy()
    anchor(boundary, coefs, axis)
    angles = initial_angles(points, boundary)
    scale = 2 / boundary_initial[2]
    damping = 1e-3
    nb = len(boundary_indices)
    columns = np.concatenate(
        (boundary_indices, np.arange(len(boundary), len(boundary) + len(family) + len(beta)))
    )
    cost = np.inf
    iteration = -1
    status = "iteration_limit"
    for iteration in range(max_iterations):
        matrix, target, cost, angles = joint_system(
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
            points,
            angles,
            axis,
            scale,
            weight,
            True,
        )
        matrix = matrix[:, columns].copy()
        if not np.isfinite(cost) or not np.all(np.isfinite(target)):
            return boundary, coefs, beta, cost, iteration + 1, "nonfinite_system", angles
        scales = np.sqrt(np.sum(matrix * matrix, axis=0))
        if not np.all(np.isfinite(scales)) or np.any(scales == 0):
            return boundary, coefs, beta, cost, iteration + 1, "zero_or_nonfinite_column", angles
        matrix /= scales
        H, g = matrix.T @ matrix, matrix.T @ target
        if np.max(np.abs(g)) <= 1e-8 * max(np.sqrt(cost), np.finfo(np.float64).eps):
            status = "gradient"
            break
        accepted, old = False, cost
        for retry in range(12):
            shifted = H.copy()
            for j in range(len(scales)):
                shifted[j, j] += damping
            step = np.linalg.solve(shifted, -g) / scales
            b, c = boundary.copy(), coefs.copy()
            b[boundary_indices] += step[:nb]
            for j in range(len(family)):
                c[family[j], orders[j]] += step[nb + j]
            be = beta + step[nb + len(family) :]
            if b[2] > 0 and b[3] > 0:
                anchor(b, c, axis)
                _, _, new_cost, new_angles = joint_system(
                    grid,
                    Rmin,
                    Zmin,
                    dR,
                    dZ,
                    b,
                    c,
                    be,
                    r,
                    theta,
                    values,
                    design,
                    flux,
                    family,
                    powers,
                    nc,
                    ns,
                    points,
                    angles,
                    axis,
                    scale,
                    weight,
                    False,
                )
                if new_cost < cost:
                    boundary, coefs, beta, cost, angles = b, c, be, new_cost, new_angles
                    damping = max(1e-12, damping / 3)
                    accepted = True
                    break
            damping *= 10
        if not accepted:
            status = "damping_exhausted"
            break
        if old - cost <= relative_tolerance * old:
            status = "relative_cost"
            break
    return boundary, coefs, beta, cost, iteration + 1, status, angles
