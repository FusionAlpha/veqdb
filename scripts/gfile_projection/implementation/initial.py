"""Contour/phase initialization only; final LCFS belongs to the joint solve."""

import math
from typing import cast

import numpy as np
from contourpy import contour_generator
from numba import njit
from scipy.interpolate import RegularGridInterpolator

from ....plasma.boundary import _phase_design_matrix, _weighted_phase_qr
from ..conversion import _axis_contour, _rectangular_nodes
from ..payload import Geqdsk


def extract_surfaces(payload: Geqdsk, config) -> tuple[np.ndarray, list[np.ndarray]]:
    """Extract ordered contours, with edge clustering in normalized poloidal flux."""
    rn, zn = _rectangular_nodes(payload)
    generator = contour_generator(x=rn, y=zn, z=payload.psi.T, line_type="Separate")
    count = max(config.Nr, 2 * int(np.max(config.counts)) + 1)
    axis = np.array([payload.Raxis, payload.Zaxis])
    sampled_axis = float(RegularGridInterpolator((rn, zn), payload.psi)(axis[None, :])[0])
    resolved_axis = max(
        0.0, (sampled_axis - payload.psi_axis) / (payload.psi_bound - payload.psi_axis)
    )
    if resolved_axis >= 1.0:
        raise ValueError("GEQDSK rectangular flux does not resolve an interior magnetic axis")
    levels = (
        resolved_axis + (1.0 - resolved_axis) * np.sin(np.linspace(0.0, np.pi / 2.0, count)) ** 2
    )
    levels[-1] = 1.0
    contours = []
    for level in levels[1:-1]:
        psi = float(payload.psi_axis + level * (payload.psi_bound - payload.psi_axis))
        contours.append(
            _axis_contour(cast(list[np.ndarray], generator.lines(psi)), axis, level=psi)
        )
    contours.append(payload.boundary)
    # Arc sampling bounds numerical fit cost independently of file resolution.
    # Retain extrema exactly: they define geometric radius, centre and elongation.
    return levels[1:], [_sample_contour(c, 8 * config.Nt) for c in contours]


def _sample_contour(points: np.ndarray, count: int) -> np.ndarray:
    if np.array_equal(points[0], points[-1]):
        points = points[:-1]
    # Canonical arc origin/orientation makes resampling independent of the
    # file's arbitrary starting vertex and traversal direction.
    area = np.sum(
        points[:, 0] * np.roll(points[:, 1], -1) - points[:, 1] * np.roll(points[:, 0], -1)
    )
    if area > 0.0:
        points = points[::-1]
    points = np.roll(points, -int(np.argmax(points[:, 0])), axis=0)
    points = np.vstack((points, points[0]))
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    keep = np.r_[True, lengths > 0.0]
    points = points[keep]
    if points.shape[0] < 4:
        raise ValueError("GEQDSK geometry fit requires nondegenerate closed contours")
    arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    extrema = [
        np.argmin(points[:, 0]),
        np.argmax(points[:, 0]),
        np.argmin(points[:, 1]),
        np.argmax(points[:, 1]),
    ]
    sample = np.unique(np.r_[np.linspace(0.0, arc[-1], count, endpoint=False), arc[extrema]])
    return np.column_stack([np.interp(sample, arc, points[:, i]) for i in (0, 1)])


@njit(cache=True, nogil=True)
def fit_surface(points: np.ndarray, c_order: int, s_order: int) -> np.ndarray:
    """Infer phase by contour topology, then solve the weighted phase system.

    Clockwise theta follows Z = Zc - a*kappa*sin(theta). The two inverse
    trigonometric branches are selected by extrema indices, independently of
    sample spacing. This avoids a phase jump when points cluster at an extremum.
    """
    R, Z = points[:, 0].copy(), points[:, 1].copy()
    n = R.size
    area = 0.0
    for i in range(n):
        area += R[i] * Z[(i + 1) % n] - Z[i] * R[(i + 1) % n]
    if area > 0.0:
        R, Z = R[::-1].copy(), Z[::-1].copy()
    R0, Z0 = (np.max(R) + np.min(R)) / 2, (np.max(Z) + np.min(Z)) / 2
    a = (np.max(R) - np.min(R)) / 2
    height = (np.max(Z) - np.min(Z)) / 2
    if a <= 0.0 or height <= 0.0:
        raise ValueError("GEQDSK geometry fit found a degenerate contour")
    bottom, top, left, right = np.argmin(Z), np.argmax(Z), np.argmin(R), np.argmax(R)
    theta, eta = np.empty(n), np.empty(n)
    for i in range(n):
        angle = math.asin(min(1.0, max(-1.0, -(Z[i] - Z0) / height)))
        theta[i] = (
            math.pi - angle if (i - bottom) % n <= (top - bottom) % n else 2 * math.pi + angle
        )
        angle = math.acos(min(1.0, max(-1.0, (R[i] - R0) / a)))
        phase = angle if (i - right) % n <= (left - right) % n else 2 * math.pi - angle
        eta[i] = theta[i] + (phase - theta[i] + math.pi) % (2 * math.pi) - math.pi
    matrix = _phase_design_matrix(theta, c_order, s_order)
    coeff = _weighted_phase_qr(theta, eta, matrix, a)
    result = np.empty(4 + coeff.size)
    result[:4] = np.array([R0, Z0, a, height / a])
    result[4:] = coeff
    return result


@njit(cache=True, nogil=True)
def fit_radial(
    r: np.ndarray,
    values: np.ndarray,
    edge: np.ndarray,
    counts: np.ndarray,
    powers: np.ndarray,
    axis: np.ndarray,
) -> np.ndarray:
    """Solve canonical enveloped Chebyshev profiles, with exact magnetic axis."""
    width = max(1, np.max(counts))
    coefficients = np.zeros((counts.size, width))
    for k in range(counts.size):
        count = counts[k]
        if count == 0:
            continue
        matrix = np.empty((r.size, count))
        target = np.empty(r.size)
        for i in range(r.size):
            radial = r[i] ** powers[k]
            x = 2 * r[i] ** 2 - 1
            matrix[i, 0] = radial * (1 - r[i] ** 2)
            if count > 1:
                matrix[i, 1] = matrix[i, 0] * x
            for j in range(2, count):
                matrix[i, j] = 2 * x * matrix[i, j - 1] - matrix[i, j - 2]
            target[i] = values[i, k] - radial * edge[k]
        if k < 2:
            # Eliminate coefficient zero using sum((-1)**j * xj) = axis shift.
            target -= matrix[:, 0] * axis[k]
            for j in range(1, count):
                matrix[:, j] -= (-1.0) ** j * matrix[:, 0]
            if count > 1:
                solution, _, rank, _ = np.linalg.lstsq(matrix[:, 1:].copy(), target, rcond=-1.0)
                if rank != count - 1:
                    raise np.linalg.LinAlgError("GEQDSK radial fit is rank deficient")
                coefficients[k, 1:count] = solution
            coefficients[k, 0] = axis[k]
            for j in range(1, count):
                coefficients[k, 0] -= (-1.0) ** j * coefficients[k, j]
        else:
            solution, _, rank, _ = np.linalg.lstsq(matrix, target, rcond=-1.0)
            if rank != count:
                raise np.linalg.LinAlgError("GEQDSK radial fit is rank deficient")
            coefficients[k, :count] = solution
    return coefficients
