"""Lobatto quadrature and spectral calculus helpers for the joint GEQDSK fit.

Lifted verbatim from ``fusionprime.base.model.geometry`` (88d0fa0):
``_golub_welsch_lobatto`` (:713) through ``_immutable_operator`` (:798),
including ``_barycentric_log_weights`` (:761).
VEQDB retains these helpers to preserve the compact codec's numerical behavior.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.linalg import eigh_tridiagonal


def _golub_welsch_lobatto(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Construct two-endpoint Gauss-Lobatto nodes and weights."""

    diagonal = np.zeros(n, dtype=np.float64)
    off_diagonal = np.empty(n - 1, dtype=np.float64)
    indices = np.arange(1, n - 1, dtype=np.float64)
    off_diagonal[:-1] = indices / np.sqrt(4.0 * indices * indices - 1.0)
    off_diagonal[-1] = np.sqrt((n - 1.0) / (2.0 * n - 3.0))
    nodes, vectors = eigh_tridiagonal(diagonal, off_diagonal)
    return nodes, 2.0 * vectors[0] ** 2


def _spectral_differentiator(nodes: np.ndarray) -> np.ndarray:
    """Construct the barycentric spectral first-derivative matrix."""

    signs, log_weights = _barycentric_log_weights(nodes)
    difference = nodes[:, None] - nodes[None, :]
    np.fill_diagonal(difference, 1.0)
    ratio = np.exp(np.clip(log_weights[None, :] - log_weights[:, None], -700.0, 700.0))
    matrix = signs[None, :] * signs[:, None] * ratio / difference
    np.fill_diagonal(matrix, 0.0)
    matrix[np.diag_indices_from(matrix)] = -np.sum(matrix, axis=1)
    return matrix


def _spectral_accumulator(nodes: np.ndarray) -> np.ndarray:
    """Construct exact nodal prefix integrals through a shifted Legendre basis."""

    n = len(nodes)
    x = 2.0 * nodes - 1.0
    full_basis = np.polynomial.legendre.legvander(x, n)
    basis = full_basis[:, :n]
    primitive = np.zeros((n, n), dtype=np.float64)
    primitive[:, 0] = 0.5 * (x + 1.0)
    degrees = np.arange(1, n, dtype=np.float64)
    primitive[:, 1:] = (
        0.5 * (full_basis[:, 2 : n + 1] - full_basis[:, : n - 1]) / (2.0 * degrees + 1.0)
    )
    lower_basis = np.polynomial.legendre.legvander(np.array([-1.0]), n)[0]
    lower = np.zeros(n, dtype=np.float64)
    lower[1:] = 0.5 * (lower_basis[2 : n + 1] - lower_basis[: n - 1]) / (2.0 * degrees + 1.0)
    rhs = (primitive - lower[None, :]).T
    condition = np.linalg.cond(basis)
    if np.isfinite(condition) and condition <= 1.0 / np.sqrt(np.finfo(np.float64).eps):
        return np.linalg.solve(basis.T, rhs).T
    return np.linalg.lstsq(basis.T, rhs, rcond=None)[0].T


def _barycentric_log_weights(nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute stable barycentric signs and logarithmic magnitudes."""

    difference = nodes[:, None] - nodes[None, :]
    mask = ~np.eye(len(nodes), dtype=bool)
    if np.any(np.isclose(difference[mask], 0.0)):
        raise ValueError("radial nodes must be distinct")
    absolute = np.abs(difference)
    signs = np.sign(difference)
    np.fill_diagonal(absolute, 1.0)
    np.fill_diagonal(signs, 1.0)
    return np.prod(signs, axis=1), -np.sum(np.log(absolute), axis=1)


@lru_cache
def _lobatto_quadrature(count: int) -> tuple[np.ndarray, np.ndarray]:
    """Build the fixed rule once per resolution, backed by immutable storage."""

    nodes, weights = _golub_welsch_lobatto(count)
    nodes = 0.5 * (nodes + 1.0)
    nodes[0], nodes[-1] = 0.0, 1.0
    return _immutable_operator(nodes), _immutable_operator(0.5 * weights)


@lru_cache
def _spectral_calculus(count: int) -> tuple[np.ndarray, np.ndarray]:
    """Share only resolution-dependent calculus, never physical field caches."""

    nodes, _ = _lobatto_quadrature(count)
    return (
        _immutable_operator(_spectral_accumulator(nodes)),
        _immutable_operator(_spectral_differentiator(nodes)),
    )


def _immutable_operator(values: np.ndarray) -> np.ndarray:
    # A bytes owner prevents callers from re-enabling writes on a shared array.
    return np.frombuffer(values.tobytes(), dtype=np.float64).reshape(values.shape)
