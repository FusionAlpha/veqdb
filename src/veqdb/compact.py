"""Encode and decode the frozen VEQDB compact JSON schema.

Public API:
- ``to_compact_json`` -- serialize roots plus display-only diagnostics at the
  requested precision (``None`` = full round-trip fidelity).
- ``from_compact_json`` -- restore a trusted-local record to standalone
  reconstruction roots without requiring the VEQ solver.

Notes:
- Port of ``fusionprime.base.io.compact_equilibrium`` (88d0fa0), via VEQPy
  revision 586afb717099fcc1eb17ac5ffca4a10a95ae6d80: key set and
  ``_JSON_FIELDS`` order, trailing-zero trimming through ``np.uint64`` views
  (``-0.0`` survives), harmonic-row popping with interior ``[0.0]`` rows, the
  ``K_max`` drop rules, the no-second-rounding decimal formatter, and
  ``orjson`` encoding. Diagnostics are computed by ``veqdb._diagnostics``.
- The compact schema stores the *raw* fit roots. Sign-regularized or
  zero-edge-gauge profiles from a solver state are never serialized.
- Restore rebuilds the roots with the record's integrated flux span
  (``psi_lcfs``); the g-file span of the original fit is not representable in
  the schema, so the restored ``psin`` gauge may differ from the pre-write
  state while every stored field and diagnostic round-trips exactly.
  Diagnostics are never read back: they are recomputed from roots on encode.
- The optional ``header`` describes the file; the reader drops it because it
  is metadata, not a reconstruction input.
"""

from __future__ import annotations

from math import isfinite
from typing import Any

import numpy as np
import orjson

from ._diagnostics import compute
from ._result import FitResult, _harmonic_order
from ._spectral import _lobatto_quadrature, _spectral_calculus

__all__ = ["from_compact_json", "to_compact_json"]

_PROFILES = ("P_psi", "FF_psi", "psi_r", "psi_rr")
_COEFFICIENTS = ("h_coeffs", "v_coeffs", "kappa_coeffs", "c_coeffs", "s_coeffs")
_BOUNDARY = ("R0", "Z0", "a", "kappa_lcfs", "c_lcfs", "s_lcfs")
_FIELDS = (
    "header",
    "Nr",
    "Nt",
    *_PROFILES,
    *_COEFFICIENTS,
    "K_max",
    *_BOUNDARY,
    "B0",
    "P0",
)
_JSON_FIELDS = (
    "header",
    "R_axis",
    "Z_axis",
    "B_axis",
    "P_axis",
    "B0",
    "P0",
    "R0",
    "Z0",
    "a",
    "betat",
    "Ip",
    "qmin",
    "q0",
    "q95",
    "Wth",
    "psi_lcfs",
    "kappa_lcfs",
    "c_lcfs",
    "s_lcfs",
    "Nr",
    "Nt",
    "K_max",
    *_COEFFICIENTS,
    *_PROFILES,
)


def to_compact_json(
    result: FitResult,
    *,
    header: str | None = "Equilibrium",
    precision: int | None = 10,
) -> bytes:
    """Encode roots and display-only diagnostics at the requested output precision.

    Defaults mirror the fusionprime ``Equilibrium.to_compact_json`` method
    (``header="Equilibrium"``, ``precision=10``).
    """

    if type(precision) is int:
        if precision < 1:
            raise ValueError("precision must be a positive significant-digit count")
    elif precision is not None:
        raise TypeError("precision must be a positive integer or None")
    data = _pack(result, header=header)
    M_max = _harmonic_order(result.c_lcfs, result.s_lcfs, result.c_coeffs, result.s_coeffs)
    if result.K_max is not None and result.K_max >= M_max:
        # min(m, K_max) equals m for every represented harmonic.
        del data["K_max"]
    for name in ("c_coeffs", "s_coeffs"):
        rows = data[name]
        while rows and rows[-1].size == 0:
            rows.pop()
        # Interior zero rows retain harmonic indices; trailing zero modes do not.
        data[name] = [row if row.size else [0.0] for row in rows]
    diagnostics = compute(result)
    q = diagnostics.q
    data.update(
        R_axis=diagnostics.R_axis,
        Z_axis=diagnostics.Z_axis,
        B_axis=diagnostics.B_axis,
        P_axis=diagnostics.P_axis,
        qmin=float(np.min(q)),
        q0=float(q[0]),
        q95=float(np.interp(0.95, diagnostics.psin, q)),
        Ip=diagnostics.Ip,
        betat=diagnostics.betat,
        Wth=diagnostics.Wth,
        psi_lcfs=diagnostics.psi_lcfs,
    )

    def convert(item: Any, key: str) -> Any:
        if isinstance(item, np.ndarray):
            item = item.tolist()
        if isinstance(item, list):
            return [convert(child, key) for child in item]
        if not isinstance(item, float):
            return item
        if not isfinite(item):
            raise ValueError(f"{key} must be finite for compact JSON export")
        if item == 0.0:
            return orjson.Fragment(b"-0.0" if np.signbit(item) else b"0")
        token = repr(item) if precision is None else format(item, f".{precision}g")
        mantissa, _, exponent = token.partition("e")
        whole, _, fraction = mantissa.lstrip("-").partition(".")
        digits = (whole + fraction).lstrip("0")
        point = len(whole) + int(exponent or "0") - (len(whole + fraction) - len(digits))
        digits = digits.rstrip("0")
        # Move the decimal point in the rounded token, without rounding a second time.
        if point > 3 or point <= -3:
            token = digits[0] + ("." + digits[1:] if len(digits) > 1 else "")
            token += f"e{point - 1}"
        elif point <= 0:
            token = "0." + "0" * -point + digits
        elif point >= len(digits):
            token = digits + "0" * (point - len(digits))
        else:
            token = digits[:point] + "." + digits[point:]
        if item < 0.0:
            token = "-" + token
        return orjson.Fragment(token.encode("ascii"))

    return orjson.dumps({name: convert(data[name], name) for name in _JSON_FIELDS if name in data})


def from_compact_json(
    payload: str | bytes | bytearray | memoryview,
    *,
    Nr: int | None = None,
    Nt: int | None = None,
) -> FitResult:
    """Restore a trusted-local compact record, optionally at another resolution."""

    try:
        data = orjson.loads(payload)
    except orjson.JSONDecodeError as error:
        raise ValueError("invalid compact Equilibrium JSON payload") from error
    if isinstance(data, dict):
        # Descriptions and diagnostics never become authoritative reconstruction inputs.
        data = {name: data[name] for name in _FIELDS if name != "header" and name in data}
    return _restore(data, Nr=Nr, Nt=Nt)


def _restore(data: Any, *, Nr: int | None, Nt: int | None) -> FitResult:
    if isinstance(data, dict):
        data = {"header": None, "K_max": None, **data}
    _keys(data, _FIELDS, "compact Equilibrium")
    if data["header"] is not None and not isinstance(data["header"], str):
        raise TypeError("header must be a string or None")
    for name in ("Nr", "Nt"):
        if type(data[name]) is not int:
            raise TypeError(f"{name} must be an integer")
        if data[name] < 4:
            raise ValueError(f"{name} must be at least 4")
    shape = {name: data[name] for name in (*_BOUNDARY, "Nr", "Nt")}
    for name in ("R0", "Z0", "a", "kappa_lcfs", "B0", "P0"):
        _number(data[name], name)
    for name in ("c_lcfs", "s_lcfs"):
        shape[name] = _array(data[name], name)
    coefficients = {}
    for name in _COEFFICIENTS:
        values = data[name]
        if name in {"c_coeffs", "s_coeffs"}:
            if not isinstance(values, list):
                raise TypeError(f"{name} must be an array of arrays")
            coefficients[name] = [_array(row, name) for row in values]
        else:
            coefficients[name] = [_array(values, name)]
    radial_count = max(1, max(row.size for rows in coefficients.values() for row in rows))
    for name, rows in coefficients.items():
        dense = np.zeros((len(rows), radial_count), dtype=np.float64)
        for index, row in enumerate(rows):
            dense[index, : row.size] = row
        shape[name] = dense if name in {"c_coeffs", "s_coeffs"} else dense[0]
    profiles = {}
    for name in _PROFILES:
        profile = _array(data[name], name)
        if profile.shape != (data["Nr"],):
            raise ValueError(f"{name} must have shape {(data['Nr'],)}, got {profile.shape}")
        profiles[name] = profile
    K_max = data["K_max"]
    if K_max is not None and (type(K_max) is not int or K_max < 2):
        raise ValueError("K_max must be null or an integer >= 2")

    record_nr = data["Nr"]
    record_nt = data["Nt"]
    target_nr = record_nr if Nr is None else Nr
    target_nt = record_nt if Nt is None else Nt
    for name, value in (("Nr", target_nr), ("Nt", target_nt)):
        if type(value) is not int or value < 4:
            raise ValueError(f"{name} must be an integer >= 4")
    if (target_nr, target_nt) != (record_nr, record_nt):
        original_r = _lobatto_quadrature(record_nr)[0]
        target_r = _lobatto_quadrature(target_nr)[0]
        if target_nr != record_nr:
            for name, values in profiles.items():
                profiles[name] = np.interp(target_r, original_r, values)

    psi_r = profiles["psi_r"]
    psi_rr = profiles["psi_rr"]
    accumulator = _spectral_calculus(target_nr)[0]
    weights = _lobatto_quadrature(target_nr)[1]
    # The record cannot represent the original g-file span; the restored state
    # carries the integrated span (== the stored psi_lcfs) as its gauge.
    span = float(np.tensordot(weights, psi_r, axes=(0, 0)))
    psin = accumulator @ psi_r / span
    psin[0], psin[-1] = 0.0, 1.0

    return FitResult(
        Nr=target_nr,
        Nt=target_nt,
        K_max=K_max,
        R0=float(shape["R0"]),
        Z0=float(shape["Z0"]),
        a=float(shape["a"]),
        span=span,
        P0=float(data["P0"]),
        B0=float(data["B0"]),
        kappa_lcfs=float(shape["kappa_lcfs"]),
        c_lcfs=np.asarray(shape["c_lcfs"], dtype=np.float64).copy(),
        s_lcfs=np.asarray(shape["s_lcfs"], dtype=np.float64).copy(),
        h_coeffs=np.asarray(shape["h_coeffs"], dtype=np.float64).copy(),
        v_coeffs=np.asarray(shape["v_coeffs"], dtype=np.float64).copy(),
        kappa_coeffs=np.asarray(shape["kappa_coeffs"], dtype=np.float64).copy(),
        c_coeffs=np.asarray(shape["c_coeffs"], dtype=np.float64).copy(),
        s_coeffs=np.asarray(shape["s_coeffs"], dtype=np.float64).copy(),
        psi_r=psi_r.copy(),
        psi_rr=psi_rr.copy(),
        P_psi=profiles["P_psi"].copy(),
        FF_psi=profiles["FF_psi"].copy(),
        psin=psin,
        detail="restored from compact record",
    )


def _pack(result: FitResult, *, header: str | None) -> dict[str, Any]:
    """Collect the stored roots without touching their precision."""

    if header is not None and not isinstance(header, str):
        raise TypeError("header must be a string or None")
    data: dict[str, Any] = {} if header is None else {"header": header}
    data.update(Nr=result.Nr, Nt=result.Nt)
    data.update({name: getattr(result, name) for name in _PROFILES})
    data.update({name: getattr(result, name) for name in (*_COEFFICIENTS, "K_max", *_BOUNDARY)})
    data.update(B0=result.B0, P0=result.P0)
    for name in _COEFFICIENTS:
        values = data[name]
        # Exact trailing zeros carry no polynomial information. Harmonic
        # rows retain their order; boundary arrays define the harmonic axis.
        data[name] = [_trim(row) for row in values] if values.ndim == 2 else _trim(values)
    if data["K_max"] is None:
        del data["K_max"]
    return data


def _trim(values: np.ndarray) -> np.ndarray:
    nonzero = np.flatnonzero(values.view(np.uint64))
    return values[: int(nonzero[-1]) + 1] if nonzero.size else values[:0]


def _keys(value: Any, expected: tuple[str, ...], name: str) -> None:
    if not isinstance(value, dict):
        raise TypeError(f"{name} must be an object")
    if set(value) != set(expected):
        raise ValueError(
            f"invalid {name} fields; missing={sorted(set(expected) - set(value))}, "
            f"unknown={sorted(set(value) - set(expected))}"
        )


def _number(value: Any, name: str) -> None:
    if type(value) not in {int, float}:
        raise TypeError(f"{name} must be a number")
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")


def _array(value: Any, name: str) -> np.ndarray:
    if isinstance(value, np.ndarray):
        if value.ndim != 1 or value.dtype != np.float64:
            raise TypeError(f"{name} must be a one-dimensional float64 array")
        result = value
    elif isinstance(value, list):
        if any(type(item) not in {int, float} for item in value):
            raise TypeError(f"{name} must contain numbers")
        result = np.asarray(value, dtype=np.float64)
    else:
        raise TypeError(f"{name} must be an array")
    # Validate finiteness once per array, rather than once per JSON number.
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result
