#!/usr/bin/env python3
"""Self-consistency check for compact-JSON records.

A record is internally consistent when its reconstruction roots reproduce its
own display diagnostics through the manuscript's Eqs. (7) and (8):

    P(r)   = P_e  - int_r^1 P_psi(t) psi_r(t) dt          (7)
    F^2(r) = F_e^2 - 2 int_r^1 FF_psi(t) psi_r(t) dt      (8)

with P_e = P0, F_e = R0 B0, evaluated at r = 0 against P_axis and
(R_axis B_axis)^2.  The integral runs over the Nr radial Gauss-Lobatto nodes,
so Gauss-Lobatto quadrature is exact for the polynomial order the grid carries.

This is the check a reader can run without any solver: it exercises the field
semantics (units, sign conventions, ordering) but not the fitting.

    python3 verify_record_consistency.py <archive-root> [--limit N] [--tol T]

Exit status is 0 only when every sampled record passes.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial import legendre as L


def lobatto(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Lobatto-Legendre nodes on [0, 1] and their quadrature weights."""
    coeffs = np.zeros(n)
    coeffs[n - 1] = 1.0
    interior = L.legroots(L.legder(coeffs)).real
    nodes = np.sort(np.r_[-1.0, 1.0, interior])
    weights = 2.0 / (n * (n - 1) * L.legval(nodes, coeffs) ** 2)
    return 0.5 * (nodes + 1.0), 0.5 * weights


def check(path: Path, tol: float) -> list[str]:
    payload = json.loads(path.read_bytes())
    failures: list[str] = []

    nr = payload["Nr"]
    if nr < 4:
        return ["Nr=%r is below the documented minimum of four" % nr]
    r, w = lobatto(nr)

    for name in ("P_psi", "FF_psi", "psi_r"):
        if len(payload[name]) != nr:
            failures.append("%s has %d samples, expected Nr=%d"
                            % (name, len(payload[name]), nr))
    if failures:
        return failures

    p_psi = np.asarray(payload["P_psi"], dtype=float)
    ff_psi = np.asarray(payload["FF_psi"], dtype=float)
    psi_r = np.asarray(payload["psi_r"], dtype=float)

    # Eq. (7) at r = 0: the first node is r = 0 for Gauss-Lobatto.
    p0 = payload["P0"] - float(np.dot(w, p_psi * psi_r))
    f_edge = payload["R0"] * payload["B0"]
    f0_sq = f_edge ** 2 - 2.0 * float(np.dot(w, ff_psi * psi_r))

    p_axis = payload["P_axis"]
    f_axis_sq = (payload["R_axis"] * payload["B_axis"]) ** 2

    p_scale = max(abs(p_axis), 1.0)
    if abs(p0 - p_axis) > tol * p_scale:
        failures.append("Eq.(7) P(0)=%.10g vs P_axis=%.10g (rel %.3g)"
                        % (p0, p_axis, abs(p0 - p_axis) / p_scale))
    f_scale = max(abs(f_axis_sq), 1.0)
    if abs(f0_sq - f_axis_sq) > tol * f_scale:
        failures.append("Eq.(8) F(0)^2=%.10g vs (R_axis B_axis)^2=%.10g (rel %.3g)"
                        % (f0_sq, f_axis_sq, abs(f0_sq - f_axis_sq) / f_scale))
    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--limit", type=int, default=0,
                        help="check at most this many records (0 = all)")
    parser.add_argument("--tol", type=float, default=1e-6,
                        help="relative tolerance on each identity")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    paths = sorted(args.archive.rglob("*.compact.json"))
    if args.limit:
        step = max(1, len(paths) // args.limit)
        paths = paths[::step][:args.limit]

    bad = 0
    worst = (0.0, "")
    for path in paths:
        failures = check(path, args.tol)
        if failures:
            bad += 1
            if not args.quiet:
                print("FAIL %s" % path.name)
                for failure in failures:
                    print("     %s" % failure)
        elif not args.quiet and len(paths) <= 20:
            print("OK   %s" % path.name)

    print("checked %d records, %d failed, tolerance %g" % (len(paths), bad, args.tol))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
