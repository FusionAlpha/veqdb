#!/usr/bin/env python3
"""Verify the published Figure 3 replay records using the standalone codec.

Run from an installed source checkout with ``python scripts/verify_figure3_pressure.py``.
The JSON report is written to standard output. Exit status is zero on success.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import BarycentricInterpolator
from verify_record_consistency import check, lobatto

from veqdb import from_compact_json, to_compact_json


def verify(directory: Path) -> dict:
    directory = directory.resolve()
    provenance = json.loads((directory / "provenance.json").read_text())
    for name, expected in provenance["file_sha256"].items():
        path = (directory / name).resolve()
        if not path.is_relative_to(directory):
            raise ValueError(f"File path leaves the validation directory: {name}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch for {name}")

    seed_path = directory / "seed.compact.json"
    seed = from_compact_json(seed_path.read_bytes())
    seed_failures = check(seed_path, 1e-6)
    if seed_failures:
        raise ValueError(f"Seed integral consistency: {seed_failures}")
    np.testing.assert_allclose(seed.P0, provenance["input_P0_Pa"], rtol=0, atol=0)

    with (directory / "states.csv").open(newline="") as handle:
        states = list(csv.DictReader(handle))
    if len(states) != 15 or len({row["state_id"] for row in states}) != 15:
        raise ValueError("Expected 15 distinct replay states")
    expected_pairs = {(kappa, shift) for kappa in (-0.15, -0.08, 0, 0.08, 0.15)
                      for shift in (-0.04, 0, 0.04)}
    actual_pairs = {(float(row["dkappa"]), float(row["dR0_m"])) for row in states}
    if actual_pairs != expected_pairs:
        raise ValueError("The replay states must cover the documented 5 by 3 scan")

    history = json.loads((directory / "historical-summary.json").read_text())
    references = {row["state_id"]: row for row in history["states"]}
    quad_nodes, quad_weights = np.polynomial.legendre.leggauss(128)
    quad_nodes = 0.5 * (quad_nodes + 1.0)
    quad_weights = 0.5 * quad_weights
    checks = []

    for row in states:
        state_id = row["state_id"]
        record_name = row["record_file"]
        if record_name not in provenance["file_sha256"]:
            raise ValueError(f"Record has no published checksum: {record_name}")
        path = directory / record_name
        payload = path.read_bytes()
        stored = json.loads(payload)
        restored = from_compact_json(payload)
        recomputed = json.loads(to_compact_json(restored, precision=None))
        failures = check(path, 1e-6)
        if failures:
            raise ValueError(f"{state_id} integral consistency: {failures}")

        radial_nodes, _ = lobatto(restored.Nr)
        pressure_derivative = BarycentricInterpolator(radial_nodes, restored.P_psi)(quad_nodes)
        flux_derivative = BarycentricInterpolator(radial_nodes, restored.psi_r)(quad_nodes)
        independent_pressure = restored.P0 - float(
            np.dot(quad_weights, pressure_derivative * flux_derivative)
        )
        reference_pressure = float(row["P_axis_Pa"])
        np.testing.assert_allclose(
            [stored["P_axis"], recomputed["P_axis"], independent_pressure],
            reference_pressure, rtol=0, atol=1e-5, err_msg=state_id,
        )
        np.testing.assert_allclose(restored.P0, float(row["P0_output_Pa"]), rtol=0, atol=1e-8)
        np.testing.assert_allclose(float(row["P0_input_Pa"]), seed.P0, rtol=0, atol=0)
        np.testing.assert_allclose(restored.R0, seed.R0 + float(row["dR0_m"]), rtol=0, atol=1e-10)
        np.testing.assert_allclose(
            restored.kappa_lcfs, seed.kappa_lcfs * (1 + float(row["dkappa"])),
            rtol=0, atol=1e-10,
        )
        reference = references[state_id]
        np.testing.assert_allclose(
            [recomputed["R_axis"], recomputed["Z_axis"]],
            [float(row["R_axis_m"]), float(row["Z_axis_m"])],
            rtol=0, atol=1e-9, err_msg=state_id,
        )
        np.testing.assert_allclose(
            [recomputed["R_axis"], recomputed["Z_axis"]],
            [reference["R_axis_m"], reference["Z_axis_m"]],
            rtol=0, atol=1e-8, err_msg=state_id,
        )
        np.testing.assert_allclose(recomputed["Ip"], reference["Ip_A"], rtol=0, atol=1e-3)
        checks.append({
            "state_id": state_id,
            "P_axis_Pa": recomputed["P_axis"],
            "independent_pressure_Pa": independent_pressure,
            "independent_pressure_difference_Pa": abs(independent_pressure - recomputed["P_axis"]),
            "R_axis_difference_m": abs(recomputed["R_axis"] - reference["R_axis_m"]),
            "Z_axis_difference_m": abs(recomputed["Z_axis"] - reference["Z_axis_m"]),
        })

    pressure_range = [min(row["P_axis_Pa"] for row in checks),
                      max(row["P_axis_Pa"] for row in checks)]
    rounded_range = [round(value / 1000, 1) for value in pressure_range]
    if rounded_range != [17.4, 23.4]:
        raise ValueError(f"Unexpected rounded pressure range: {rounded_range}")
    return {
        "passed": True,
        "replay_date": provenance["replay_date"],
        "replay_states": len(checks),
        "seed_checked": True,
        "P_axis_range_Pa": pressure_range,
        "rounded_P_axis_range_kPa": rounded_range,
        "max_independent_pressure_difference_Pa": max(
            row["independent_pressure_difference_Pa"] for row in checks
        ),
        "max_historical_axis_difference_m": max(
            max(row["R_axis_difference_m"], row["Z_axis_difference_m"]) for row in checks
        ),
        "historical_pressure_status": history["historical_pressure_status"],
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory", nargs="?", type=Path,
        default=Path(__file__).resolve().parents[1] / "data/validation/figure3",
    )
    args = parser.parse_args()
    try:
        report = verify(args.directory)
    except (AssertionError, KeyError, OSError, ValueError) as error:
        print(f"Figure 3 verification failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
