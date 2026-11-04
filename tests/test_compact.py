from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import veqdb
from veqdb.compact import from_compact_json, to_compact_json

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "compact"
DATA = ROOT / "data" / "figure-data"
SOURCE_REVISION = "586afb717099fcc1eb17ac5ffca4a10a95ae6d80"
SOURCE_FILES = (
    "veqpy/fitting/compact.py",
    "veqpy/fitting/diagnostics.py",
    "veqpy/fitting/state.py",
    "veqpy/fitting/_spectral.py",
    "veqpy/fitting/_geometry.py",
)


def fixture_bytes(name: str) -> bytes:
    return (FIXTURES / name).read_bytes()


def test_public_exports() -> None:
    assert veqdb.from_compact_json is from_compact_json
    assert veqdb.to_compact_json is to_compact_json
    assert isinstance(from_compact_json(fixture_bytes("east-precision-full.compact.json")), veqdb.FitResult)


@pytest.mark.parametrize("case", ["east", "mast-u", "iter"])
def test_precision_10_matches_source_golden(case: str) -> None:
    source = (DATA / f"crossdev_{case}_p100.compact.json").read_bytes()
    actual = to_compact_json(from_compact_json(source), precision=10)
    assert actual == fixture_bytes(f"{case}-precision-10.compact.json")


@pytest.mark.parametrize("precision", [None, 12])
def test_east_precisions_match_source_golden(precision: int | None) -> None:
    source = (DATA / "crossdev_east_p100.compact.json").read_bytes()
    actual = to_compact_json(from_compact_json(source), precision=precision)
    label = "full" if precision is None else str(precision)
    assert actual == fixture_bytes(f"east-precision-{label}.compact.json")


def test_full_precision_round_trip_is_byte_identical() -> None:
    payload = fixture_bytes("east-precision-full.compact.json")
    assert to_compact_json(from_compact_json(payload), precision=None) == payload


def test_resolution_override_matches_source_golden() -> None:
    source = (DATA / "crossdev_east_p100.compact.json").read_bytes()
    restored = from_compact_json(source, Nr=40, Nt=24)
    assert restored.Nr == 40
    assert restored.Nt == 24
    assert to_compact_json(restored, precision=10) == fixture_bytes(
        "east-nr40-nt24-precision-10.compact.json"
    )


def test_decode_restores_raw_roots_and_integrated_flux() -> None:
    payload = fixture_bytes("east-precision-full.compact.json")
    record = json.loads(payload)
    result = from_compact_json(payload)
    for name in ("P_psi", "FF_psi", "psi_r", "psi_rr", "psin"):
        values = getattr(result, name)
        assert values.shape == (record["Nr"],)
        assert values.dtype == np.float64
    np.testing.assert_array_equal(result.P_psi, record["P_psi"])
    np.testing.assert_array_equal(result.FF_psi, record["FF_psi"])
    np.testing.assert_array_equal(result.psi_r, record["psi_r"])
    np.testing.assert_array_equal(result.psi_rr, record["psi_rr"])
    assert result.span == pytest.approx(record["psi_lcfs"])
    assert result.psin[0] == 0.0
    assert result.psin[-1] == 1.0
    assert result.detail == "restored from compact record"


def test_stored_diagnostics_and_unknown_fields_are_not_authoritative() -> None:
    payload = fixture_bytes("east-precision-full.compact.json")
    record = json.loads(payload)
    for name in ("R_axis", "B_axis", "betat", "Ip", "qmin", "q0", "q95", "Wth"):
        record[name] = 12345.0
    record["unknown_future_field"] = {"ignored": True}
    restored = from_compact_json(json.dumps(record))
    assert to_compact_json(restored, precision=None) == payload


def test_signed_zero_and_harmonic_trimming() -> None:
    result = from_compact_json(fixture_bytes("east-precision-full.compact.json"))
    h_coeffs = np.concatenate((result.h_coeffs, np.array([-0.0, 0.0])))
    c_coeffs = np.vstack(
        (
            result.c_coeffs[:-1],
            np.zeros((1, result.c_coeffs.shape[1])),
            result.c_coeffs[-1:],
            np.zeros((2, result.c_coeffs.shape[1])),
        )
    )
    modified = replace(result, h_coeffs=h_coeffs, c_coeffs=c_coeffs)
    payload = to_compact_json(modified, precision=None)
    record = json.loads(payload)
    assert math.copysign(1.0, record["h_coeffs"][-1]) == -1.0
    assert record["h_coeffs"][-1] == 0.0
    assert record["c_coeffs"][-2] == [0.0]
    assert len(record["c_coeffs"]) == result.c_coeffs.shape[0] + 1


def test_k_max_omission_and_negative_b0() -> None:
    result = from_compact_json(fixture_bytes("east-precision-full.compact.json"))
    redundant = json.loads(to_compact_json(replace(result, K_max=10_000)))
    assert "K_max" not in redundant
    active = json.loads(to_compact_json(replace(result, K_max=2)))
    assert active["K_max"] == 2
    negative = json.loads(to_compact_json(replace(result, B0=-result.B0), precision=None))
    assert negative["B0"] == -result.B0
    assert math.copysign(1.0, negative["B_axis"]) == math.copysign(1.0, negative["B0"])


@pytest.mark.parametrize("precision,error", [(0, ValueError), (-1, ValueError), (1.5, TypeError), (True, TypeError)])
def test_invalid_precision(precision: object, error: type[Exception]) -> None:
    result = from_compact_json(fixture_bytes("east-precision-full.compact.json"))
    with pytest.raises(error):
        to_compact_json(result, precision=precision)  # type: ignore[arg-type]


def test_invalid_json_and_top_level_type() -> None:
    with pytest.raises(ValueError, match="invalid compact Equilibrium JSON"):
        from_compact_json(b"{")
    with pytest.raises(TypeError, match="must be an object"):
        from_compact_json(b"[]")


@pytest.mark.parametrize(
    "key,value,error",
    [
        ("Nr", "32", TypeError),
        ("Nr", 3, ValueError),
        ("K_max", 1, ValueError),
        ("P_psi", [1.0, 2.0], ValueError),
        ("R0", "1.0", TypeError),
    ],
)
def test_invalid_record_fields(key: str, value: object, error: type[Exception]) -> None:
    record = json.loads(fixture_bytes("east-precision-full.compact.json"))
    record[key] = value
    with pytest.raises(error):
        from_compact_json(json.dumps(record))


def test_missing_reconstruction_field_is_rejected() -> None:
    record = json.loads(fixture_bytes("east-precision-full.compact.json"))
    del record["psi_r"]
    with pytest.raises(ValueError, match="missing=.*psi_r"):
        from_compact_json(json.dumps(record))


@pytest.mark.parametrize("kwargs", [{"Nr": 3}, {"Nt": 2.0}])
def test_invalid_resolution_override(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValueError, match="must be an integer >= 4"):
        from_compact_json(fixture_bytes("east-precision-full.compact.json"), **kwargs)  # type: ignore[arg-type]


@pytest.mark.skipif(not os.environ.get("VEQPY_PARITY_ROOT"), reason="set VEQPY_PARITY_ROOT for source parity")
def test_explicit_source_parity() -> None:
    source_root = os.environ["VEQPY_PARITY_ROOT"]
    checked = subprocess.run(
        ["git", "-C", source_root, "diff", "--quiet", SOURCE_REVISION, "--", *SOURCE_FILES],
        check=False,
    )
    if checked.returncode:
        pytest.fail("VEQPY_PARITY_ROOT relevant files differ from the pinned source revision")
    sys.path.insert(0, source_root)
    from veqpy.fitting.compact import from_compact_json as source_decode
    from veqpy.fitting.compact import to_compact_json as source_encode

    for case in ("east", "mast-u", "iter"):
        payload = (DATA / f"crossdev_{case}_p100.compact.json").read_bytes()
        for precision in (None, 10, 12):
            expected = source_encode(source_decode(payload), precision=precision)
            actual = to_compact_json(from_compact_json(payload), precision=precision)
            assert actual == expected
    payload = (DATA / "crossdev_east_p100.compact.json").read_bytes()
    expected = source_encode(source_decode(payload, Nr=40, Nt=24), precision=10)
    actual = to_compact_json(from_compact_json(payload, Nr=40, Nt=24), precision=10)
    assert actual == expected
