#!/usr/bin/env python3
"""Verify same-case PNG/JSON coverage and hashes without extracting archives.

Run from any directory: python3 scripts/verify_png_records.py
This checks file associations, not image rendering or equilibrium topology.
"""
from __future__ import annotations

import hashlib
import json
import sys
import tarfile
from contextlib import ExitStack
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = "manifests/VEQDB-png-partial-20260920.records.json"


def verify(root: Path) -> tuple[int, int, int]:
    manifest = json.loads((root / MANIFEST).read_bytes())
    if manifest["schema"] != "veqdb-png-records-v1":
        raise ValueError("unsupported manifest schema")
    entries = manifest["entries"]
    if not entries:
        raise ValueError("manifest contains no PNG entries")

    with ExitStack() as stack:
        archives: dict[str, tarfile.TarFile] = {}

        def archive(name: str) -> tarfile.TarFile:
            if name not in archives:
                archives[name] = stack.enter_context(tarfile.open(root / name))
            return archives[name]

        def member_bytes(name: str, path: str) -> bytes:
            member = archive(name).getmember(path)
            if not member.isfile():
                raise ValueError(f"not a regular archive member: {path}")
            stream = archive(name).extractfile(member)
            if stream is None:
                raise ValueError(f"cannot read archive member: {path}")
            with stream:
                return stream.read()

        png_archive = manifest["png_archive"]
        png_paths = [m.name for m in archive(png_archive).getmembers()
                     if m.isfile() and m.name.endswith(".png")]
        listed = [entry["png_path"] for entry in entries]
        if len(set(png_paths)) != len(png_paths):
            raise ValueError("PNG archive contains duplicate paths")
        if len(set(listed)) != len(listed):
            raise ValueError("manifest contains duplicate PNG paths")
        if set(listed) != set(png_paths):
            raise ValueError("manifest PNG coverage differs from archive: "
                             f"{sorted(set(listed) ^ set(png_paths))}")

        archived = standalone = 0
        for entry in entries:
            png_path, record_path = entry["png_path"], entry["record_path"]
            if (not record_path.endswith(".compact.json") or
                    Path(png_path).stem !=
                    Path(record_path).name.removesuffix(".compact.json")):
                raise ValueError(f"case identifiers differ: {png_path}, {record_path}")
            png = member_bytes(png_archive, png_path)
            if entry["record_archive"] is None:
                record = (root / record_path).read_bytes()
                standalone += 1
            else:
                record = member_bytes(entry["record_archive"], record_path)
                archived += 1
            for label, data in (("png", png), ("record", record)):
                if hashlib.sha256(data).hexdigest() != entry[f"{label}_sha256"]:
                    raise ValueError(f"{label} SHA-256 mismatch for {png_path}")
            payload = json.loads(record)
            if not isinstance(payload, dict) or not {
                "c_lcfs", "s_lcfs", "P_psi", "FF_psi", "psi_r", "psi_rr"
            }.issubset(payload):
                raise ValueError(f"not a compact equilibrium record: {record_path}")
    return len(entries), archived, standalone


def main() -> int:
    try:
        total, archived, standalone = verify(ROOT)
    except (OSError, ValueError, KeyError, TypeError, tarfile.TarError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1
    print(f"PASS: {total} PNGs paired; {archived} archive records, "
          f"{standalone} supplementary records; all hashes match")
    return 0


if __name__ == "__main__":
    sys.exit(main())
