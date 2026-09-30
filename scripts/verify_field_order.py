#!/usr/bin/env python3
"""Check the field order and field count of compact-JSON records.

The expected order is the one documented in FORMAT.md, which restates Table A.7
of the manuscript. Only `K_max` is optional, and it is always adjacent to `Nt`.

    python3 verify_field_order.py <archive-root>

Exit status is 0 only when every record matches.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

# FORMAT.md section 2, in the writer's output order.
EXPECTED = (
    "header",
    "R_axis", "Z_axis", "B_axis", "P_axis",
    "B0", "P0", "R0", "Z0", "a",
    "betat", "Ip", "qmin", "q0", "q95", "Wth", "psi_lcfs",
    "kappa_lcfs", "c_lcfs", "s_lcfs",
    "Nr", "Nt", "K_max",
    "h_coeffs", "v_coeffs", "kappa_coeffs", "c_coeffs", "s_coeffs",
    "P_psi", "FF_psi", "psi_r", "psi_rr",
)

WITHOUT_KMAX = tuple(name for name in EXPECTED if name != "K_max")
ALLOWED = {EXPECTED: 32, WITHOUT_KMAX: 31}

# Roles, as documented. Used only for the summary line.
RECONSTRUCTION = 19   # 17 marked Reconstruction + Nr + Nt
DIAGNOSTIC = 11
DESCRIPTION = 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()

    orders: Counter = Counter()
    bad: list[tuple[str, list[str]]] = []
    total = 0
    for path in sorted(args.archive.rglob("*.compact.json")):
        total += 1
        keys = tuple(json.loads(path.read_bytes()).keys())
        orders[keys] += 1
        if keys not in ALLOWED:
            bad.append((str(path), list(keys)))

    print("records checked : %d" % total)
    for keys, count in orders.most_common():
        print("  %7d x %2d fields" % (count, len(keys)))

    if total and len(orders) == 1 and not bad:
        keys = next(iter(orders))
        expected_roles = RECONSTRUCTION + DIAGNOSTIC + DESCRIPTION
        if len(keys) == expected_roles:
            print("OK  single key order; %d fields (%d reconstruction + %d "
                  "diagnostic + %d description)"
                  % (len(keys), RECONSTRUCTION, DIAGNOSTIC, DESCRIPTION))
        else:
            print("OK  single key order; %d fields" % len(keys))

    if bad:
        print("\nFAIL  %d record(s) deviate from the documented order" % len(bad))
        for name, keys in bad[:20]:
            print("  %s" % name)
            print("    %s" % keys)
        return 1

    print("\nPASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
