# Compact JSON field reference

Every equilibrium distributed in this repository is one flat JSON object. This
document describes that format: the object shape, the meaning and units of each
field, the encoding rules, and what a reader has to do to use a record.

It restates Appendix A (Table A.7) of the accompanying manuscript, which is the
normative reference for the format. The statements below are checked against the
distributed records by the scripts in `scripts/`; see §5.

## 1. Object

- One flat object per record. The only nested structure is inside the
  coefficient and profile arrays.
- **31 fields** by default. `K_max` is conditional: when it is present the
  object has **32** fields.
- The field order in §2 is the writer's output order and is normative. All
  13,291 records here use it, with no exceptions.
- Records are stored as
  `<device>/<device>__GAQ_<id>_<index>.compact.json`, where `<index>` is
  zero-padded to four digits and unique within a device.

## 2. Fields

The **Role** column is what a reader may do with the field:

| Role | Count | Meaning |
| --- | --- | --- |
| Reconstruction | 17 | Required input to reconstruct the equilibrium. |
| Configuration | 2 | Discrete grid parameters; required input. |
| Diagnostic | 11 | Derived from the reconstruction roots. Recompute; do not read back. |
| Description | 1 | Free-form label. Ignored by a reader. |
| Conditional | 1 | `K_max`; present only when it changes a represented harmonic. |

Units are SI.

| # | Field | Role | Meaning |
| --- | --- | --- | --- |
| 1 | `header` | Description | Optional human-readable string; the public default is `"Equilibrium"`. Not a reconstruction input. |
| 2 | `R_axis` | Diagnostic | Magnetic-axis major-radius coordinate, in m. |
| 3 | `Z_axis` | Diagnostic | Magnetic-axis vertical coordinate, in m. |
| 4 | `B_axis` | Diagnostic | Signed toroidal magnetic field at the magnetic axis, in T. |
| 5 | `P_axis` | Diagnostic | Pressure at the magnetic axis, in Pa. |
| 6 | `B0` | Reconstruction | Signed reference toroidal field, in T; anchors the LCFS field function through `F(1) = R0 * B0`. |
| 7 | `P0` | Reconstruction | LCFS pressure `P(1)`, in Pa. |
| 8 | `R0` | Reconstruction | Major-radius coordinate of the boundary reference frame, in m. |
| 9 | `Z0` | Reconstruction | Vertical coordinate of the boundary reference frame, in m. |
| 10 | `a` | Reconstruction | Minor-radius scale of the boundary reference frame, in m. |
| 11 | `betat` | Diagnostic | Dimensionless toroidal beta, `2 * mu0 * <P> / B0^2`. |
| 12 | `Ip` | Diagnostic | Total plasma current, in A. The G-EQDSK comparison records use the explicit positive-current convention. |
| 13 | `qmin` | Diagnostic | Signed minimum of the safety-factor values on the stored radial grid. |
| 14 | `q0` | Diagnostic | Signed safety factor at the magnetic-axis node, `q[0]`. |
| 15 | `q95` | Diagnostic | Safety factor linearly interpolated at normalized poloidal flux `psi_n = 0.95`. |
| 16 | `Wth` | Diagnostic | Thermal energy `(3/2) * integral(P dV)`, in J. |
| 17 | `psi_lcfs` | Diagnostic | LCFS poloidal flux `psi[-1]` in the axis-zero gauge, in Wb/rad. The comparison records use an outward-increasing flux convention. |
| 18 | `kappa_lcfs` | Reconstruction | Dimensionless LCFS elongation `kappa(1)`. |
| 19 | `c_lcfs` | Reconstruction | Dimensionless LCFS cosine phase coefficients `[c0, c1, ...]`. |
| 20 | `s_lcfs` | Reconstruction | Dimensionless LCFS sine phase coefficients `[s1, s2, ...]`. |
| 21 | `Nr` | Configuration | Number of radial Gauss-Lobatto nodes; an integer of at least four. |
| 22 | `Nt` | Configuration | Number of uniformly spaced periodic poloidal nodes; an integer of at least four and sufficient to resolve the retained harmonics. |
| 23 | `K_max` | Conditional | Optional cap `K_m = min(m, K_max)`. Omitted when unset or when `K_max >= M_max`, because the cap then changes no represented harmonic. Inserted immediately after `Nt`. |
| 24 | `h_coeffs` | Reconstruction | One-dimensional Chebyshev coefficients for the interior horizontal-shift profile `h(r)`. |
| 25 | `v_coeffs` | Reconstruction | One-dimensional Chebyshev coefficients for the interior vertical-shift profile `v(r)`. |
| 26 | `kappa_coeffs` | Reconstruction | One-dimensional Chebyshev coefficients for the interior correction to `kappa(r)` about `kappa_lcfs`. |
| 27 | `c_coeffs` | Reconstruction | Array of Chebyshev-coefficient rows for the interior cosine phase profiles, ordered `c0, c1, ...`. |
| 28 | `s_coeffs` | Reconstruction | Array of Chebyshev-coefficient rows for the interior sine phase profiles, ordered `s1, s2, ...`. |
| 29 | `P_psi` | Reconstruction | Samples of `P_psi = dP/dpsi` on the `Nr` radial nodes, in `A/m^3` = `Pa/(Wb/rad)`. |
| 30 | `FF_psi` | Reconstruction | Samples of `FF_psi = d(F^2/2)/dpsi` on the `Nr` radial nodes, in T. |
| 31 | `psi_r` | Reconstruction | Samples of `psi_r = dpsi/dr` on the `Nr` radial nodes, in Wb/rad. |
| 32 | `psi_rr` | Reconstruction | Independent samples of `psi_rr = d^2psi/dr^2` on the `Nr` radial nodes, in Wb/rad. |

In the 13,291 records distributed here, rows 1-22 and 24-32 are present;
row 23 is absent. In the general format, row 1 (`header`) may also be omitted
when the writer receives `header=None`; `K_max` (row 23) is conditional.
Numbering above is for reference only and is not part of the format.

## 3. Encoding rules

- **Precision.** The writer's default is `precision=10`, at most ten significant
  digits. `precision=None` preserves the source `float64` values in decimal form.
  The records in this repository were written with twelve significant digits.
- Both modes require finite floating-point values and preserve signed zero.
- Exact trailing positive zeros are omitted from coefficient vectors.
- The cosine and sine coefficient families additionally omit trailing
  all-positive-zero harmonic rows. An all-zero row between retained harmonics is
  written as `[0.0]`, so that its index is not lost.
- On input the reader constructs `float64` arrays and restores omitted trailing
  entries by zero-padding to the common radial coefficient width.
- `c_lcfs` and `c_coeffs` use cosine orders `c0, c1, ...`; `s_lcfs` and
  `s_coeffs` use sine orders `s1, s2, ...`.
- Passing `header=None` omits only the description.
- `Nr` and `Nt` are JSON integers. In the records distributed here, zero-valued
  `P0` and `Z0` are written as `0`, not `0.0`. Other integral-valued floating-point
  fields can use exponent notation (for example, a large `Wth` or `P_psi`).
  A reader must not depend on the JSON number type or textual notation of a
  floating-point field.

## 4. Using a record

A reader needs the 19 reconstruction fields: the 17 rows marked
`Reconstruction`, plus `Nr` and `Nt`. When `K_max` is present it is a twentieth
input.

The 11 `Diagnostic` fields are derived from those roots and are not
authoritative. A reader recomputes them; the stored values exist so that a
record can be inspected without a solver. This is what lets diagnostic
definitions change in future without changing the meaning of an archived
equilibrium.

Missing diagnostics therefore do not prevent reconstruction, whereas every
applicable reconstruction field is required. Unrecognized extra fields are
ignored.

Two identities tie the roots to the diagnostics and are the cheapest way to
confirm that a record is internally consistent:

```
P(r)   = P0      - integral_r^1 P_psi(t) * psi_r(t) dt      (Eq. 7)
F^2(r) = (R0*B0)^2 - 2 * integral_r^1 FF_psi(t) * psi_r(t) dt   (Eq. 8)
```

Evaluated at `r = 0`, the first must reproduce `P_axis` and the second must
reproduce `(R_axis * B_axis)^2`.

## 5. Verifying a record

The downloads are covered by `SHA256SUMS`:

```sh
sha256sum -c SHA256SUMS
```

Individual records are covered by `manifests/VEQDB-equilibria-20260920.sha256`,
which carries one SHA-256 per record with the archive's top-level directory as
the path prefix:

```sh
tar xzf VEQDB-equilibria-20260920.tar.gz
sha256sum -c manifests/VEQDB-equilibria-20260920.sha256
```

The identities of §4 are checked by `scripts/verify_record_consistency.py`, which
needs only the records and NumPy:

```sh
python3 scripts/verify_record_consistency.py /path/to/VEQDB-equilibria-20260920
```

`scripts/verify_field_order.py` checks the field order and field count of §2
against the records:

```sh
python3 scripts/verify_field_order.py /path/to/VEQDB-equilibria-20260920
```

## 6. Scope of this document

This document defines the format. It does not ship a reader: no implementation
of the reconstruction step is distributed in this repository. The identities in
§4 are a consistency check on the stored fields, not a reconstruction — they
verify units, sign conventions and ordering, not the fitting that produced the
coefficients.
