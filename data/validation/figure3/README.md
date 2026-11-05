# Figure 3 pressure verification

This bundle supports the 17.4–23.4 kPa on-axis pressure range reported for the
15-state EAST scan in the manuscript. The compact records were generated in
a verification replay on 5 November 2026. The primary database contains
13,291 sampled equilibria and three supplementary comparison records.

## Files

- `seed.compact.json` contains the archived EAST seed derived from
  `g048059.03650`, discharge 48059 at 3.65 s.
- `records/` contains 15 full-precision compact records from the verification
  replay, retaining the seed input edge pressure of 29.6742825 Pa.
- `states.csv` maps scan coordinates, input and output edge pressure, axis
  pressure, magnetic-axis coordinates, plasma current and solver residual to
  each replay record. Its column names specify SI units.
- `historical-summary.json` contains the corresponding 15 rows selected from
  the original 45-row scan summary, plus source identifiers and the original
  summary checksum. The original scan repeated each geometric state at three
  common source-scaling levels.
- `provenance.json` records the seed and solver revisions, generation settings,
  environment and file checksums.
- `verification.json` records the output of the public verification command.

## Run the check

From the repository root, install the package with Python 3.12 or later and run
the standalone check. Dependencies are declared in `pyproject.toml`.

```sh
python -m pip install .
python scripts/verify_figure3_pressure.py
```

The check verifies checksums and the 5-by-3 scan, decodes all records through
the public VEQDB codec, recomputes their diagnostics, and independently
integrates the reconstructed pressure gradient with 128-point Gauss–Legendre
quadrature. It also compares the magnetic-axis coordinates and plasma current
with the historical summary. A successful check exits with status zero and
prints a JSON report.

The pressure definition is

```text
P_axis = P0 - integral_0^1 P_psi(r) * psi_r(r) dr
```

The replay gives 17.445862–23.398263 kPa, which rounds to the manuscript range.
`P0_input_Pa` records the prescribed edge pressure. `P0_output_Pa` records its
value after the solver's common pressure normalization at fixed plasma current.

## Historical pressure diagnostic

The original summary contains values of 34.9–46.8 kPa in its `P_axis_Pa`
field. Those values are retained here as `historical_P_axis_Pa` and marked
`unreconciled_legacy_diagnostic`. Their original writing step remains to be
traced. The manuscript pressure check uses the replay records and the integral
above. Magnetic-axis coordinates from the replay agree with the historical
summary within 6.4e-11 m.

The archived solver revision and complete replay settings are identified in
`provenance.json`. Re-running the fixed-boundary solves requires that solver
environment. The command above uses the publicly distributed VEQDB codec,
NumPy and SciPy to verify the stored replay states.

The EAST source repository, fixed revision and file path are listed in
[`SOURCE_DATA.md`](../../../SOURCE_DATA.md). Applicable data and source terms
are recorded in [`LICENSE`](../../../LICENSE) and
[`THIRD_PARTY_NOTICES.md`](../../../THIRD_PARTY_NOTICES.md).
