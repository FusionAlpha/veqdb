# v0.1.1 — 5 November 2026

Publish the Figure 3 pressure verification material accompanying the
Scientific Data submission.

- Archive the EAST seed derived from `g048059.03650` and 15 full-precision
  compact states from the 5 November verification replay.
- Provide the state table, solver and input provenance, source checksums and
  an explicitly labelled excerpt of the historical scan summary.
- Supply `scripts/verify_figure3_pressure.py` to check record diagnostics,
  independently integrate the pressure gradient and compare magnetic-axis
  coordinates with the historical scan.
- Verify the manuscript range of 17.4–23.4 kPa. The replay range is
  17.445862–23.398263 kPa with the archived seed edge-pressure input.

The original summary's approximately doubled pressure diagnostic is retained
as `historical_P_axis_Pa` with status `unreconciled_legacy_diagnostic`.
Tracing its original writing step requires the summary-generation script.
Re-running the fixed-boundary solves requires the historical solver environment
identified in the provenance. The public verification command runs with the
VEQDB package dependencies.

The primary release comprises 13,291 sampled equilibria and three supplementary
comparison records. The Figure 3 replay states are a separate validation bundle.

# v0.1.0 — 4 November 2026

First tagged GitHub release of the VEQDB standalone compact-record codec and
dataset repository. The package version remains 0.1.0.

- Decode and re-encode compact equilibrium records with `veqdb.from_compact_json`
  and `veqdb.to_compact_json`, without importing VEQPy or Fusionprime.
- Retain exact reconstruction inputs and rounded golden encodings while allowing
  floating-point roundoff in full-precision diagnostic reference comparisons
  (`rtol=1e-12`, `atol=1e-18`). These are test tolerances.
- Document the Figshare DOI, data and codec licence scopes, and the pinned
  source repositories for the EAST, CHEASE and EFIT manuscript inputs.
- Match the figure script's panel spacing and shared-axis presentation to
  the manuscript.

The original 13,291-record archive and three supplementary reference records
are unchanged. The codec implementation is unchanged from commit `4ef67f3`.

The source archive includes the package, tests, data archives, checksum
manifests and source documentation. Python 3.12 or newer is required; see
the README for installation and verification commands.

Full manuscript figure regeneration additionally requires the frozen inputs
and VEQPy package context described in the README. The optional direct VEQPy
source-parity test requires its separate source environment.

Dataset DOI: https://doi.org/10.6084/m9.figshare.34485114
