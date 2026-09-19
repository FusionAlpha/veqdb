# G-EQDSK Projection Snapshot

This directory preserves the numerical projection implementation used to
produce the G-EQDSK-derived manuscript records. It fits a G-EQDSK flux map and
LCFS jointly to the MXH--Chebyshev representation, then applies rank,
geometry, flux monotonicity and physical-reconstruction acceptance checks.

The implementation is intentionally kept as source provenance rather than
advertised as a second public API. Its model-layer dependencies are being moved
to VEQPy, which is the solver and G-EQDSK I/O boundary for VEQDB. The frozen
reference inputs and accepted compact records in this repository permit
inspection and replotting without re-running this fitting path.
