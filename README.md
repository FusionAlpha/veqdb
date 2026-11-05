# VEQDB

VEQDB is the release repository for a compact database of fixed-boundary
tokamak equilibria represented with MXH--Chebyshev geometry and physical-profile
roots. It accompanies the manuscript *A compact and reconstructible multi-device
tokamak equilibrium database*. An earlier version is available as
[arXiv:2609.23296](https://arxiv.org/abs/2609.23296).

The initial release contains 13,291 accepted equilibria from 267 conventional
and spherical tokamak devices. The compact-JSON record set occupies 41.24 MB
before archive compression, and 17.47 MB in the distributed archive.

Three additional G-EQDSK-derived reference records are distributed separately
in `data/figure-data/`. Together with the 13,291 GAQ parameter-sampled records
in the original archive, they provide 13,294 records. The original archive and
its per-record checksum manifest are unchanged.

The dataset is archived in [Figshare](https://doi.org/10.6084/m9.figshare.34485114)
under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Cite this DOI
when using the dataset. This repository provides the maintained data, standalone
compact-record codec and verification scripts.

A separate Figure 3 verification bundle provides the archived EAST seed and
15 replay states generated on 5 November 2026, with provenance and a standalone
pressure check in [`data/validation/figure3/`](data/validation/figure3/).

## Contents

| Path | Contents |
| --- | --- |
| `VEQDB-equilibria-20260920.tar.gz` | The release record set: 13,291 compact-JSON records, one file per accepted case, under 267 device directories. |
| `VEQDB-png-partial-20260920.tar.gz` | 45 rendered equilibrium images covering 16 devices. Partial and illustrative; not a complete image set. |
| `data/figure-data/` | Three supplementary reference records for the EAST, MAST-U and ITER `crossdev_*_p100` images. |
| `data/validation/figure3/` | Archived EAST seed, 15 pressure-replay records, state table, historical-summary excerpt and provenance. |
| `tokamak_data_2026-09-09.xlsx` | Device-survey workbook. |
| `FORMAT.md` | The compact-JSON format specification: key table, units, encoding rules and the verification commands. |
| `SHA256SUMS` | SHA-256 of each distributed download and supplementary reference record. |
| `manifests/VEQDB-equilibria-20260920.sha256` | One SHA-256 per record, 13,291 lines, with the archive's top-level directory as the path prefix. |
| `manifests/VEQDB-png-partial-20260920.records.json` | All 45 PNG-to-record associations, with file checksums and supplement provenance. |
| `scripts/` | Manuscript figure source, style definitions, the G-EQDSK projection snapshot, and record and PNG-pairing verification scripts. |
| `src/veqdb/` | Standalone compact-record decoder and encoder, including diagnostic recomputation. |
| `tests/` | Golden byte-parity, round-trip, validation and resolution-override tests. |
| `SOURCE_DATA.md` | Source repositories, fixed revisions and paths for the external EAST, CHEASE and EFIT inputs. |
| `data/source-provenance.json` | Machine-readable source identifiers, raw-file sizes and SHA-256 checksums. |

`VEQDB-equilibria-20260920.tar.gz` remains the authoritative initial record
set. Downloading that archive alone does not include the three supplementary
reference records.

## Reference Records and PNG Pairing

Of the 45 images in the PNG archive, 42 have a same-case record in the original
equilibrium archive. The other three map to these supplementary files:

| PNG path within `VEQDB-png-partial-20260920/` | JSON path in this repository |
| --- | --- |
| `EAST/crossdev_east_p100.png` | [`data/figure-data/crossdev_east_p100.compact.json`](data/figure-data/crossdev_east_p100.compact.json) |
| `MAST_Upgrade/crossdev_mast-u_p100.png` | [`data/figure-data/crossdev_mast-u_p100.compact.json`](data/figure-data/crossdev_mast-u_p100.compact.json) |
| `ITER_ITER-FEAT/crossdev_iter_p100.png` | [`data/figure-data/crossdev_iter_p100.compact.json`](data/figure-data/crossdev_iter_p100.compact.json) |

These records are copied without numerical changes from `FusionAlpha/veqdb-dev`
revision `eba15100d5c1027cc7ee8f7a01c917d4bf9bd94c`, at the same paths.
They are the `p100` reference-fit inputs used by `scripts/build_figures.py`,
with 100 active interior fitting parameters and 15th-order LCFS harmonics.
The EAST source is the experimental lower-single-null EFIT equilibrium for
discharge 48059 at 3.65 s. MAST-U and ITER are FreeGSNKE model references,
not experimental reconstructions; the ITER conversion uses the closed inset
contour at normalized flux 0.999, as described in the manuscript’s cross-device reconstruction section.

The compact records describe the fitted LCFS and plasma interior. They do not
supply the exterior field, divertor legs or strike points, and do not certify
preservation of an exact X-point topology. The manifest establishes complete
same-case file coverage and checksums; a pixel-identical rerender of the
existing PNGs has not been verified.

## Not Included In This Repository

The following are referenced by the manuscript and by `scripts/`, but are **not
distributed in this repository**:

- the per-device release archives, the per-device device manifests and the
  release summary;
- the device-survey sweep configuration;
- the retained G-EQDSK reference files; public source identifiers for EAST,
  CHEASE and EFIT are recorded in `SOURCE_DATA.md`;
- the full frozen figure-data collections and figure inputs used by
  `scripts/build_figures.py`; the Figure 3 verification bundle is documented
  separately above;
- the solver.

Consequently `scripts/build_figures.py` additionally requires the frozen inputs
listed above and the VEQPy package context used by the projection snapshot.
The `settings` function is present in
`scripts/gfile_projection/implementation/__init__.py`; that snapshot imports
solver components outside the standalone VEQDB package. The record-verification
scripts below and the compact-record codec have no such dependency.

## Verify The Figure 3 Pressure Range

After installing this repository, run

```sh
python scripts/verify_figure3_pressure.py
```

The command checks all 15 replay records through the public codec and
independent pressure integration, then prints a JSON report. The bundle README
explains the historical pressure diagnostic and the solver provenance.

## Reading A Record

Each record is a flat JSON object. Its key table, units and encoding rules are
specified in [`FORMAT.md`](FORMAT.md), which restates Table A.7 of the
manuscript and is self-contained: a reader needs nothing outside this
repository to interpret a record.

## Verifying A Record

The distributed downloads are covered by [`SHA256SUMS`](SHA256SUMS):

```sh
sha256sum -c SHA256SUMS
```

Individual records are covered by
[`manifests/VEQDB-equilibria-20260920.sha256`](manifests/VEQDB-equilibria-20260920.sha256),
one SHA-256 per record with the archive's top-level directory as the path
prefix:

```sh
tar xzf VEQDB-equilibria-20260920.tar.gz
sha256sum -c manifests/VEQDB-equilibria-20260920.sha256
```

Two checks run against the extracted records themselves, and need only the
records and NumPy:

```sh
python3 scripts/verify_field_order.py VEQDB-equilibria-20260920
python3 scripts/verify_record_consistency.py VEQDB-equilibria-20260920
```

`verify_field_order.py` checks the field order and field count of `FORMAT.md`
§2 against every record; `verify_record_consistency.py` checks the identities
of `FORMAT.md` §4, which exercise units, sign conventions and ordering but not
the fitting that produced the coefficients. `FORMAT.md` §5 carries the same
commands.

The same checks apply to the three supplementary reference records. The
PNG-pairing check reads the archives directly, so no extraction is needed:

```sh
python3 scripts/verify_field_order.py data/figure-data
python3 scripts/verify_record_consistency.py data/figure-data
python3 scripts/verify_png_records.py
```

The pairing check requires only Python's standard library. It verifies that
every PNG has exactly one manifest entry and an existing same-case JSON
record, and checks both file hashes. It does not render images or validate
the physical topology.

## Reconstruction

The standalone `veqdb.compact` interface decodes and re-encodes compact records
without importing VEQPy or Fusionprime. Install the package from this repository
with Python 3.12 or newer:

```sh
python3 -m pip install .
```

```python
from pathlib import Path
from veqdb import from_compact_json, to_compact_json

result = from_compact_json(
    Path("data/figure-data/crossdev_east_p100.compact.json").read_bytes()
)
payload = to_compact_json(result, precision=None)
```

`from_compact_json` returns a `FitResult` containing the raw reconstruction
roots, radial profiles and geometry coefficients. It is a compact-record
representation, not a full solver `Equilibrium` object. The optional `Nr`
argument interpolates the four radial profiles onto new Lobatto nodes, and
`Nt` sets the angular evaluation resolution. Both leave the shape coefficients
unchanged. The encoder recomputes the diagnostic
fields from the roots; stored diagnostics and unknown JSON fields are not
authoritative reconstruction inputs.

The default output precision is 10 significant digits. Use `precision=None`
to retain full floating-point precision. Re-encoding bytes produced in the same
numeric runtime is byte-stable. Recomputed diagnostics may differ in the last
floating-point digits across Python, library or platform builds, so the
full-precision golden comparison uses `rtol=1e-12` and `atol=1e-18` in each
diagnostic field's documented units. Reconstruction inputs, field order and
signed zero are still compared exactly, as are the 10- and 12-digit golden
encodings. These are test comparison tolerances; the codec's calculations,
stored records and default output precision are unchanged. The solver and
G-EQDSK projection are still separate: `scripts/gfile_projection/` retains the
projection implementation as source provenance and is not a second public
solver.

## Data Scope

VEQDB stores fixed-boundary plasma equilibria. It does not reconstruct the
exterior vacuum magnetic field and does not replace device-specific equilibrium
reconstruction.

## License and Citation

The VEQDB dataset is available under CC BY 4.0 through the
[Figshare deposit](https://doi.org/10.6084/m9.figshare.34485114).
The standalone codec retains its BSD 3-Clause licence in
[`LICENSES/veqpy-BSD-3-Clause.txt`](LICENSES/veqpy-BSD-3-Clause.txt).
See [`LICENSE`](LICENSE) for scope and
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for upstream attribution.

The external EAST, CHEASE and EFIT inputs are identified by source repository,
fixed revision and file path in [`SOURCE_DATA.md`](SOURCE_DATA.md).
Their raw-file checksums are recorded in
[`data/source-provenance.json`](data/source-provenance.json).
