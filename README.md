# VEQDB

VEQDB is the release repository for a compact database of fixed-boundary
tokamak equilibria represented with MXH--Chebyshev geometry and physical-profile
roots. It accompanies the manuscript *VEQDB: A Compact and Reconstructible
Multi-Device Tokamak Equilibrium Database,
https://arxiv.org/abs/2609.23296*.

The initial release contains 13,291 accepted equilibria from 267 conventional
and spherical tokamak devices. The compact-JSON record set occupies 41.24 MB
before archive compression, and 17.47 MB in the distributed archive.

Three additional G-EQDSK-derived reference records are distributed separately
in `data/figure-data/`. Together with the 13,291 GAQ parameter-sampled records
in the original archive, they provide 13,294 records. The original archive and
its per-record checksum manifest are unchanged.

## Contents

| Path | Contents |
| --- | --- |
| `VEQDB-equilibria-20260920.tar.gz` | The release record set: 13,291 compact-JSON records, one file per accepted case, under 267 device directories. |
| `VEQDB-png-partial-20260920.tar.gz` | 45 rendered equilibrium images covering 16 devices. Partial and illustrative; not a complete image set. |
| `data/figure-data/` | Three supplementary reference records for the EAST, MAST-U and ITER `crossdev_*_p100` images. |
| `tokamak_data_2026-09-09.xlsx` | Device-survey workbook. |
| `FORMAT.md` | The compact-JSON format specification: key table, units, encoding rules and the verification commands. |
| `SHA256SUMS` | SHA-256 of each distributed download and supplementary reference record. |
| `manifests/VEQDB-equilibria-20260920.sha256` | One SHA-256 per record, 13,291 lines, with the archive's top-level directory as the path prefix. |
| `manifests/VEQDB-png-partial-20260920.records.json` | All 45 PNG-to-record associations, with file checksums and supplement provenance. |
| `scripts/` | Manuscript figure source, style definitions, the G-EQDSK projection snapshot, and record and PNG-pairing verification scripts. |

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
contour at normalized flux 0.999, as described in manuscript section 4.2.

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
- the retained G-EQDSK reference files and their source metadata;
- the remaining frozen figure-data records and figure inputs used by
  `scripts/build_figures.py`, apart from the three supplementary records above;
- the solver, and the compact-record codec.

Consequently `scripts/build_figures.py` cannot run from this repository alone:
it reads the figure inputs listed above, and it imports a `gfile_projection`
`settings` module that is not present here. The two record-verification scripts
named under Reading A Record below have no such dependency.

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

Reconstructing equilibria from these records requires the VEQPy compact-record
interface. [VEQPy](https://github.com/FusionAlpha/veqpy) is the solver and
G-EQDSK I/O boundary for VEQDB. That interface is **not present** in the VEQPy
revision this release was produced with,
`2b64c9a1bffddc25c176bd52231bc3b481715c6c`, so the pinned revision cannot
reconstruct a distributed record. The pinned revision does provide G-EQDSK
reading and writing.

The compact-record interface is published at
<https://hub.veloalpha.cn/gitea/suyuexinghen/veqpy-pkg>, revision
`586afb717099fcc1eb17ac5ffca4a10a95ae6d80`:

```sh
git clone https://hub.veloalpha.cn/gitea/suyuexinghen/veqpy-pkg.git
```

That revision decodes all 13,291 records in this release; the decoder entry
point is `veqpy.fitting.compact`.

The G-EQDSK projection implementation is retained in
`scripts/gfile_projection/` as source provenance so that its acceptance logic
remains auditable. It is not presented as a second public solver.

## Data Scope

VEQDB stores fixed-boundary plasma equilibria. It does not reconstruct the
exterior vacuum magnetic field and does not replace device-specific equilibrium
reconstruction.

## License and Citation

This repository contains no license file, and citation metadata has not been
added. No license is granted by this repository.
