# VEQDB

VEQDB is the release repository for a compact database of fixed-boundary
tokamak equilibria represented with MXH--Chebyshev geometry and physical-profile
roots. It accompanies the manuscript *VEQDB: A Compact and Reconstructible
Multi-Device Tokamak Equilibrium Database,
https://arxiv.org/abs/2609.23296*.

The initial release contains 13,291 accepted equilibria from 267 conventional
and spherical tokamak devices. The compact-JSON record set occupies 41.24 MB
before archive compression, and 17.47 MB in the distributed archive.

## Contents

| Path | Contents |
| --- | --- |
| `VEQDB-equilibria-20260920.tar.gz` | The release record set: 13,291 compact-JSON records, one file per accepted case, under 267 device directories. |
| `VEQDB-png-partial-20260920.tar.gz` | 45 rendered equilibrium images covering 16 devices. Partial and illustrative; not a complete image set. |
| `tokamak_data_2026-09-09.xlsx` | Device-survey workbook. |
| `scripts/` | Manuscript figure source, style definitions and the G-EQDSK projection snapshot. |

`VEQDB-equilibria-20260920.tar.gz` is the authoritative public record set.

## Not Included In This Repository

The following are referenced by the manuscript and by `scripts/`, but are **not
distributed in this repository**:

- the per-device release archives, the per-device device manifests and the
  release summary;
- the device-survey sweep configuration;
- the retained G-EQDSK reference files and their source metadata;
- the frozen figure-data records and figure inputs used by
  `scripts/build_figures.py`;
- the solver, and the compact-record codec.

Consequently `scripts/build_figures.py` cannot run from this repository alone:
it reads the figure inputs listed above, and it imports a `gfile_projection`
`settings` module that is not present here.

## Reading A Record

Each record is a flat JSON object. This repository does **not** include a format
specification, a key table or a unit table, and it does not include a checksum
manifest, so a reader cannot verify a record against a published reference using
this repository alone. The key semantics are defined in the manuscript.

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
