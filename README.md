# VEQDB

VEQDB is the release and reproducibility repository for a compact database of
fixed-boundary tokamak equilibria represented with MXH--Chebyshev geometry and
physical-profile roots. It accompanies the manuscript *VEQDB: A Compact and 
Reconstructible Multi-Device Tokamak Equilibrium Database, 
https://arxiv.org/abs/2609.23296*.

The initial release contains 13,291 accepted equilibria from 267 conventional
and spherical tokamak devices. The distributed compact-JSON archives occupy
18 MB; the same JSON records occupy 41 MB before archive compression.

## Contents

| Path | Contents |
| --- | --- |
| `data/release/` | 267 compact-record archives, device manifests and the release summary. |
| `data/device-table/` | Device-survey workbook and the GAQ sweep configuration. |
| `data/gfile/` | Retained G-EQDSK references, source metadata and 257 x 257 examples. |
| `data/figure-data/` | Frozen numerical inputs used to render the manuscript figures. |
| `data/figure-inputs/` | Frozen TCV and START compact records used by the manuscript figures. |
| `scripts/` | Manuscript figure source, style definitions and the G-EQDSK projection snapshot. |

`data/release/` is the authoritative public record set. The other directories
provide provenance and reproducibility inputs; they are not a second release
of the same equilibria.

## Reproduce the Release Checks

The release archives are compact JSON only. Each `*.manifest.json` identifies
the corresponding archive and records its accepted cases. The source data and
figure inputs can be inspected without a solver. Reconstructing equilibria,
reading or writing G-EQDSK, and regenerating a solved case use
[VEQPy](https://github.com/FusionAlpha/veqpy), pinned for this release at
commit `2b64c9a1bffddc25c176bd52231bc3b481715c6c`.

The manuscript plotting code is preserved in `scripts/build_figures.py`.
It can re-render figures from the frozen inputs after the VEQPy compact-record
and projection interfaces are installed. The G-EQDSK projection implementation
is retained in `scripts/gfile_projection/` so that its acceptance logic remains
auditable while that interface is promoted into VEQPy. The source file records
the current boundary: the snapshot is not presented as a second public solver.

## Data Scope

VEQDB stores fixed-boundary plasma equilibria. It does not reconstruct the
exterior vacuum magnetic field and does not replace device-specific equilibrium
reconstruction. Each retained G-EQDSK input remains traceable through source
metadata and the figure-data provenance fields.

## License and Citation

License and citation metadata will be added with the public release. Until
then, this development repository is intended for manuscript review and
reproducibility preparation.
