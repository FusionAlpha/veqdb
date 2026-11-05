# External equilibrium sources

The following public files are the external EAST, CHEASE and EFIT inputs
identified in the manuscript. Links pin the repository revisions used for
source identification.

| Input | Source repository and file | Use |
| --- | --- | --- |
| EAST `g048059.03650` | [Youjun Hu, GTAW](https://github.com/Youjunhu/GTAW/blob/7ae4929f2b4e56637f0462a9a186d92edb17f316/gfile/g048059.03650) | Figure 2 EAST reference, Figure 3 seed and the supplementary EAST compact record; discharge 48059 at 3.65 s |
| `CHEASE.geqdsk` | [FusionAlpha, VEQPy](https://github.com/FusionAlpha/veqpy/blob/2b64c9a1bffddc25c176bd52231bc3b481715c6c/data/CHEASE.geqdsk) | H-mode conversion benchmark |
| `EFIT.geqdsk` | [FusionAlpha, VEQPy](https://github.com/FusionAlpha/veqpy/blob/2b64c9a1bffddc25c176bd52231bc3b481715c6c/data/EFIT.geqdsk) | X-point conversion benchmark; file-header identifier 165908 at 3.00 s |

Repository revisions, download URLs, raw-file sizes and SHA-256 checksums are
also provided in [data/source-provenance.json](data/source-provenance.json).
These checksums identify the upstream files before the sign transformations
and projection steps described in the manuscript.

The GTAW source carries a GNU GPL software statement without a specified
version. The public VEQPy repository carries a BSD 3-Clause notice. Copies
and attribution are retained in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)
and [LICENSES/](LICENSES/).

The original files can be obtained from the linked repositories. This source
list covers the three identified external inputs; the README describes the
other inputs required for full figure regeneration.

The Figure 3 seed, the 5 November 2026 pressure replay records, and the
corresponding verification instructions are archived under
[`data/validation/figure3/`](data/validation/figure3/).
