# Third-party notices

The compact codec under `src/veqdb/` was migrated from
`veqpy.fitting.compact` and its fitting diagnostics dependencies at
<https://hub.veloalpha.cn/gitea/suyuexinghen/veqpy-pkg>, revision
`586afb717099fcc1eb17ac5ffca4a10a95ae6d80`.

The VEQPy BSD 3-Clause license and copyright notice are retained in
`LICENSES/veqpy-BSD-3-Clause.txt`. This notice applies to the ported source
code; it does not change the licensing of the VEQDB data, figures, or other
repository contents.

The source modules document their narrower provenance: the compact schema and
writer were ported from `fusionprime.base.io.compact_equilibrium` (88d0fa0),
and the numerical kernels were lifted from the Fusionprime geometry and
equilibrium modules at the same revision. Their original attribution and
provenance comments are preserved in the migrated files.

## External equilibrium sources

The source repositories, fixed revisions and file paths for the EAST, CHEASE
and EFIT inputs are listed in [SOURCE_DATA.md](SOURCE_DATA.md).

- **GTAW, Youjun Hu**, revision
  `7ae4929f2b4e56637f0462a9a186d92edb17f316`: the original GNU GPL software
  statement and author attribution from `main.f90` are retained in
  [`LICENSES/GTAW-software-notice.txt`](LICENSES/GTAW-software-notice.txt).
  The statement does not specify a GPL version.
- **VEQPy, Ruohan Zhang**, revision
  `2b64c9a1bffddc25c176bd52231bc3b481715c6c`: the original repository notice
  is retained in
  [`LICENSES/VEQPy-upstream-BSD-3-Clause.txt`](LICENSES/VEQPy-upstream-BSD-3-Clause.txt).
  This is separate from the migrated codec's notice above.

Third-party materials retain their applicable upstream terms; the copied
software notices preserve their original attribution and scope.
