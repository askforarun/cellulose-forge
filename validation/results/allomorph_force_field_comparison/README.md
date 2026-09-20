# Matched GLYCAM and experimental OPLS-CM5 NPT comparison

Periodic I-alpha, I-beta, II, and III-I bulk crystals were independently built
with GLYCAM and the experimental OPLS-CM5 backend. Every system underwent PME
minimization, 10 ps NVT at 300 K, and 10 ps isotropic C-rescale NPT at 300 K and
1 bar. III-I used 6x4x4 cells; the other allomorphs used 4x4x4 cells.

All eight workflows completed without LINCS warnings or non-finite energies.
The table reports means over the final 5 ps of NPT:

| Allomorph | Force field | T (K) | P (bar) | Density (kg m-3) | Volume (nm3) |
|---|---|---:|---:|---:|---:|
| I-alpha | GLYCAM | 298.722 | -285.382 | 1499.441 | 22.9840 |
| I-alpha | OPLS-CM5 | 301.981 | 325.280 | 1573.863 | 21.8975 |
| I-beta | GLYCAM | 300.433 | 53.728 | 1527.919 | 45.1109 |
| I-beta | OPLS-CM5 | 301.529 | 94.612 | 1597.593 | 43.1442 |
| II | GLYCAM | 299.955 | -35.235 | 1483.899 | 46.4491 |
| II | OPLS-CM5 | 302.049 | 72.792 | 1541.889 | 44.7029 |
| III-I | GLYCAM | 298.991 | -74.982 | 1466.087 | 35.2602 |
| III-I | OPLS-CM5 | 301.810 | 2.629 | 1504.104 | 34.3695 |

The OPLS-AA defaults use combination rule 3 with `fudgeLJ = 0.5` and
`fudgeQQ = 0.5`. The converted GLYCAM topologies use combination rule 2 with
`fudgeLJ = 1.0` and `fudgeQQ = 1.0`, consistent with the explicit GLYCAM
1--4 treatment. Both values are recorded separately in `summary.json`.

Across these short runs, OPLS-CM5 produced 2.6--5.0% higher late-window density
than GLYCAM. This is an engineering observation, not a converged force-field
comparison. Ten picoseconds of isotropic pressure coupling cannot establish
equilibrium lattice parameters, anisotropic response, or physical accuracy.
Longer anisotropic NPT simulations with uncertainty analysis and comparison to
experimental lattice data are required before drawing scientific conclusions.
