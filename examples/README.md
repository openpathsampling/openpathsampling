# Examples

This directory contains examples of various path sampling approaches as
applied to both toy models and simple real systems (mostly solvated alanine
dipeptide). 

* `alanine_dipeptide_tps`: Transition path sampling (both fixed length and
flexible length) on alanine dipeptide as a 2-state system. Compare with
[Bolhuis, Dellago, Chandler PNAS
2000](http://dx.doi.org/10.1073/pnas.100127697).
* `alanine_dipeptide_mstis`: Multiple state transition interface sampling
(including both multiple replica and single replica approaches) for alanine
dipeptide as a 6-state system. Compare with [Du, Marino, Bolhuis JCP
2011](http://dx.doi.org/10.1063/1.3644344).
* `toy_model_mstis`: Multiple state transition interface sampling (including
both multiple replica and single replica approaches) for a simple 3-state 2D
model.
* `toy_model_mistis`: Multiple interface set transition interface sampling for
a simple 3-state model.
* `tests`: Tests of various OpenPathSampling modules
* `misc`: Miscellaneous example code
* `resources`: Additional files used by the examples (initial PDB structures,
logging configurations, etc)
