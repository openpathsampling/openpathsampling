# Alanine dipeptide TPS example

Running flexible- and fixed-length TPS on alanine dipeptide. This sequence
is split into two tracks: flexible-length and fixed-length. The first and
last files are common to both tracks, but the others are only needed for one
or the other.

These are related to, but not an exact reproduction of, the work done in
[Bolhuis, Dellago, and Chandler. PNAS **97**, 5877,
(2000)](http://dx.doi.org/10.1073/pnas.100127697).

- [`alanine_dipeptide_tps_first_traj.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_tps_first_traj.ipynb)
  Obtaining a transition trajectory using a high temperature integrator, and
  equilibrating that into a lower temperature path.

- [`alanine_dipeptide_tps_run.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_tps_run.ipynb)
  (Flexible-length track) Running flexible-length TPS (production run).

- [`alanine_dipeptide_fixed_tps_traj.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_fixed_tps_traj.ipynb)
  (Fixed-length track) Extending the transition trajectory to satsify the
  fixed-length path ensemble.

- [`alanine_dipeptide_fixed_tps_run.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_fixed_tps_run.ipynb)
  (Fixed-length track) Running fixed-length TPS (production run).

- [`alanine_dipeptide_tps_analysis.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_tps_analysis.ipynb)
  Analysis of both fixed and flexible TPS runs. Includes a number of
  advanced analysis techniques.


