# Alanine dipeptide TPS example

Running flexible- and fixed-length TPS on alanine dipeptide. This sequence
is split into two tracks: flexible-length and fixed-length. The first file
is common to both, and the final file, which compares the results of the
two, requires running both. However, the others are separate

If you only want to do flexible length TPS (sequence `a`), you should run
the files with numbers `1`, `2a`, and `3a`. If you only want to fixed length
TPS (sequence `b`), you should run the files `1`, `1b`, `2b`, and `3b`. To
run the file with number `4`, you must have run both `2a` and `2b`.

These are related to, but not an exact reproduction of, the work done in
[Bolhuis, Dellago, and Chandler. PNAS **97**, 5877,
(2000)](http://dx.doi.org/10.1073/pnas.100127697).

- [`AD_tps_1_trajectory.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_tps/AD_tps_1_trajectory.ipynb)
  (Both tracks) Obtaining a transition trajectory using a high temperature
  integrator, and equilibrating that into a lower temperature path.

- [`AD_tps_1b_trajectory.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_tps/AD_tps_1b_trajectory.ipynb)
  (Fixed-length track) Extending the transition trajectory to satsify the
  fixed-length path ensemble.

- [`AD_tps_2a_run_flex.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_tps/AD_tps_2a_run_flex.ipynb)
  (Flexible-length track) Running flexible-length TPS (production run).

- [`AD_tps_2b_run_fixed.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_tps/AD_tps_2b_run_fixed.ipynb)
  (Fixed-length track) Running fixed-length TPS (production run).

- [`AD_tps_3a_analyis_flex.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_tps/AD_tps_3a_analyis_flex.ipynb)
  (Flexible-length track) Analyzing the results of the flexible-length TPS
  run.

- [`AD_tps_3b_analyis_fixed.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_tps/AD_tps_3b_analyis_fixed.ipynb)
  (Flexible-length track) Analyzing the results of the flexible-length TPS
  run. Very similar to `3a`.

- [`AD_tps_4_advanced.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_tps/AD_tps_4_advanced.ipynb)
  (Both tracks) Advanced analysis techniques, much of which compares the
  behavior of the two approaches. 
