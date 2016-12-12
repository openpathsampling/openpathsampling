# Alanine dipeptide committor simulation

This sequence of notebooks includes an approach for calculating many
committor points after performing the TPS simulation of alanine dipeptide.
One point to consider is how one selects the initial points from which to
shoot those committor simulations. This particular example was developed in
the context of trying to create a way to approximate the committor, so the
important thing was to get several independent estimates of the committor,
and the statistical weights within the path ensemble were not important.
Because of this, we use an approach based on directly accessing the
`storage.snapshots`. To get the correct statistical weights for the paths,
you would need to loop over the `storage.steps`, and extract frames from the
trajectories in the `active` sample set.

* [`1_select_snapshots.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/misc/alanine_dipeptide_committor/1_select_snapshots.ipynb):
  Script to select initial snapshots.  Requires the output from the alanine
  dipeptide TPS example as input.
* [`2_committor_simulation.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/misc/alanine_dipeptide_committor/2_committor_simulation.ipynb):
  Run the committor for the selected snapshots.
* [`3_committor_analysis.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/misc/alanine_dipeptide_committor/3_committor_analysis.ipynb):
  Convert the simulation data to a simplified form which can be used for
  several different kinds of analysis. Note that this is a bit of a
  peculiarity of the specific project from which this example comes.
  Normally you would directly study the `ShootingPointAnalysis` object, as
  is done in the next notebook.
* [`4_analysis_help.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/misc/alanine_dipeptide_committor/4_analysis_help.ipynb)
  Analyze the final results. This shows how to build a
  `ShootingPointAnalysis` object from the previous output (which isn't
  noramlly needed) and then how to analyze that object.
