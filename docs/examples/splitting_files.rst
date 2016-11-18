Splitting Files
===============

OpenPathSampling saves all the information about the simulation, including
the coordinates and velocities of every snapshot. This makes it possible to
perform many different analyses later, even analyses that hadn't been
expected before the sampling.

However, this also means that the files can be very large, and frequently we
*don't* need all the coordinate and velocity data. This example will show
how to split the file into two: a large file with the coordinates and
velocities, and a smaller file with only the information needed to run the
main analysis. This allows you to copy the smaller file to a local drive and
perform the analysis interactively.

This particular example extends the :ref:`toy MSTIS example <toy-mstis>`. It
shows how to split the file, and then shows that the analysis still works.

-----

.. notebook:: examples/toy_model_mstis/toy_mstis_A1_split.ipynb
   :skip_exceptions:

-----

.. notebook:: examples/toy_model_mstis/toy_mstis_A2_split_analysis.ipynb
   :skip_exceptions:
