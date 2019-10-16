.. _how_use:

How I use OpenPathSampling
==========================

*David W.H. Swenson*

Usually these docs are written in the plural, and usually there isn't a
specific byline for it. That's because most things in the docs are just
facts, and are obviously agreed by all authors. However, this is about the
specific workflow that *I* like to use, and others may have different
opinions. That said, when I develop OPS, it is with my workflow in mind. So
this workflow is particularly well-supported.

The essence of my workflow is based on the idea of splitting the simulation
intro three stages: (1) setup; (2) sampling; and (3) analysis. In my
workflow, setup and analysis are done in an interactive Python environment
(e.g., Jupyter notebooks), while sampling is done in a simple script. There
are a few ideas/guidelines that lead to my workflow:

* OPS has a powerful storage system that helps you track provenance of data.
  Whenever possible use it! Load objects from storage instead of creating
  new ones. 
* Setting up a path sampling simulation can require care and attention. Set
  up simulations in an interactive environment (e.g., Jupyter notebook) and
  use the tools that OPS provides to perform sanity checks along the way.
* Analysis is best done interactively, using a tool like a Jupyter notebook
  that enables visualization of graphs and also stores the process in cells.
  This removes the overhead of redoing previous parts of the analysis every
  time a new idea to explore comes up. This also makes it easy to track the
  development of ideas, since the notebook acts as a sort of log.

The result is that my workflow tends to be something like this:

1. **Setup:** I perform the setup (defining engines, CVs, states/interfaces,
   networks, and move schemes) in a Jupyter notebook. Importantly, I provide
   these objects with names (which, for most objects, means using the
   ``.named()`` method at creation). I save all of these objects to a file,
   often called ``setup.nc``.
2. **Sampling:** The sampling process may run for a long time, and often on
   a remote machine. So I create a regular (non-interactive) Python script
   to run this. Note that this script is very simple: it basically just
   needs to create the path simulator, so frequently it just involves
   loading things from the ``setup.nc`` and plugging them into the path
   simulator class. In some cases (e.g., changing GPU device index) you may
   need to modify the engine, or even create a new engine and therefore a
   new move scheme. But it might be better to anticipate this and put
   multiple objects (with different names) in the ``setup.nc``.
3. **Analysis:** I do the analysis in Jupyter notebooks, for the reasons
   discussed above.

One feature that I often see users overlooking is the ability of OPS to
re-use simulation objects. For example, if you want to do a committor
analysis after performing TPS, you can load the load the states from a file
instead of redefining them, and this *guarantees* that your state
definitions are identical. The OPS storage subsystem has tools that enable
us to ensure that objects are identical across multiple simulations, even
with results stored in many different files. I choose to do this with a
single common ``setup.py``, although you could also load the states from the
output file of the TPS simulation.

****

Our examples tend to partially reflect this approach. In practice, I usually
put equilibration with the sampling, while it is often with the setup in the
examples. Also, the examples use a notebook for the sampling (to fit more
easily with the other stages, and to mix in text descriptions).
