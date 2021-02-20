.. _dev_general_advice:

General Advice for Developers
=============================

If you want to add to OpenPathSampling, we're excited to have your
contribution. We consider OPS to be an :ref:`ecosystem <ops_ecosystem>` of
tools for path sampling. That means that you can contribute either by
directly contributing to the code, or by maintaining your addition as a
separate plugin.

One thing we recommend against is modifying your copy of the core code to do
something useful but different (and breaking other things). That makes it
impossible for you to share your improvements with others, and makes it hard
for you to keep your changes up to date with improvements we make.

There are several great ways to contribute to OpenPathSampling and share
your hard work with the rest of the world. Which one is best for you will
depend on both how 

Fix things
----------

We're not perfect! Like every code, ours has bugs, or small usability
improvements that should be made. Please, help us out! If you want to help,
but don't know where to start, try looking at the issues marked `"good first
issue"`_ on GitHub. Also, might consider improving our documentation. Like
most open source projects, our developers are focused on the code, and would
love to have contributions to the docs.

.. _"good first issue": https://github.com/openpathsampling/openpathsampling/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22

Contribute to Core
------------------

The OPS Core is the material in the `OpenPathSampling GitHub repository
<https://github.com/openpathsampling/openpathsampling>`_. This gets
installed whenever someone installs OPS. When we add major new functionality
to core, we say we're willing to take on the responsibility to maintain that
code. As such, we have stricter requirements on things that are added to
core.

If you want your contribution added to core, we require thorough tests: your
code should have 100% coverage. We also strongly prefer that you provide an
example of whatever you have implemented. We'll also be strict about our
coding style, and may ask you to restructure parts of your code to make it
more maintainable in the long run. We know this can be a burden, but since
we're taking on the responsibility for the code's upkeep, we need to know
that we can do that.

If you're not familiar with writing software tests, one of the OPS devs has
given a tutorial on `software testing in scientific programming
<https://training.e-cam2020.eu/files/5ca5ba6fe4b0fed490544a7b?dataset=5ca5b9dfe4b0fed490544a56&space=5ca35151e4b0fed490540623>`_.

For what it's worth, `code review will help you become a better programmer
<https://simpleprogrammer.com/why-code-reviews-make-better-code-teams/>`_.
Sure, it will take some effort, but that will pay off in the long run.

Create a plugin
---------------

If you don't want to go through the extra requirements to be added to core,
we're still happy to help you share your work with others! You can create a
separate repository that makes your code available. Remember, because of the
way Python allows you to import multiple libraries into a script, your code
can talk to OPS without needing to be in the OPS core library. We're happy
to announce your plugin as part of our page on :ref:`ops_ecosystem`. Submit
a pull request to the OPS library, changing the file at
``docs/ecosystem.rst`` to include a link to your package. We can't vouch for
these packages, but we can help raise their visibility.
