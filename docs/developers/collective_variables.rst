.. _dev-collective-variables:

Collective Variables
====================

The collective variable in OpenPathSampling is essentially just a wrapper
around a function. Developers can extend collective variables in a couple
ways. First, is the very simple approach of creating a particular generally
useful function which you might like to wrap in an OPS CV. Second is the
more complicated idea of creating a general wrapper for a class functions,
e.g., from some molecular dynamics analysis library.

The first examples


Wrapping a function in a CV
---------------------------


Creating a wrapper class for a library
--------------------------------------

OPS doesn't bother defining any collective variables internally, because we
count on other codes to have a wide range and fast implementations. To make
it very easy to wrap functions from your favorite libraries in OPS CVs, we
create custom classes for various libraries. In particular, OPS supports
simple wrapping of MDTraj functions and of pyEmma "featurizers." If you'd
like to extend OPS to support simple wrapping of some other library, this
section is for you.
