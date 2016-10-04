.. _which-network:

===========================
Which network should I use?
===========================

As of the release of OpenPathSampling 1.0, even people who are familiar with
path sampling are pretty new to the various multiple state approaches --- so
don't feel bad if that seems a bit overwhelming.

The most common network setups are described in the :ref:`Common Setups
<common-setups>` documentation. You should be familiar with the material in
that before continuing here.  This document mainly covers the differences
between :class:`MSTISNetworks <.MSTISNetwork>` and :class:`MISTISNetworks
<.MISTISNetwork>` and when to use each in a multiple-state TIS. The common
setups documentation contains sections for both :ref:`unidirectional 2-state
TIS <unidirectional-TIS>` and :ref:`bidirectional 2-state TIS
<bidirectional-TIS>`.

Once we get into the world of multiple states, it gets more complicated.
This is partly inevitable, but partly also because multiple state transition
sampling methods are relatively recent, dating back to the work of `Rogal
and Bolhuis`_.

.. _Rogal and Bolhuis: http://dx.doi.org/10.1063/1.3029696

The main approaches to multiple states in TIS are Multiple State TIS (MSTIS)
and `Multiple Interface Set TIS (MISTIS)`_. Both are fully supported by
OpenPathSampling. In MSTIS, there is one order parameter associated with
each state. In MISTIS, there is one order parameter associated with each
transition.

.. _Multiple Interface Set TIS (MISTIS): http://dx.doi.org/10.1063/1.4890037

This means that, depending on the nature of the network of transitions to be
studied, one or the other of these approaches might be better.

The advantages of the MISTIS approach are:

* it's easier to study networks where two transitions from a given state are
  very different
* you can selectively increase the sampling of rarer transitions
* you can easily study a subset of the total transition network

The advantages of the MSTIS approach are:

* it requires far fewer interfaces
* it is easier to set up
* the flux can be calculated from the minus interface 
* you can use adaptive methods to discover new states (not yet implemented
  in OPS)

As you might guess from these lists, if your system is easy to study, it
might be better to use an MSTIS network. However, the MISTIS network
provides more flexibility when dealing with complex reaction networks.

You can also create a custom network. Information on how to do that is
included in the documentation on :ref:`Transitions and Networks
<transitions-and-networks>`.
