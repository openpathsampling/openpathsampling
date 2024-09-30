Input parameters for compile files
==================================

.. ifconfig:: HAS_OPS_CLI is True

    The following sections give information on the input file parameters for
    each type of object. The first group of categories are objects that can
    have top-level entries in the input file. For examples ``engines`` is
    one of the main heading in the input.

    .. toctree::
       :maxdepth: 1
       :glob:

       input/engines
       input/collective_variables
       input/volumes
       input/networks
       input/move_schemes

    The remaining categories of objects do not have top-level headings. They
    can only occur when an of this type is inside another object.

    .. toctree::
       :maxdepth: 1

       input/move_strategies
       input/shooting_point_selectors.rst
       input/interface_set.rst
       input/missing_type_information.rst

.. ifconfig:: HAS_OPS_CLI is False

    The OpenPathSampling CLI does not seem to be installed in this
    environment, so we were unable to build the detailed usage information
    for the ``compile`` command.
