# Audiences for OpenPathSampling


## Users

At the basic level, we have users. These are the people who will be running
our code every day. The will primarily use the code in standard, predictable
ways: directly using one of the supported MD engines to run calculations of
one of the supported calculation types. 

We should assume that users are unfamiliar with programming, and may be
intimidated by the fact that run scripts are technically Python code.
Because of this, we should make the common, simple tasks as easy as possible
to use.

## Contributors

Contributors are comfortable with Python programming, but aren't interested
in diving into the details of how our code is implemented. However, they
might want to improve our code by add new modules which fall into
predictable categories: support for new MD engines, new order parameters,
new analysis tools, new calculation types, and perhaps even new ensembles.

To best help contributors, it is important that our code be carefully
modularized. A contributor who just wants to add a new MD engine should not
need to anything about other modules, except to the extent that it is
necessary to implement a new engine.

For the modules that we expect contributors to develop, we should clearly
document the API that needs to be implemented, and if possible, give minimum
examples (e.g., the ToyDynamics as an example of an MD engine.).

## Developers

At the developer level, we want the flexibility to implement whatever we
want. However, we can assume that developers are deeply involved in our
code, and don't need


## Trade-offs

In general, I think we should structure trade-offs such that 
