# Experimental features

The API for things in `openpathsampling.experimental` should be considered
fluid. This subpackage exists in order to make cutting-edge features available
to advanced users without requiring that these features be mature enough
that the API has solidified.

For most of OPS, we subscribe to [semantic versioning](http://semver.org).
The `experimental` subpackage provides an exception to that: the API here
can change at any time.

Once we're confident in the API for features in the `experimental`
subpackage, we'll move them into the main package. If an `experimental`
feature requires replacing part of the main package (thus changing the API),
it will be held in the `experimental` subpackage until the next major
version release of OPS.
