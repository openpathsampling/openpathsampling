# Adding new CV wrappers

OpenPathSampling is intended to make your life easier by interfacing with your
preferred tool for calculating collective variables. Several of those are
included in the `cv_wrappers` subpackage. To add a new CV wrapper, create a
file (module) with your classes, and import any public classes/functions in the
`cv_wrappers/__init__.py` file.
