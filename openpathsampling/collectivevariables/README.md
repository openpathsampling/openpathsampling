# Adding new CV wrappers

OpenPathSampling is intended to make your life easier by interfacing with your
preferred tool for calculating collective variables. Several of those are
included in the `collectivevariables` subpackage. To add a new CV wrapper,
create a file (module) with your classes, and import any public
classes/functions in the `collectivevariables/__init__.py` file.
