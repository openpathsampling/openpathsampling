from __future__ import print_function
import sys
import warnings
from collections import namedtuple
from functools import wraps
from inspect import isclass

if sys.version_info > (3,):
    basestring = str

numpydoc_deprecation = """

.. deprecate:: {deprecated_in}
    {problem} {remedy}
"""

class Deprecation(object):
    """
    Parameters
    -----------
    problem : str
        description of what is being deprecated
    remedy : str
        description of how to avoid the deprecation warning
    remove_version : tuple
        version in which the deprecation will cause errors
    deprecated_in : tuple
        version when the deprecation was started
    warn_once : bool
        if True (default) only raise warning once per deprecation object
    """
    def __init__(self, problem, remedy, remove_version, deprecated_in,
                 warn_once=True):
        self.problem = problem
        self.remedy = remedy
        self.remove_version = remove_version
        self.deprecated_in = deprecated_in
        self.has_warned = False
        self.warn_once = warn_once
        self.str_replace = {
            'problem': self.problem,
            'remedy': self.remedy,
            'deprecated_in': version_tuple_to_string(self.deprecated_in),
            'version': version_tuple_to_string(self.remove_version),
            'OPS': "OpenPathSampling"
        }

    @property
    def message(self):
        return self.format_string("{problem} {remedy}")

    def format_string(self, string):
        result = string
        while any('{' + key + '}' in result for key in self.str_replace):
            result = result.format(**self.str_replace)
        return result
        # result = string.format(**self.str_replace).format(**self.str_replace)
        # result = string.format(problem=self.problem,
                               # remedy=self.remedy,
                               # deprecated_in=deprecated_in_str)
        # result = result.format(problem=self.problem,
                               # remedy=self.remedy,
                               # deprecated_in=deprecated_in_str,
                               # version=remove_version_str,
                               # OPS="OpenPathSampling")

    def warn(self):
        if not (self.has_warned and self.warn_once):
            warnings.warn(self.message, DeprecationWarning, stacklevel=2)
            self.has_warned = True

    def docstring_message(self, style='numpydoc'):
        string = {'numpydoc': numpydoc_deprecation}[style]
        return self.format_string(string)

    def __str__(self):  # pragma: no cover
        return "DEPRECATION: " + self.message


def update_docstring(thing_with_docstring, deprecation):
    # TODO: make a better version of this. Should come immediately after the
    # short description
    docs = thing_with_docstring.__doc__
    if thing_with_docstring.__doc__ is None:
        docs = ""
    else:
        docs = thing_with_docstring.__doc__
    return docs + deprecation.docstring_message()

def version_tuple_to_string(version_tuple):
    """
    Parameters
    ----------
    version_tuple : tuple of int
        the version, e.g, (1, 0) gives 1.0; (0, 2, 3) gives 0.2.3
    """
    return ".".join([str(x) for x in version_tuple])


# DEPRECATED THINGS SLATED FOR REMOVAL IN 2.0

OPENMMTOOLS_VERSION = Deprecation(
    problem="{OPS} {version} will require OpenMMTools 0.15 or later.",
    remedy="Please update OpenMMTools.",
    remove_version=(2, 0),
    deprecated_in=(0, 9, 6)
)

SAMPLE_DETAILS = Deprecation(
    problem="SampleDetails will be removed in {OPS} {version}.",
    remedy="Use generic Details class instead.",
    remove_version=(2, 0),
    deprecated_in=(0, 9, 3)
)

MOVE_DETAILS = Deprecation(
    problem="MoveDetails will be removed in {OPS} {version}.",
    remedy="Use generic Details class instead.",
    remove_version=(2, 0),
    deprecated_in=(0, 9, 3)
)

SAVE_RELOAD_OLD_TPS_NETWORK = Deprecation(
    problem="Old TPS networks will not be reloaded in {OPS} {version}.",
    remedy="This file may not work with {OPS} {version}.",
    remove_version=(2, 0),
    deprecated_in=(0, 9 ,3)
)

# has_deprecation and deprecate inspired by:
# https://stackoverflow.com/a/47441572/4205735
def has_deprecations(cls):
    for obj in [cls] + list(vars(cls).values()):
        if callable(obj) and hasattr(obj, '__new_docstring'):
            obj.__doc__ = obj.__new_docstring
            del obj.__new_docstring
    return cls

def deprecate(deprecation):
    """Decorator to deprecate a class/method

    Note
    ----
        Properties can be particularly challenging. First, you must put the
        @property decorator outside (above) the @deprecate decorator.
        Second, this does not (yet) change the docstrings of properties.
        However, it will raise a warning when the property is used.
    """
    def decorator(dep_obj):
        dep_obj.__new_docstring = update_docstring(dep_obj, deprecation)
        wrap_class = isclass(dep_obj)
        to_wrap = dep_obj.__init__ if wrap_class else dep_obj

        @wraps(to_wrap)
        def wrapper(*args, **kwargs):
            deprecation.warn()
            return to_wrap(*args, **kwargs)

        if wrap_class:
            dep_obj.__init__ = wrapper
            return dep_obj
        else:
            return wrapper
    return decorator

def list_deprecations(version=None, deprecations=None):
    """List deprecations that should have been removed by ``version``"""
    if deprecations is None:
        module = sys.modules[__name__]
        module_objects = [getattr(module, name) for name in dir(module)]
        deprecations = [obj for obj in module_objects
                        if isinstance(obj, Deprecation)]
    if version is None:
        version = "99.99"
    if version is not None:
        if isinstance(version, basestring):
            version = tuple(map(int, version.split('.')))
        deprecations = [d for d in deprecations if version >= d.remove_version]
    return deprecations


def print_deprecations(version=None):  # pragma: no cover
    # useful when preparing a release
    deprecations = list_deprecations(version)
    for dep in deprecations:
        print(dep)

