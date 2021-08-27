"""
Snapshot code writing utilities.

Here we use the FeaturesCollection that has been attached to the snapshot
class in order to create a few methods for the class. In particular, we
write:

* ``__init__``
* ``copy_with_replacement``

The main reason we write these (as opposed to writing code in the superclass
that can be universal) is that we want to provide proper signatures and
docstrings for these methods for introspection and documentation. Indeed,
we keep the real logic for ``copy_with_replacement`` in the superclass.
"""

import inspect

def compile_function(func_name, code):
    code = compile("\n".join(code),
                   filename=f"<source:{func_name}>",
                   mode='exec')
    ns = {}
    exec(code, ns)
    return ns[func_name]

def _make_arguments(parameters):
    """Code for the (non-self) arguments for a function"""
    no_defaults = [p for p in parameters
                   if p.default is inspect.Parameter.empty]
    defaults = [p for p in parameters
                if p.default is not inspect.Parameter.empty]
    code = ", ".join([p.name for p in no_defaults])
    if len(defaults):
        code += ", "
    code += ", ".join([f"{p.name}={p.default_as_code}" for p in defaults])
    return code

def _call_super(method, parameters, assign_to=None):
    """Call super with ``parameters`` as pass-through parameters"""
    params = ", ".join([f"{p.name}={p.name}" for p in parameters])
    code = f"super().{method}({params})"
    if assign_to is not None:
        code = f"{assign_to} = {code}"
    return code

def _set_attribute(attr):
    """Code to ``self.attr`` to the name ``attr``"""
    return f"self.{attr} = {attr}"

def _def_function(func_name, parameters):
    """Code for function def"""
    return f"def {func_name}(self, {_make_arguments(parameters)}):"

INDENT = "    "

def make_init_code(features):
    """Code for the ``__init__`` method for a snapshot"""
    fname = '__init__'
    code = [
        _def_function(fname, features.parameters),
        INDENT + _call_super(fname, []),
    ]
    code += [
        INDENT + _set_attribute(p.name) for p in features.parameters
    ]
    return code

def make_copy_with_replacement_code(features):
    """Code for the ``copy_with_replacement`` method for a snapshot"""
    fname = 'copy_with_replacement'
    parameters = [p.copy_with_replacement(default=None)
                  for p in features.parameters]
    code = [
        _def_function(fname, parameters),
        INDENT + _call_super(fname, parameters, assign_to='new'),
        INDENT + "return new"
    ]
    return code

def make_init(features):
    code = make_init_code(features)
    return compile_function('__init__', code)

def make_copy_with_replacement(features):
    code = make_copy_with_replacement_code(features)
    return compile_function('copy_with_replacement', code)
