import types
import dis
import dill
from .serialization_helpers import has_uuid, do_import

GLOBALS_ERROR_MESSAGE="""
The function you tried to save uses information from global scope, which
can't be saved. This most frequently occurs when you use a package imported
outside the function. For more information, see:
http://openpathsampling.org/latest/topics/creating_collective_variables.html

The following globals need to be defined in your function, or passed in as a
parameter to your function:
"""

UNKNOWN_MODULES_ERROR_MESSAGE="""
The function you tried to save imports from unregistered modules, and your
code requires only registered modules, in order to prevent users from
loading modules that won't be available in other environements. You can
either register the unknown modules with ???, or you can turn this off
entirely with ???.

The following modules are not registered:
"""

class CallableCodec(object):
    """JSON codec for callables.

    Parameters
    ----------
    settings : Dict
        Dictionary with settings for how the codec should behave. Entries
        are:

        * ``'safemode'`` (bool): If True, this codec will not deserialize
          functions. Use this when working with untrusted data.
        * ``'required_modules'`` (List[str]): names of modules that can be
          expected to be present in the user's environment.
        * ``'only_allow_required_modules'`` (bool): If True, forbid
          serialization of any callable that requires importing a module
          that isn't in the ``'required_modules'`` list.

    """
    def __init__(self, settings=None):
        defaults = {
            'only_allow_required_modules': False,
            'required_modules': [],
            'safemode': False
        }
        if settings is None:
            settings = defaults
        else:
            for key in defaults:
                # use the original if it exists, otherwise use default
                settings[key] = settings.get(key, defaults[key])

        self.settings = settings

    @staticmethod
    def _error_message(message, problem_names):
        output = ""
        if problem_names:
            bulleted = "\n".join([" " + m for m in problem_names])
            output = message + bulleted
        return output

    @property
    def only_allow_required_modules(self):
        return self.settings['only_allow_required_modules']

    @only_allow_required_modules.setter
    def only_allow_required_modules(self, value):
        self.settings['only_allow_required_modules'] = value

    @property
    def required_modules(self):
        return self.settings['required_modules']

    @required_modules.setter
    def required_modules(self, value):
        self.settings['required_modules'] = value

    @property
    def safemode(self):
        return self.settings['safemode']

    @safemode.setter
    def safemode(self, value):
        self.settings['safemode'] = value

    def default(self, obj):
        only_req = self.only_allow_required_modules  # convenience
        if not has_uuid(obj) and callable(obj):
            errors = ""
            dct = {}
            # Case 1: this is a function from one of our requirements
            root_mod = obj.__module__.split('.')[0]
            if root_mod in self.required_modules:
                return {
                    '__module__': obj.__module__,
                    '__callable_name__': obj.__name__
                }
            elif obj.__module__ != "__main__" and only_req:
                # else implies root_mod not in required
                errors += self._error_message("bar" + UNKNOWN_MODULES_ERROR_MESSAGE,
                                              [obj.__module__])

            # Case 2: arbitrary function
            if obj.__module__ == "__main__":
                all_globals = dill.detect.globalvars(obj)
            else:
                all_globals = {}
            errors += self._error_message(GLOBALS_ERROR_MESSAGE,
                                          set(all_globals.keys()))

            if obj.__module__ == "__main__" and only_req:
                imported = [instr.argval
                            for instr in dis.get_instructions(obj)
                            if "IMPORT" in instr.opname]
                forbidden_imports = [mod for mod in imported
                                     if mod not in self.required_modules]
                errors += self._error_message(UNKNOWN_MODULES_ERROR_MESSAGE,
                                              forbidden_imports)

            if errors:
                raise RuntimeError("Cannot store function! \n\n" + errors)

            return {
                '__callable_name__': obj.__name__,
                '_dilled': dill.dumps(obj)
            }
        return obj

    def object_hook(self, dct):
        if '__callable_name__' in dct:
            if self.safemode:
                func = None
            elif '__module__' in dct:
                func = do_import(dct['__module__'], dct['__callable_name__'])
            elif '_dilled' in dct:
                func =  dill.loads(dct['_dilled'])
            else:  # pragma: no cover
                raise RuntimeError("Error reloading ",
                                   dct['__callable_name__'])
            return func

        return dct
