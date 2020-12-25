import pytest
from .callable_codec import *
import dill

import numpy as np

_global_var = 0

def _known_module_func(self):
    import foo
    return foo.bar()

def _globals_using_function(self):
    return _global_var

class TestCallableCodec(object):
    def setup(self):
        self.codec = CallableCodec()
        self.functions = {
            'generic': lambda x: x.xyz[0][0],
            'known_module': np.sum,
            'known_submodule': np.linalg.matrix_power,
            'use_known_module': _known_module_func,
            'use_globals': _globals_using_function,
        }
        for func in ['generic', 'use_known_module', 'use_globals']:
            self.functions[func].__module__ = "__main__"

        dilled = {key: dill.dumps(func)
                  for key, func in self.functions.items()}
        self.dcts = {
            'generic': {
                '__callable_name__': '<lambda>',
                '_dilled': dilled['generic']
            },
            'known_module': {
                '__module__': 'numpy',
                '__callable_name__': 'sum'
            },
            'known_submodule': {
                '__module__': 'numpy.linalg',
                '__callable_name__': 'matrix_power'
            },
            'use_known_module': {
                '__callable_name__': '_known_module_func',
                '_dilled': dilled['use_known_module']
            }
            # use_globals never gets serialized
        }

    @pytest.mark.parametrize(
        'func,known_modules,allow_unknown',
        [
            ('generic', [], True),
            ('known_module', ['numpy'], True),
            ('use_known_module', ['foo'], True),
            ('use_known_module', [], False),
        ],
        ids=[
            'generic',
            'known_module_func',
            'use_known_module-allowed',
            'use_unknown_module-allowed',
        ]
    )
    def test_default_allowed(self, func, known_modules, allow_unknown):
        self.codec.required_modules = known_modules
        self.codec.only_allow_required_modules = allow_unknown
        results = self.codec.default(self.functions[func])
        # up to here was the smoke test; now we validate
        assert results == self.dcts[func]

    def test_default_error_use_unknown(self):
        self.codec.required_modules = []
        self.codec.only_allow_required_modules = True
        with pytest.raises(RuntimeError):
            results = self.codec.default(self.functions['use_known_module'])

    def test_default_error_unknown_module_function(self):
        self.codec.required_modules = []
        self.codec.only_allow_required_modules = True
        with pytest.raises(RuntimeError):
            results = self.codec.default(self.functions['known_module'])

    def test_default_error_use_globals(self):
        with pytest.raises(RuntimeError):
            results = self.codec.default(self.functions['use_globals'])

    def test_default_not_mine(self):
        obj = 'foo'
        assert self.codec.default(obj) is obj

    # TODO: tests for object_hook
    @pytest.mark.parametrize('safe', ['safemode', 'normal'])
    @pytest.mark.parametrize('func', ['generic', 'known_module',
                                      'known_submodule',
                                      'use_known_module'])
    def test_object_hook(self, func, safe):
        safemode = {'safemode': True, 'normal': False}[safe]
        self.codec.known_modules = ['openpathsampling', 'numpy']
        self.codec.safemode = safemode
        result = self.codec.object_hook(self.dcts[func])
        # up to here is smoke test; now we validate
        expected = None if safemode else self.functions[func]

        untestable = ['generic', 'use_known_module']
        # generic can't be tested for recreastion; use_known_module changes
        # the module and therefore can't be tested my memory location
        if not safemode and func not in untestable:
            assert result == expected

    def test_object_hook_not_mine(self):
        dct = {'foo': 'bar'}
        assert self.codec.object_hook(dct) is dct

    def test_settings_properties_setters(self):
        defaults = {'only_allow_required_modules': False,
                    'required_modules': [],
                    'safemode': False}
        assert self.codec.settings == defaults
        custom = {'only_allow_required_modules': True,
                  'required_modules': ['numpy'],
                  'safemode': True}
        codec = CallableCodec(custom)
        assert codec.settings == custom
        for key, value in custom.items():
            setattr(self.codec, key, value)

        assert self.codec.settings == custom
