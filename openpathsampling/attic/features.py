__author__ = 'jan-hendrikprinz'

import simtk.unit as u
import functools
import re

class VariableFeature(object):

    @classmethod
    def add(cls, storage):
        self = cls()
        storage.features.append(self)

    def is_available(self):
        return True

    def fnc_type_check(self, var_type):
        return False

    def fnc_type_to_nc_type(self, var_type):
        return None

    def fnc_creation_wrapper(self, fnc_creation):
        def _wrapper(*args, **kwargs):
            fnc_creation(*args, **kwargs)

        return _wrapper

    def fnc_restoration_wrapper(self, fnc_restoration):
        def _wrapper(*args, **kwargs):
            fnc_restoration(*args, **kwargs)

        return _wrapper

    def fnc_type_delegates(self):
        setter = lambda x : x
        getter = lambda x : x
        return getter, setter


class SimtkUnits(VariableFeature):
    def add_simtk_support(self, var_name, unit):

        ncvar = self.variables[var_name]

        unit_instance = u.Unit({})
        symbol = 'none'

        if isinstance(unit, u.Unit):
            unit_instance = unit
            symbol = unit_instance.get_symbol()
        elif isinstance(unit, u.BaseUnit):
            unit_instance = u.Unit({unit : 1.0})
            symbol = unit_instance.get_symbol()
        elif type(unit) is str and hasattr(u, unit):
            unit_instance = getattr(u, unit)
            symbol = unit_instance.get_symbol()
        elif type(unit) is str and unit is not None:
            raise NotImplementedError('Unit by abbreviated string representation is not yet supported')

        json_unit = self.simplifier.unit_to_json(unit_instance)

        # store the unit in the dict inside the Storage object
        self.units[var_name] = unit_instance

        # Define units for a float variable
        setattr(ncvar,      'unit_simtk', json_unit)
        setattr(ncvar,      'unit', symbol)


    @staticmethod
    def compose(*functions):

        funcs = [func for func in functions if func is not None]

        print len(funcs)

        def _do_lambda():
            return functools.reduce(lambda f, g: lambda x: f(g(x)), funcs, lambda x: x)

        return _do_lambda()

    def parse_var_type(self, var_type):
        inner_type = None
        outer_fncs = []

        if '(' in var_type:
            red_var_type = var_type
            first_group = True
            while '(' in red_var_type:
                outer_params = list(re.search('([a-zA-Z\._]+)\([ ]*(?:([a-zA-Z0-9_\.]+)[ ]*,)*[ ]*([a-zA-Z0-9_\.]+)[ ]*\)', red_var_type).groups())
                outer_type = outer_params[0]
                if first_group:
                    if outer_params[1] is None:
                        inner_type = outer_params[2]
                        outer_params = [outer_params[0]] + outer_params[2:]
                    else:
                        inner_type = outer_params[1]
                    first_group = False
                outer_args = list(outer_params)[2:]
                outer_fncs.append([outer_type, outer_args])
                red_var_type = re.sub('(([a-zA-Z\._]+)\([ ]*(?:([a-zA-Z0-9_\.]+)[ ]*,)*[ ]*([a-zA-Z0-9_\.]+)[ ]*\))', 'inner', red_var_type)
        else:
            inner_type = re.search('[ ]*([a-zA-Z0-9_\.]+)[ ]*', var_type).group()

        inner = self.create_type_delegate(inner_type)

        if not inner['atomic']:
            raise ValueError('Type of inner function must be numeric!')

        get_list = [inner['getter']]
        set_list = [inner['setter']]

        for outer_type, outer_args in outer_fncs:
            outer = self.create_type_delegate(outer_type, *outer_args)

            if outer['atomic']:
                raise ValueError('Type of outer function must be')
            _get = outer['getter']
            _set = outer['setter']

            get_list.insert(0,_get)
            set_list.append(_set)

        getter = self.compose(*get_list)
        setter = self.compose(*set_list)

        nc_type = self.var_type_to_nc_type(inner_type)

        return {
            'getter' : getter,
            'setter' : setter,
            'nc_type' : nc_type,
            'outer' : outer_fncs,
            'inner' : inner_type,
            'get_list' : get_list,
            'set_list' : set_list
        }