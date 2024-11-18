import openpathsampling as paths
from ..simstore.dict_serialization_helpers import (
    tuple_keys_to_dict, tuple_keys_from_dict
)

import importlib
def import_class(full_classname_string):
    splitter = full_classname_string.split('.')
    module_name = ".".join(splitter[:1])
    class_name = splitter[-1]
    module = importlib.import_module(module_name)
    cls = getattr(module, class_name)
    return cls


def callable_cv_from_dict(cls, dct):
    kwargs = dct.pop('kwargs')
    dct.update(kwargs)
    obj = cls(**dct)
    cv_callable = obj.cv_callable
    if callable(cv_callable):
        return obj

    try:
        cv_callable['_marshal'] = cv_callable['_marshal']['bytes']
    except:
        pass

    cv_callable = paths.netcdfplus.ObjectJSON.callable_from_dict(cv_callable)
    obj.cv_callable = cv_callable
    return obj


def function_pseudo_attribute_to_dict(obj):
    dct = super(paths.netcdfplus.FunctionPseudoAttribute, obj).to_dict()
    cls = dct['key_class']
    dct['key_class'] = cls.__module__ + '.' + cls.__name__
    return dct

def function_pseudo_attribute_from_dict(cls, dct):
    key_class = import_class(dct['key_class'])
    dct['key_class'] = key_class
    return cls.from_dict(dct)

def from_dict_attr_to_class(from_dict, attr_name):
    def inner(cls, dct):
        class_ = import_class(dct[attr_name])
        dct[attr_name] = class_
        return from_dict(dct)
    return inner

_PREPATCH = {
    cls: {'to': getattr(cls, 'to_dict'),
          'from': getattr(cls, 'from_dict')}
    for cls in [
        paths.netcdfplus.FunctionPseudoAttribute,
        paths.TPSNetwork,
        paths.MISTISNetwork,
        paths.CallableCV,
    ]
}

_IS_PATCHED_SAVING = False
_IS_PATCHED_LOADING = False

def monkey_patch_saving(paths, with_old_cvs=True):
    global _IS_PATCHED_SAVING
    if _IS_PATCHED_SAVING:
        return paths

    if with_old_cvs:
        paths.netcdfplus.FunctionPseudoAttribute.to_dict = \
                function_pseudo_attribute_to_dict

    paths.TPSNetwork.to_dict = tuple_keys_to_dict(
        paths.TPSNetwork.to_dict, 'transitions'
    )
    paths.MISTISNetwork.to_dict = tuple_keys_to_dict(
        paths.MISTISNetwork.to_dict, 'input_transitions'
    )
    _IS_PATCHED_SAVING = True
    return paths

def monkey_patch_loading(paths, with_old_cvs=True):
    global _IS_PATCHED_LOADING
    if _IS_PATCHED_LOADING:
        return paths

    if with_old_cvs:
        paths.CallableCV.from_dict = classmethod(callable_cv_from_dict)
        paths.netcdfplus.FunctionPseudoAttribute.from_dict = \
                classmethod(from_dict_attr_to_class(
                    paths.netcdfplus.FunctionPseudoAttribute.from_dict,
                    attr_name='key_class'
                ))
                # classmethod(function_pseudo_attribute_from_dict)

    paths.TPSNetwork.from_dict = classmethod(tuple_keys_from_dict(
        paths.TPSNetwork.from_dict, 'transitions'
    ))
    paths.MISTISNetwork.from_dict = classmethod(tuple_keys_from_dict(
        paths.MISTISNetwork.from_dict, 'input_transitions'
    ))
    _IS_PATCHED_LOADING = True
    return paths

def monkey_patch_all(paths, with_old_cvs=True):
    paths = monkey_patch_saving(paths, with_old_cvs)
    paths = monkey_patch_loading(paths, with_old_cvs)
    return paths

def unpatch(paths, with_old_cvs=True):
    global _IS_PATCHED_SAVING
    global _IS_PATCHED_LOADING
    if not (_IS_PATCHED_LOADING or _IS_PATCHED_SAVING):
        return paths

    old_cv_modules = [
        paths.collectivevariable,
        paths.collectivevariables,
    ]
    unpatch_modules = (
        [paths.netcdfplus]
        + (old_cv_modules if with_old_cvs else [])
        + [paths]
    )

    for cls, old in _PREPATCH.items():
        if cls.to_dict != old['to']:
            cls.to_dict = old['to']
        if cls.from_dict != old['from']:
            cls.from_dict = old['from']

    # we loop over the modules twice as quick-and-dirty solution to avoid
    # figuring the DAG topological order for imports (should actually in
    # order leaf to root)
    for module in unpatch_modules * 2:
        importlib.reload(module)

    _IS_PATCHED_SAVING = False
    _IS_PATCHED_LOADING = False
    return paths
