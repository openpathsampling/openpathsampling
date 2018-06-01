import openpathsampling as paths

def callable_cv_from_dict(cls, dct):
    kwargs = dct.pop('kwargs')
    dct.update(kwargs)
    obj = cls(**dct)
    cv_callable = paths.netcdfplus.ObjectJSON.callable_from_dict(obj.cv_callable)
    obj.cv_callable = cv_callable
    return obj


