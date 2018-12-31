from . import tools

class ClassIsSomething(object):
    def __init__(self, check_method):
        self.check_method = check_method
        self._true_set = set()
        self._false_set = set()

    def __call__(self, obj):
        key = obj.__class__
        if key in self._false_set:
            result = False
        elif key in self._true_set:
            result = True
        else:
            result = self.check_method(obj)
            set_for_obj = {True: self._true_set,
                           False: self._false_set}[result]
            set_for_obj.add(key)
        return result


_is_iterable_method = lambda obj: \
        tools.is_iterable(obj) and not tools.is_numpy_iterable(obj)

is_storage_iterable = ClassIsSomething(_is_iterable_method)
is_storage_mappable = ClassIsSomething(tools.is_mappable)
