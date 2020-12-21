from . import tools

class ClassIsSomething(object):
    """Method to test whether a class exhibits a given attribute.

    The idea here is that the class is expected to be immutable, so whether
    or not it exhibits the desired attribute is something we can cache. For
    example, a class will remain iterable, and looking up whether a given
    object's class is in the set of classes known to be iterable is much
    faster than using a function that checks whether it is iterable.
    """
    def __init__(self, check_method):
        self.check_method = check_method
        self._true_set = set()
        self._false_set = set()

    def force_true(self, cls):
        self._false_set.discard(cls)
        self._true_set.add(cls)

    def force_false(self, cls):
        self._true_set.discard(cls)
        self._false_set.add(cls)

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


def _is_iterable_method(obj):
    return tools.is_iterable(obj) and not tools.is_numpy_iterable(obj)

is_storage_iterable = ClassIsSomething(_is_iterable_method)
is_storage_mappable = ClassIsSomething(tools.is_mappable)
