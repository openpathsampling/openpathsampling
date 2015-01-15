def storable(super_class):
    super_class.default_storage = None
    _base_cls_name = super_class.__name__

    def _init(self, *args, **kwargs):
        super_class._init(self, *args, **kwargs)

        if 'idx' in kwargs:
            self.idx = kwargs['idx']
        else:
            self.idx = dict()

        if 'storage' in kwargs:
            self.default_storage = kwargs['storage']

        super_class.base_cls_name = _base_cls_name

    def _save(self, storage=None):
        if storage is None:
            storage = self.default_storage
        storage.save(self)

    super_class._init = super_class.__init__
    super_class.__init__ = _init

    super_class.save = _save

    # register as a base_class for storable objects

    def __descendents__(clazz, cls=None):
        if cls is None:
            cls = super_class

        d = dict()
        d.update(_recurse_subclasses(cls))

        return d

    def _recurse_subclasses(cls):
        d = dict()
        for cls in super_class.__subclasses__():
            if cls is not object and cls is not None:
                d[cls.__name__] = cls
                subs = cls.__subclasses__()
                if len(subs) > 0:
                    d.update(super_class.subclass_dict(super_class, cls))
        return d

    super_class.__descendents__ = classmethod(__descendents__)

    return super_class


# Tree support
def nestable(super_class):
    super_class.nestable = True
    return super_class

# Register a class to be creatable. Basically just a dict to match a classname to the actual class
# This is mainly for security so that we do not have to use globals to find classes

class_list = dict()

def creatable(super_class):
    class_list[super_class.__name__] = super_class

    super_class.creatable = True
    if not hasattr(super_class, '_excluded_attr'):
        super_class._excluded_attr = []

    if not hasattr(super_class, 'to_dict'):
        def _to_dict(self):
            excluded_keys = ['idx']
            return {key: value for key, value in self.__dict__.iteritems() if key not in excluded_keys and key not in self._excluded_attr and not key.startswith('_')}

        super_class.to_dict = _to_dict

    if not hasattr(super_class, 'from_dict'):
        def _from_dict(cls, my_dict = None):
            if my_dict is None:
                my_dict={}

            return cls(**my_dict)

        super_class.from_dict = classmethod(_from_dict)

    return super_class

class LoadedObject(object):
    pass

def dictable(super_class):
    """
    A class decorator that marks a class to be storable in the storage using a LoadedObject class.
    This object will have the same class name, the same dict, but none of the functions and will not have
    been initialized. If you want real objects use @creatable
    :param super_class: The class to be decorated
    :return: The decorated class
    """
    class_list[super_class.__name__] = super_class

    super_class.dictable = True
    if not hasattr(super_class, '_excluded_attr'):
        super_class._excluded_attr = []

    if not hasattr(super_class, 'to_dict'):
        def _to_dict(self):
            excluded_keys = ['idx']
            return {key: value for key, value in self.__dict__.iteritems() if key not in excluded_keys and key not in self._excluded_attr and not key.startswith('_')}

        super_class.to_dict = _to_dict

    if not hasattr(super_class, 'from_dict'):
        def _from_dict(cls, my_dict = None):
            if my_dict is None:
                my_dict={}

            obj = LoadedObject()

            for key, value in my_dict.iteritems():
                setattr(obj, key, value)

            setattr(obj, 'cls', cls.__name__)

            return obj

        super_class.from_dict = classmethod(_from_dict)

    return super_class
