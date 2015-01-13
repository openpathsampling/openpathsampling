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
def treebase(super_class):
    def _init(self, *args, **kwargs):
        super_class._init(self, *args, **kwargs)

        if 'idx' in kwargs:
            self.idx = kwargs['idx']
        else:
            self.idx = dict()

        if 'storage' in kwargs:
            self.default_storage = kwargs['storage']

    def _save(self, storage=None):
        if storage is None:
            storage = self.default_storage
        storage.save(self)

    super_class._init = super_class.__init__
    super_class.__init__ = _init

    super_class.save = _save

    def _example_cls_function(clazz, cls=None):
        if cls is None:
            cls = super_class

        return

    super_class.__descendents__ = classmethod(_example_cls_function)

    return super_class

# Register a class to be creatable. Basically just a dict to match a classname to the actual class
# This is mainly for security so that we do not have to use globals to find classes

class_list = dict()

def creatable(super_class):
    class_list[super_class.__name__] = super_class

    super_class.creatable = True

    if not hasattr(super_class, 'to_dict'):
        def _to_dict(self):
            excluded_keys = ['idx']
            return {key: value for key, value in self.__dict__.iteritems() if key not in excluded_keys}

        super_class.to_dict = _to_dict

    if not hasattr(super_class, 'from_dict'):
        def _from_dict(cls, my_dict = None):
            if my_dict is None:
                my_dict={}
            return cls(**my_dict)

        super_class.from_dict = classmethod(_from_dict)


    return super_class
