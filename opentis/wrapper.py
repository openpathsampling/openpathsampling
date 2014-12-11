def storable(super_class):
    super_class.default_storage = None

    def _init(self, *args, **kwargs):
        super_class._initialize_netCDF(self, *args, **kwargs)

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

    super_class._initialize_netCDF = super_class.__init__
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

def jsonserializable(super_class):
    if not hasattr(super_class, 'to_dict'):
        def _to_dict(self):
            return { '_cls' : self.__class__.__name__, '_args' : self.__dict__ }

        super_class.to_dict = _to_dict

    if not hasattr(super_class, 'from_dict'):

        class Object(object):
            pass

        def _from_dict(self, d):
            obj = super_class.__new__()
            for key, value in d.args.iteritems:
                setattr(obj, key, value)

            return obj

        super_class.from_dict = _from_dict

    return super_class