import inspect

import logging
logger = logging.getLogger(__name__)

class DelayedLoader(object):
    def __get__(self, instance, owner):
        if instance is not None:
            obj = instance._lazy[self]
            if type(obj) is tuple:
                (store, idx) = obj
                return store[idx]
            elif hasattr(obj, '__subject__'):
                return obj.__subject__
            else:
                return obj
        else:
            return self

    def __set__(self, instance, value):
        if type(value) is tuple:
            instance._lazy[self] = value
        else:
            instance._lazy[self] = value

class StorableObject(object):
    """Mixin that allows an object to carry a .name property that can be saved

    It is not allowed to rename object once it has been given a name. Also
    storage usually sets the name to empty if an object has not been named
    before. This means that you cannot name an object, after is has been saved.
    """

    _base = None

    def idx(self, store):
        return store.index.get(self, None)

    def __init__(self):
        self._lazy = dict()
        self._uid = ''

    @property
    def cls(self):
        return self.__class__.__name__

    def save(self, storage):
        storage.save(self)

    @classmethod
    def base(cls):
        if cls._base is None:
            if cls is not StorableObject and cls is not StorableNamedObject:
                if StorableObject in cls.__bases__ or StorableNamedObject in cls.__bases__:
                    cls._base = cls
                else:
                    cls._base = cls.__base__.base()

        return cls._base

    @property
    def base_cls_name(self):
        return self.base().__name__

    @property
    def base_cls(self):
        return self.base()

    @classmethod
    def descendants(cls):
        return cls.__subclasses__() + \
               [g for s in cls.__subclasses__()  for g in s.descendants()]

    @staticmethod
    def objects():
        """
        Returns a dictionary of all subclasses
        """
        subclasses = StorableObject.descendants()

        return { subclass.__name__ : subclass for subclass in subclasses }

    @classmethod
    def args(cls):
        try:
            args = inspect.getargspec(cls.__init__)
        except TypeError:
            return []
        return args[0]

    _excluded_attr = []
    _exclude_private_attr = True
    _restore_non_initial_attr = True
    _restore_name = True

    def to_dict(self):
        excluded_keys = ['idx', 'json', 'identifier']
        return {
            key: value for key, value in self.__dict__.iteritems()
            if key not in excluded_keys
            and key not in self._excluded_attr
            and not (key.startswith('_') and self._exclude_private_attr)
        }

    def __hash__(self):
        return id(self) / 8

    @classmethod
    def from_dict(cls, dct):
        if dct is None:
            dct = {}
        try:
            init_dct = dct
            non_init_dct = {}
            if hasattr(cls, 'args'):
                args = cls.args()
                init_dct = {key: dct[key] for key in dct if key in args}
                non_init_dct = {key: dct[key] for key in dct if key not in args}

            obj = cls(**init_dct)

            if cls._restore_non_initial_attr:
                if len(non_init_dct) > 0:
                    for key, value in non_init_dct.iteritems():
                        setattr(obj, key, value)
            else:
                if cls._restore_name:
                    if 'name' in dct:
                        obj.name = dct['name']

        except TypeError as e:
            print dct
            print cls.__name__
            print e
            print args
            print init_dct
            print non_init_dct

        return obj

class StorableNamedObject(StorableObject):
    """Mixin that allows an object to carry a .name property that can be saved

    It is not allowed to rename object once it has been given a name. Also
    storage usually sets the name to empty if an object has not been named
    before. This means that you cannot name an object, after is has been saved.
    """

    def __init__(self):
        super(StorableNamedObject, self).__init__()
        self._name = ''
        self._name_fixed = False

    @property
    def default_name(self):
        return '[' + self.__class__.__name__ + ']'

    def fix_name(self):
        self._name_fixed = True

    @property
    def name(self):
        if self._name == '':
            return self.default_name
        else:
            return self._name

    @name.setter
    def name(self, name):
        if self._name_fixed:
            raise ValueError('Objects cannot be renamed to "%s" after is has been saved, it is already named "%s"' % ( name, self._name))
        else:
            self._name = name

    def named(self, name):
        """Set the name

        This is only for syntactic sugar and allow for chained generation

        Examples
        --------
        >>> import openpathsampling as p
        >>> full = p.FullVolume().named('myFullVolume')
        """
#        copied_object = copy.copy(self)
#        copied_object._name = name
#        if hasattr(copied_object, 'idx'):
#            copied_object.idx = dict()

        self._name = name

        return self

