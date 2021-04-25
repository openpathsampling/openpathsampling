import inspect
import logging
import weakref
import uuid
from types import MethodType

import sys
if sys.version_info > (3, ):
    long = int

logger = logging.getLogger(__name__)

try:
    getfullargspec = inspect.getfullargspec
except AttributeError:
    getfullargspec = inspect.getargspec


class StorableObject(object):
    """Mixin that allows objects of the class to to be stored using netCDF+

    """

    _weak_cache = weakref.WeakKeyDictionary()
    _weak_index = 0

    _base = None
    _args = None

    observe_objects = False

    INSTANCE_UUID = list(uuid.uuid1().fields[:-1])
    CREATION_COUNT = 0
    ACTIVE_LONG = int(uuid.UUID(
            fields=tuple(
                INSTANCE_UUID +
                [CREATION_COUNT]
            )
        ))

    @staticmethod
    def get_uuid():
        StorableObject.ACTIVE_LONG += 2
        return StorableObject.ACTIVE_LONG

    def reverse_uuid(self):
        return self.__uuid__ ^ 1

    @staticmethod
    def ruuid(uid):
        return uid ^ 1

    def __init__(self):
        self.__uuid__ = StorableObject.get_uuid()

    @staticmethod
    def set_observer(active):
        """
        (De-)Activate observing creation of storable objects

        This can be used to track which storable objects are still alive and
        hence look for memory leaks and inspect caching. Use
        :meth:`openpathsampling.netcdfplus.base.StorableObject.count_weaks`
        to get the current summary of created objects

        Parameters
        ----------
        active : bool
            if `True` then observing is enabled. `False` disables observing.
            Per default observing is disabled.

        See Also
        --------
        :meth:`openpathsampling.netcdfplus.base.StorableObject.count_weaks`

        """
        if StorableObject.observe_objects is active:
            return

        if active:
            # activate and add __init__

            def _init(self):
                StorableObject._weak_cache[self] = StorableObject._weak_index
                StorableObject._weak_index += 1

            StorableObject.__init__ = MethodType(_init, None, StorableObject)
            StorableObject.observe_objects = True

        if not active:
            del StorableObject.__init__

    @staticmethod
    def count_weaks():
        """
        Return number of objects subclassed from StorableObject still in memory

        This includes objects not yet recycled by the garbage collector.

        Returns
        -------
        dict of str : int
            the dictionary which assigns the base class name of each references
            objects the integer number of objects still present

        """
        summary = dict()
        complete = list(StorableObject._weak_cache)
        for obj in complete:
            name = obj.base_cls_name
            summary[name] = summary.get(name, 0) + 1

        return summary

    def idx(self, store):
        """
        Return the index which is used for the object in the given store.

        Once you store a storable object in a store it gets assigned a unique
        number that can be used to retrieve the object back from the store. This
        function will ask the given store if the object is stored if so what
        the used index is.

        Parameters
        ----------
        store : :class:`openpathsampling.netcdfplus.ObjectStore`
            the store in which to ask for the index

        Returns
        -------
        int or None
            the integer index for the object of it exists or `None` else

        """
        if hasattr(store, 'index'):
            return store.index.get(self, None)
        else:
            return store.idx(self)

    @property
    def cls(self):
        """
        Return the class name as a string

        Returns
        -------
        str
            the class name

        """
        return self.__class__.__name__

    @classmethod
    def base(cls):
        """
        Return the most parent class actually derived from StorableObject

        Important to determine which store should be used for storage

        Returns
        -------
        type
            the base class
        """
        if cls._base is None:
            if cls is not StorableObject and cls is not StorableNamedObject:
                if StorableObject in cls.__bases__ \
                        or StorableNamedObject in cls.__bases__:
                    cls._base = cls
                else:
                    if hasattr(cls.__base__, 'base'):
                        cls._base = cls.__base__.base()
                    else:
                        cls._base = cls

        return cls._base

    def __hash__(self):
        return self.__uuid__ & 1152921504606846975

    def __eq__(self, other):
        if self is other:
            return True

        if hasattr(other, '__uuid__'):
            return self.__uuid__ == other.__uuid__

        return NotImplemented

    @property
    def base_cls_name(self):
        """
        Return the name of the base class

        Returns
        -------
        str
            the string representation of the base class

        """
        return self.base().__name__

    @property
    def base_cls(self):
        """
        Return the base class

        Returns
        -------
        type
            the base class

        See Also
        --------
        :func:`base()`

        """
        return self.base()

    @classmethod
    def descendants(cls):
        """
        Return a list of all subclassed objects

        Returns
        -------
        list of type
            list of subclasses of a storable object
        """
        return cls.__subclasses__() + \
            [g for s in cls.__subclasses__() for g in s.descendants()]

    @staticmethod
    def objects():
        """
        Returns a dictionary of all storable objects

        Returns
        -------
        dict of str : type
            a dictionary of all subclassed objects from StorableObject.
            The name points to the class
        """
        subclasses = StorableObject.descendants()
        return {subclass.__name__: subclass for subclass in subclasses
                if not subclass.__module__.startswith(
                    'openpathsampling.experimental.storage'
                )}

    @classmethod
    def args(cls):
        """
        Return a list of args of the `__init__` function of a class

        Returns
        -------
        list of str
            the list of argument names. No information about defaults is
            included.

        """
        try:
            args = getfullargspec(cls.__init__)
        except TypeError:
            return []
        return args[0]

    _excluded_attr = []
    _included_attr = []
    _exclude_private_attr = True
    _restore_non_initial_attr = True
    _restore_name = True

    def to_dict(self):
        """
        Convert object into a dictionary representation

        Used to convert the dictionary into JSON string for serialization

        Returns
        -------
        dict
            the dictionary representing the (immutable) state of the object

        """
        excluded_keys = ['idx', 'json', 'identifier']
        keys_to_store = {
            key for key in self.__dict__
            if key in self._included_attr or (
                key not in excluded_keys and
                key not in self._excluded_attr and
                not (key.startswith('_') and self._exclude_private_attr)
            )
        }
        return {
            key: self.__dict__[key] for key in keys_to_store
        }

    @classmethod
    def from_dict(cls, dct):
        """
        Reconstruct an object from a dictionary representaiton

        Parameters
        ----------
        dct : dict
            the dictionary containing a state representaion of the class.

        Returns
        -------
        :class:`openpathsampling.netcdfplus.StorableObject`
            the reconstructed storable object
        """
        if dct is None:
            dct = {}

        if hasattr(cls, 'args'):
            args = cls.args()
            init_dct = {key: dct[key] for key in dct if key in args}
            try:
                obj = cls(**init_dct)

                if cls._restore_non_initial_attr:
                    non_init_dct = {
                        key: dct[key] for key in dct if key not in args}

                    if len(non_init_dct) > 0:
                        for key, value in non_init_dct.items():
                            setattr(obj, key, value)

                return obj

            except TypeError as e:
                if hasattr(cls, 'args'):
                    err = (
                        'Could not reconstruct the object of class `%s`. '
                        '\nStored parameters: %s \n'
                        '\nCall parameters: %s \n'
                        '\nSignature parameters: %s \n'
                        '\nActual message: %s'
                    ) % (
                        cls.__name__,
                        str(dct),
                        str(init_dct),
                        str(cls.args),
                        str(e)
                    )
                    raise TypeError(err)
                else:
                    raise

        else:
            return cls(**dct)


class StorableNamedObject(StorableObject):
    """Mixin that allows an object to carry a .name property that can be saved

    It is not allowed to rename an object once it has been given a name. Also
    storage usually sets the name to empty if an object has not been named
    before. This means that you cannot name an object, after is has been saved.
    """

    def __init__(self):
        super(StorableNamedObject, self).__init__()
        self._name = ''
        self._name_fixed = False

    @property
    def default_name(self):
        """
        Return the default name.

        Usually derived from the objects class

        Returns
        -------
        str
            the default name

        """
        return '[' + self.__class__.__name__ + ']'

    def fix_name(self):
        """
        Set the objects name to be immutable.

        Usually called after load and save to fix the stored state.
        """
        self._name_fixed = True

    @property
    def name(self):
        """
        Return the current name of the object.

        If no name has been set a default generated name is returned.

        Returns
        -------
        str
            the name of the object
        """
        if self._name == '':
            return self.default_name
        else:
            return self._name

    @name.setter
    def name(self, name):
        if self._name_fixed:
            raise ValueError((
                'Objects cannot be renamed to `%s` after is has been saved, '
                'it is already named `%s`') %
                (name, self._name))
        else:
            if name != self._name:
                self._name = name
                logger.debug(
                    'Nameable object is renamed from `%s` to `%s`' %
                    (self._name, name))

    @property
    def is_named(self):
        """True if this object has a custom name.

        This distinguishes default algorithmic names from assigned names.
        """
        return self._name != ""

    def named(self, name):
        """Name an unnamed object.

        This only renames the object if it does not yet have a name. It can
        be used to chain the naming onto the object creation. It should also
        be used when naming things algorithmically: directly setting the
        .name attribute could override a user-defined name.

        Parameters
        ----------
        name : str
            the name to be used for the object. Can only be set once

        Examples
        --------
        >>> import openpathsampling as p
        >>> full = p.FullVolume().named('myFullVolume')

        """
        if not self.is_named and not self._name_fixed:
            self.name = name
        return self


def create_to_dict(keys_to_store):
    def to_dict(self):
        return {key: getattr(self, key) for key in keys_to_store}

    return to_dict
