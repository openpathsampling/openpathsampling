from openpathsampling.netcdfplus.base import StorableObject

from .object import ObjectStore

import logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class VariableStore(ObjectStore):
    def __init__(self, content_class, var_names):
        super(VariableStore, self).__init__(
            content_class,
            json=False
        )

        # TODO: determine var_names automatically from content_class
        # problem is that some decorators, e.g. using delayed loader
        # hide the actual __init__ signature and so we cannot determine
        # what variables to store. Could be 2.0

        if not issubclass(content_class, StorableObject):
            raise ValueError(('Content_class %s must be subclassed from '
                             'StorableObject') % content_class.__name__)

        var_names_class = self.content_class.args()[1:]

        # backwards compatibility. Reorder the stored var_names to comply
        # with the signature. Optional variables need to be at the end!!
        var_names_new = []
        for name in var_names_class:
            if name in var_names:
                var_names_new.append(name)
            else:
                break

        logger.info(
            'Creates VariableStore with variables %s and instatiated with %s' %
            (str(var_names_new), str(var_names))
        )

        self.var_names = var_names_new
        self._cached_all = False

    def to_dict(self):
        return {
            'content_class': self.content_class,
            'var_names': self.var_names
        }

    def _save(self, obj, idx):
        for var in self.var_names:
            self.write(var, idx, obj)

    def _load(self, idx):
        # kwargs = {var: self.vars[var][idx] for var in self.var_names}
        args = [self.vars[var][idx] for var in self.var_names]
        return self.content_class(*args)

    def initialize(self):
        super(VariableStore, self).initialize()

        # Add here the stores to be imported
        # self.create_variable('name', 'var_type')

    def all(self):
        self.cache_all()
        return self

    def cache_all(self, part=None):
        """Load all samples as fast as possible into the cache

        Parameters
        ----------
        part : list of int or `None`
            If `None` (default) all samples will be loaded. Otherwise the
            list of indices in `part` will be loaded into the cache

        """
        max_length = self.cache.size[0]
        max_length = len(self) if max_length < 0 else max_length

        if part is None:
            length = min(len(self), max_length)
            part = range(length)
        else:
            part = sorted(list(set(part)))
            length = min(len(part), max_length)
            part = part[:length]

        if not part:
            return

        # just in case we saved the var_names in another order and so we are
        # backwards compatible

        if not self._cached_all:
            data = zip(*[
                self.vars[var][part]
                for var in self.var_names
            ])

            [self.add_to_cache(idx, v) for idx, v in zip(part, data)]

            self._cached_all = True

    def add_to_cache(self, idx, data):
        if idx not in self.cache:
            # attr = {var: self.vars[var].getter(data[nn])
            #         for nn, var in enumerate(self.var_names)}
            obj = self.content_class(*data)
            self._get_id(idx, obj)

            # self.index[obj.__uuid__] = idx
            self.cache[idx] = obj
