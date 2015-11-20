from base import StorableObject
from openpathsampling.netcdfplus import ObjectStore

class ExternalFile(StorableObject):
    """A Reference to an external file
    """

    def __init__(self, url):
        super(ExternalFile, self).__init__()
        self.url = url


class ExternalFileStore(ObjectStore):
    def __init__(self, file_class):
        super(ExternalFileStore, self).__init__(
            file_class,
            json=False
        )

        self.file_class = file_class
        self._cached_all = False

    def to_dict(self):
        return {}

    def _save(self, mcstep, idx):
        self.vars['url'][idx] = mcstep.url

    def _load(self, idx):
        return self.file_class(
            url=self.vars['url'][idx],
        )

    def _init(self, units=None):
        super(ExternalFileStore, self)._init()

        # New short-hand definition
        self.init_variable('url', 'str', chunksizes=(1,))

    def all(self):
        self.cache_all()
        return self

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            storage = self.storage

            urls = storage.variables[self.prefix + '_url'][:]

            [self.add_to_cache(idx, *v) for idx, v in enumerate(zip(urls))]

            self._cached_all = True

    def add_to_cache(self, idx, url):
        if idx not in self.cache:
            obj = ExternalFile(url=url)

            self.index[obj] = idx
            self.cache[idx] = obj