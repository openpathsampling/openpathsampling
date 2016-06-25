import openpathsampling as paths

def split_md_storage(self, filename):
    """
    Creates two copies of the current storage. One containing trajectories and the other the rest

    Currently this makes only sense for storage with snapshots that use additional stores.
    Otherwise we need to store the full snapshots for CVs anyway and nothing is gained.

    """

    import openpathsampling.engines.openmm as peng

    storage_from = paths.Storage(
        filename=filename,
        mode='r'
    )

    if not isinstance(storage.template, peng.Snapshot):
        raise RuntimeError('Split only makes sense (for now) for storages with openmm.Snapshots')

    filename_base = '.'.join(self.filename.split('.')[:-1])

    filename_main = filename_base + '_main.nc'
    filename_data = filename_base + '_frames.nc'

    storage_main = paths.Storage(filename=filename_main, template=self.template, mode='w')
    storage_data = paths.Storage(filename=filename_data, template=self.template, mode='w')

    map(storage_data.trajectories.save, self.trajectories)
    map(storage_main.kinetics.remember, self.trajectories)

    for storage_name in [
        'steps',
        'pathmovers', 'topologies', 'networks', 'details', 'trajectories',
        'shootingpointselectors', 'engines', 'volumes',
        'samplesets', 'ensembles', 'transitions', 'pathmovechanges',
        'samples', 'pathsimulators', 'cvs'
    ]:
        map(getattr(storage_main, storage_name).save, getattr(self, storage_name))

    storage_main.close()
    storage_data.close()