import openpathsampling as paths


def split_md_storage(filename, update_cvs=True):
    """
    Creates two copies of the current storage. One containing trajectories and the other the rest

    Currently this makes only sense for storage with snapshots that use additional stores.
    Otherwise we need to store the full snapshots for CVs anyway and nothing is gained.

    """

    import openpathsampling.engines.openmm as peng

    storage_from = paths.AnalysisStorage(
        filename=filename
    )

    if not issubclass(storage_from.snapshots.snapshot.snapshot_class, peng.Snapshot):
        raise RuntimeError('Split only makes sense (for now) for storages with openmm.Snapshots')

    filename_base = '.'.join(filename.split('.')[:-1])

    filename_main = filename_base + '_main.nc'
    filename_data = filename_base + '_frames.nc'

    storage_main = paths.Storage(filename=filename_main, mode='w')
    storage_data = paths.Storage(filename=filename_data, mode='w')

    # save trajectories to data
    map(storage_data.trajectories.save, storage_from.trajectories)

    # mark all kinetics and statics as already saved in main to prevent auto save there
    map(storage_main.kinetics.remember, storage_from.kinetics)
    map(storage_main.statics.remember, storage_from.statics)

    for storage_name in [
        'steps',
        'pathmovers', 'topologies', 'networks', 'details', 'trajectories',
        'shootingpointselectors', 'engines', 'volumes',
        'samplesets', 'ensembles', 'transitions', 'pathmovechanges',
        'samples', 'pathsimulators', 'cvs'
    ]:
        map(getattr(storage_main, storage_name).save, getattr(storage_from, storage_name))

    storage_main.statics.index.clear()
    storage_main.kinetics.index.clear()
    storage_main.statics.cache.clear()
    storage_main.kinetics.cache.clear()

    storage_main.snapshots.cache.clear()

    if update_cvs:
        # make we now how to load snapshots if necessary
        storage_main.fallback = storage_from
        cvs = storage_from.cvs[:]

        for cv in cvs:
            storage_main.cvs.set_cache_store(cv)
#            storage_main.cvs.create_cache(cv)

        for cv in cvs:
            q = storage_from.snapshots.all()[10:11]
            s = list.__getitem__(q, 0)
            cv(q)
            storage_main.cvs.sync(cv)

    storage_main.close()
    storage_data.close()


def join_md_storage(filename_main, filename_data=None):
    if filename_data is None:
        filename_data = filename_main[:-7] + 'frames.nc'

    filename_to = filename_main[:-7] + 'joined.nc'

    storage_main = paths.Storage(
        filename=filename_main,
        mode='r'
    )

    storage_data = paths.Storage(
        filename=filename_data,
        mode='r'
    )

    storage_to = paths.Storage(
        filename_to,
        mode='w'
    )

    map(storage_to.trajectories.save, storage_data.trajectories)

    for storage_name in [
        'steps',
        'pathmovers', 'topologies', 'networks', 'details', 'trajectories',
        'shootingpointselectors', 'engines', 'volumes',
        'samplesets', 'ensembles', 'transitions', 'pathmovechanges',
        'samples', 'pathsimulators', 'cvs'
    ]:
        map(getattr(storage_to, storage_name).save, getattr(storage_main, storage_name))

    storage_data.clone()
    storage_main.clone()
    storage_to.clone()
