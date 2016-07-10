import openpathsampling as paths


def split_md_storage(filename, update_cvs=True):
    """
    Split storage into two files; trajectories and the rest

    Currently this makes only sense for storage with snapshots that use
    additional stores. Otherwise we need to store the full snapshots for
    CVs anyway and nothing is gained.

    """

    storage_from = paths.AnalysisStorage(
        filename=filename
    )

    filename_base = '.'.join(filename.split('.')[:-1])

    filename_main = filename_base + '_main.nc'
    filename_data = filename_base + '_frames.nc'

    # `use_uuid=True`, otherwise we cannot later recombine the two!
    storage_main = paths.Storage(filename=filename_main, mode='w')
    storage_data = paths.Storage(filename=filename_data, mode='w')

    # this will tell the data storage not to save snapshots only a reference
    storage_data.snapshots.only_mention = True

    # save trajectories to data
    map(storage_data.trajectories.save, storage_from.trajectories)

    for storage_name in [
        'steps',
        'pathmovers', 'topologies', 'networks', 'details', 'trajectories',
        'shootingpointselectors', 'engines', 'volumes',
        'samplesets', 'ensembles', 'transitions', 'pathmovechanges',
        'samples', 'pathsimulators', 'cvs'
    ]:
        map(
            getattr(storage_main, storage_name).save,
            getattr(storage_from, storage_name)
        )

    if update_cvs:
        # tell the main store to load from the data storage if we need to
        storage_main.fallback = storage_from
        cvs = storage_from.cvs[:]

        for cv in cvs:
            storage_main.cvs.set_cache_store(cv)

        for cv in cvs:
            q = storage_from.snapshots.all()
            data = cv(q)
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
