import openpathsampling as paths


def split_md_storage(filename):
    """
    Split storage into two files; trajectories and the rest

    Currently this makes only sense for storage with snapshots that use
    additional stores. Otherwise we need to store the full snapshots for
    CVs anyway and nothing is gained.

    """

    st_from = paths.AnalysisStorage(
        filename=filename
    )

    filename_base = '.'.join(filename.split('.')[:-1])

    filename_main = filename_base + '_main.nc'
    filename_data = filename_base + '_frames.nc'

    # `use_uuid=True`, otherwise we cannot later recombine the two!
    st_main = paths.Storage(filename=filename_main, mode='w')
    st_traj = paths.Storage(filename=filename_data, mode='w')

    st_main.snapshots.save(st_from.snapshots[0])
    st_traj.snapshots.save(st_from.snapshots[0])

    # this will tell the data storage not to save snapshots only a reference
    st_main.snapshots.only_mention = True

    # save trajectories to data
    map(st_traj.trajectories.save, st_from.trajectories)
    q = st_from.snapshots.all()
    cvs = st_from.cvs

    [cv(q) for cv in cvs]

    map(st_main.cvs.save, st_from.cvs)
    map(st_main.trajectories.mention, st_from.trajectories)

    for storage_name in [
        'steps', 'pathmovers', 'topologies', 'networks', 'details',
        'shootingpointselectors', 'engines', 'volumes', 'samples',
        'samplesets', 'ensembles', 'transitions', 'movechanges',
        'pathsimulators', 'cvs', 'interfacesets', 'msouters'
    ]:
        map(
            getattr(st_main, storage_name).save,
            getattr(st_from, storage_name)
        )

    st_main.close()
    st_traj.close()
    st_from.close()


def join_md_storage(filename_main, filename_data=None):
    if filename_data is None:
        filename_data = filename_main[:-7] + 'frames.nc'

    filename_to = filename_main[:-7] + 'joined.nc'

    st_main = paths.Storage(
        filename=filename_main,
        mode='r'
    )

    st_traj = paths.Storage(
        filename=filename_data,
        mode='r'
    )

    st_to = paths.Storage(
        filename_to,
        mode='w'
    )

    map(st_to.trajectories.save, st_traj.trajectories)

    for storage_name in [
        'steps',
        'pathmovers', 'topologies', 'networks', 'details', 'trajectories',
        'shootingpointselectors', 'engines', 'volumes',
        'samplesets', 'ensembles', 'transitions', 'movechanges',
        'samples', 'pathsimulators', 'cvs', 'interfacesets', 'msouters'
    ]:
        map(
            getattr(st_to, storage_name).save,
            getattr(st_main, storage_name))

    st_traj.close()
    st_main.close()
    st_to.close()
