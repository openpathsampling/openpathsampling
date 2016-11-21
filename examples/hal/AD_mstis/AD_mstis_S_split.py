import openpathsampling as paths
import os.path

import argparse
from tqdm import tqdm


def strip_snapshots(input_storage, output_storage=None):
    # =============================================================================
    # SPLITTING A SIMULATION
    # =============================================================================

    print """SPLITTING A SIMULATION"""

    base_filename = os.path.basename(input_storage)
    path_to_file = os.path.dirname(os.path.realpath(input_storage))


    print """1. Cache for reading"""

    storage = paths.Storage(input_storage)

    parts = base_filename.split('.')
    parts[-2] += '_strip'

    if output_storage is None:
        strip_filename = '.'.join(parts)
    else:
        strip_filename = output_storage

    # st_traj = paths.Storage('mstis_traj.nc', 'w')
    # st_data = paths.Storage('mstis_data.nc', 'w')
    st_strip = paths.Storage(os.path.join(path_to_file, strip_filename), 'w')

    st_strip.fallback = storage

    # save single snappshot
    st_strip.snapshots.save(storage.snapshots[0])
    # st_traj.snapshots.save(storage.snapshots[0])

    # cache proxies for all snapshots
    q = storage.snapshots.all()

    cvs = storage.cvs

    # make sure all values for CVs are precomputed
    print """2. Precompute CVs"""
    for cv in tqdm(cvs):
        _ = cv(q)

    # this will also switch the storage cache to the new file

    # save the cvs in the new file
    print """3. Copy CVs"""
    for cv in tqdm(cvs):
        st_strip.cvs.save(cv)

    # copy all trajectories with shallow copies of snapshots
    print """4. Copy trajectories w/o snapshots"""
    for traj in tqdm(storage.trajectories):
        st_strip.trajectories.mention(traj)

    # _ = map(st_strip.trajectories.mention, storage.trajectories)

    # _ = map(st_traj.trajectories.save, storage.trajectories)

    # save all the steps / all we want for analysis
    print """5. Copy steps"""
    for step in tqdm(storage.steps):
        st_strip.steps.save(step)

    # Some output
    print 'Original file:', storage.file_size_str
    print 'Stripped file:', st_strip.file_size_str
    # print 'Traj file:', st_traj.file_size_str
    print 'So we saved about %2.0f %%' % (
        (1.0 - st_strip.file_size / float(storage.file_size)) * 100.0)

    print """6. Closing"""

    st_strip.close()
    # st_traj.close()
    storage.close()

    print """7. Done"""


# ==============================================================================
#  MAIN
# ==============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create a copy of an OPS Storage without snapshots')

    parser.add_argument(
        'file',
        metavar='file.ipynb',
        help='the notebook to be checked',
        type=str)

    parser.add_argument(
        '-o', '--output', dest='output',
        type=str, default=300,
        help='The output file name. If not specified `_split` is appended.')

    args = parser.parse_args()
    inp_file = args.file
    out_file = args.output

    strip_snapshots(inp_file, out_file)
