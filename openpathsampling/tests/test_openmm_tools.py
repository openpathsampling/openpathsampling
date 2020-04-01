import string
import sys

import pytest

import openpathsampling as paths

if sys.version_info > (3,):
    maketrans = str.maketrans
else:
    maketrans = string.maketrans

# can't import * here bc something is named Test*!
from openpathsampling.engines.openmm.tools import (
    reduce_trajectory_box_vectors,
    reduced_box_vectors,
    load_trr
)

from .test_helpers import data_filename

from openpathsampling.engines.openmm import Snapshot
from openpathsampling.integration_tools import HAS_OPENMM, HAS_MDTRAJ

try:
    from simtk.unit import nanometer as nm
    from simtk.unit import picosecond as ps
except ImportError:
    HAS_SIMTK_UNIT = False
else:
    HAS_SIMTK_UNIT = True

import numpy.testing as npt
import numpy as np

try:
    from mdtraj.utils import box_vectors_to_lengths_and_angles
except ImportError:
    HAS_MDTRAJ = False
else:
    HAS_MDTRAJ = True

class BoxVectorWarning(RuntimeWarning):
    pass

# reduced form requires pair[0] >= 2.0*abs(pair[1])
_reduced_form_pairs = [('a_x', 'b_x'), ('a_x', 'c_x'), ('b_y', 'c_y')]

def get_index_pair(v, s):
    s = s.translate(maketrans("abcxyz", "012012"))
    vect, elem = tuple(map(int, s.split('_')))
    return vect, elem

def get_element(v, s):
    vect, elem = get_index_pair(v, s)
    return v[vect][elem]

def check_reduced_box_vectors(v):
    if v[0][0] <= 0.0 or v[1][1] <= 0.0 or v[2][2] <= 0.0:
        raise BoxVectorWarning("A diagonal element is not greater than zero")
    for left, right in _reduced_form_pairs:
        elem_l = get_element(v, left)
        elem_r = get_element(v, right)
        if elem_l < abs(2.0 * elem_r):
            raise BoxVectorWarning("Box vectors not reduced: "
                                   + left + " < |2.0*" + right
                                   + '|;  i.e., ' + str(elem_l)
                                   + ' < ' + str(abs(2.0 * elem_r)))


def mock_snapshot_with_box_vector(box):
    if not HAS_OPENMM:
        pytest.skip()
    return Snapshot.construct(coordinates=np.array([[]]) * nm,
                              velocities=np.array([[]]) * nm/ps,
                              box_vectors=box)


def test_reduced_box_vectors():
    pytest.importorskip("mdtraj")
    box = np.array([[ 6.70596027,  0.,          0.        ],
                    [ 0.,          6.70596027,  0.        ],
                    [ 3.35299015,  3.35299015,  4.74183893]])
    # TODO: replace this with a proper raises test in pytest
    try:
        check_reduced_box_vectors(box)
    except BoxVectorWarning:
        pass
    else:
        raise AssertionError("Box already reduced")

    if not HAS_OPENMM:
        pytest.skip()
    snap = mock_snapshot_with_box_vector(box * nm)
    reduced_box = reduced_box_vectors(snap).value_in_unit(nm)
    check_reduced_box_vectors(reduced_box)
    assert not all(box.flatten() == reduced_box.flatten())

    orig_lengths_angles = box_vectors_to_lengths_and_angles(*box)
    red_lengths_angles = box_vectors_to_lengths_and_angles(*reduced_box)

    orig_lengths = orig_lengths_angles[:3]
    red_lengths = red_lengths_angles[:3]
    npt.assert_allclose(orig_lengths, red_lengths, atol=1e-4)


def test_reduce_trajectory_box_vectors():
    box = np.array([[ 6.70596027,  0.,          0.        ],
                    [ 0.,          6.70596027,  0.        ],
                    [ 3.35299015,  3.35299015,  4.74183893]])
    snap_1 = mock_snapshot_with_box_vector(box)
    reduced_box = reduced_box_vectors(snap_1).value_in_unit(nm)
    snap_2 = mock_snapshot_with_box_vector(reduced_box)
    traj = paths.Trajectory([snap_1, snap_2])

    red_traj = reduce_trajectory_box_vectors(traj)
    orig_box = traj.box_vectors
    red_box = red_traj.box_vectors

    assert not all(orig_box[0].flatten() == red_box[0].flatten())
    npt.assert_allclose(red_box[0], red_box[1])
    for reduced in red_box:
        check_reduced_box_vectors(reduced.value_in_unit(nm))


def test_load_trr_with_velocities():
    if not (HAS_MDTRAJ and HAS_OPENMM):
        pytest.skip()
    box_vect_dir = "reduce_box_vects"
    gro = data_filename(box_vect_dir + "/dna.gro")
    trr = data_filename(box_vect_dir + "/dna.trr")
    traj = load_trr(trr, top=gro, velocities=True)
    for snap in traj:
        check_reduced_box_vectors(snap.box_vectors.value_in_unit(nm))
        assert np.count_nonzero(snap.velocities.value_in_unit(nm / ps)) > 0


def test_load_trr_no_velocities():
    if not (HAS_MDTRAJ and HAS_OPENMM):
        pytest.skip()
    box_vect_dir = "reduce_box_vects"
    gro = data_filename(box_vect_dir + "/dna.gro")
    trr = data_filename(box_vect_dir + "/dna.trr")
    traj = load_trr(trr, top=gro, velocities=False)
    for snap in traj:
        check_reduced_box_vectors(snap.box_vectors.value_in_unit(nm))
        assert np.count_nonzero(snap.velocities.value_in_unit(nm / ps)) == 0
