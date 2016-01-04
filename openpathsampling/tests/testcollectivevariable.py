"""
@author David W.H. Swenson
"""

from test_helpers import data_filename, assert_close_unit

import mdtraj as md
import numpy as np

import openpathsampling.collectivevariable as op
import openpathsampling as paths

from msmbuilder.featurizer import AtomPairsFeaturizer

import time

class testCV_Function(object):

    def setUp(self):
        self.mdtraj = md.load(data_filename("ala_small_traj.pdb"))
        self.traj = paths.tools.trajectory_from_mdtraj(self.mdtraj)

    def teardown(self):
        pass

    def test_dihedral_op(self):
        """ Create a dihedral order parameter """
        psi_atoms = [6,8,14,16]
        dihedral_op = op.CV_MDTraj_Function("psi", md.compute_dihedrals, indices= [psi_atoms])

        md_dihed = md.compute_dihedrals(self.mdtraj, indices=[psi_atoms])
        my_dihed =  dihedral_op(self.traj)

        np.testing.assert_allclose(md_dihed.reshape(md_dihed.shape[:-1]), my_dihed, rtol=10**-6, atol=10**-10)

    def test_atom_pair_featurizer(self):
        """ Create an atom pair collectivevariable using MSMSBuilder3 """

        atom_pairs = [[0,1], [10,14]]
        atom_pair_op = op.CV_MSMB_Featurizer("atom_pairs", AtomPairsFeaturizer, pair_indices=atom_pairs)

        md_distances = md.compute_distances(self.mdtraj, atom_pairs)

        my_distances = atom_pair_op(self.traj)

        np.testing.assert_allclose(md_distances, my_distances, rtol=10**-6, atol=10**-10)

    def test_return_parameters_from_template(self):

        atom_pairs = [[0,1], [10,14]]
        atom_pair_op = op.CV_MSMB_Featurizer("atom_pairs", AtomPairsFeaturizer, pair_indices=atom_pairs)

        # little trick. We just predent the atom_pairs_op is a function we want to use
        # it cannot be stored though, but for from_template it is enough

        params = atom_pair_op.return_parameters_from_template(self.traj[0])

        assert params['cv_return_type'] == 'numpy.float32'
        assert params['cv_return_simtk_unit'] is None
        assert params['cv_return_shape'] == tuple([2])

    def test_register_as_attribute(self):
        cv = paths.CV_Function('x', lambda x: x.coordinates[:, 0])

        snap = self.traj[0]

        cv.register_as_attribute()

        ref_value = cv(snap)
        test_value = snap.x

        assert_close_unit(ref_value, test_value)

