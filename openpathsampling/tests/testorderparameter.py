"""
@author David W.H. Swenson
"""

from test_helpers import SimulationDuckPunch, data_filename

import mdtraj as md
import numpy as np

import openpathsampling.trajectory as trajectory
from openpathsampling.storage import Storage
import openpathsampling.collectivevariable as op

from msmbuilder.featurizer import AtomPairsFeaturizer

class testCV_Function(object):

    def setUp(self):
        # setUp is just reading in some alanine dipeptide frames: this is an
        # ugly hack
        self.storage = Storage(
                               filename=data_filename("ala_small_traj.nc"),
                               mode="a"
                              )

        topol = md.load(data_filename("ala_small_traj.pdb")).top.to_openmm()
        self.simulation = SimulationDuckPunch(topol, None)

        trajectory.Trajectory.simulator = self
        trajectory.Trajectory.storage = self.storage


    def teardown(self):
        pass

    def test_dihedral_op(self):
        """ Create a dihedral order parameter """
        psi_atoms = [6,8,14,16]
        dihedral_op = op.CV_MD_Function("psi", md.compute_dihedrals,
                                    indices=[psi_atoms])

        mdtraj_version = self.storage.trajectories.load(0).md()
        md_dihed = md.compute_dihedrals(mdtraj_version, indices=[psi_atoms])
        traj = self.storage.trajectories.load(0)

        my_dihed =  dihedral_op( traj )

        np.testing.assert_allclose(md_dihed, my_dihed)

    def test_atom_pair_featurizer(self):
        """ Create an atom pair collectivevariable using MSMSBuilder3 """

        atom_pairs = [[0,1], [10,14]]

        atom_pair_featurizer = AtomPairsFeaturizer(atom_pairs)
        atom_pair_op = op.CV_Featurizer("atom_pairs", atom_pair_featurizer)

        mdtraj_version = self.storage.trajectories.load(0).md()
        md_distances = md.compute_distances(mdtraj_version, atom_pairs)
        traj = self.storage.trajectories.load(0)

        my_distances = atom_pair_op( traj )

        np.testing.assert_allclose(md_distances, my_distances)