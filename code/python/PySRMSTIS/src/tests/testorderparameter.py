"""
@author David W.H. Swenson
"""

import os
import sys
from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import Skip, SkipTest
from duckpunching import AtomCounter, SimulationDuckPunch

import mdtraj as md

sys.path.append(os.path.abspath('../'))
import snapshot
import trajectory
from storage import TrajectoryStorage
import orderparameter as op


class testOP_Function(object):

    def setUp(self):
        # setUp is just reading in some alanine dipeptide frames: this is an
        # ugly hack
        storage = TrajectoryStorage(    
                                    #topology="../data/Alanine_solvated.pdb",
                                    topology=None,
                                    filename="../data/trajectory.nc",
                                    mode="restore"  
                                    )
        topol = md.load("../data/Alanine_solvated.pdb").top.to_openmm()
        self.simulation = SimulationDuckPunch(topol, None)

        trajectory.Trajectory.simulator=self
        trajectory.Trajectory.storage = storage



    def teardown(self):
        pass

    def test_dihedral_op(self):
        """ Create a dihedral order parameter """
        phi_atoms = [7,9,15,17]

        dihedral_op = op.OP_Function("phi", md.compute_dihedrals,
                                    atom_indices=[phi_atoms])
        print dihedral_op( trajectory.Trajectory.load(1).md() )
