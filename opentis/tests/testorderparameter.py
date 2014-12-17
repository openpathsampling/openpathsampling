"""
@author David W.H. Swenson
"""

import os
import sys
from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import Skip, SkipTest
from test_helpers import AtomCounter, SimulationDuckPunch, data_filename

from nose.tools import assert_equal
import mdtraj as md

import opentis.trajectory as trajectory
from opentis.storage import Storage
import opentis.orderparameter as op


class testOP_Function(object):

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
        dihedral_op = op.OP_Function("psi", md.compute_dihedrals,
                                    trajdatafmt="mdtraj",
                                    indices=[psi_atoms])

        mdtraj_version = self.storage.trajectory.load(0).md()
        md_dihed = md.compute_dihedrals(mdtraj_version, indices=[psi_atoms])
        traj = self.storage.trajectory.load(0)

        my_dihed =  dihedral_op( traj )
        for (truth, beauty) in zip(md_dihed, my_dihed):
            assert_equal(truth, beauty)
