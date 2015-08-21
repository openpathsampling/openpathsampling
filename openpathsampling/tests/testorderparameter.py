"""
@author David W.H. Swenson
"""

from test_helpers import SimulationDuckPunch, data_filename

import mdtraj as md
import numpy as np

import openpathsampling.collectivevariable as op
import openpathsampling as paths
from simtk import unit as u
import time

from msmbuilder.featurizer import AtomPairsFeaturizer

def setUp():
    class Object():
        pass
    # Use the standard Alanine to generate snapshots to store for higher testing
    global this

    this = Object()

    this.options = {'temperature' : 300.0 * u.kelvin,
               'collision_rate' : 1.0 / u.picoseconds,
               'timestep' : 1.0 * u.femtoseconds,
               'nsteps_per_frame' : 1,
               'n_frames_max' : 5,
               'start_time' : time.time(),
               'fn_initial_pdb' : data_filename("ala_small_traj.pdb"),
               'platform' : 'fastest',
               'solute_indices' : range(22),
               'forcefield_solute' : 'amber96.xml',
               'forcefield_solvent' : 'tip3p.xml'
              }

    # create a template snapshot
    this.template_snapshot = paths.snapshot_from_pdb(data_filename("ala_small_traj.pdb"))

    # and an openmm engine
    this.engine = paths.OpenMMEngine(options=this.options, template=this.template_snapshot)
    this.engine.initialized = True

    # run a small trajectory of a few steps that can be used to save, etc...
    this.traj = this.engine.generate(this.template_snapshot, running=[paths.LengthEnsemble(5).can_append])
    this.mdtraj = this.traj.md()

class testCV_Function(object):

    def setUp(self):
        # reuse objects everytime
        for key, value in this.__dict__.iteritems():
            setattr(self, key, value)


    def teardown(self):
        pass

    def test_dihedral_op(self):
        """ Create a dihedral order parameter """
        psi_atoms = [6,8,14,16]
        dihedral_op = op.CV_MD_Function("psi", md.compute_dihedrals,
                                    indices=[psi_atoms])

        md_dihed = md.compute_dihedrals(self.mdtraj, indices=[psi_atoms])
        my_dihed =  dihedral_op(self.traj)

        np.testing.assert_allclose(md_dihed.reshape(md_dihed.shape[:-1]), my_dihed)

    def test_atom_pair_featurizer(self):
        """ Create an atom pair collectivevariable using MSMSBuilder3 """

        atom_pairs = [[0,1], [10,14]]
        atom_pair_op = op.CV_Featurizer("atom_pairs", AtomPairsFeaturizer, pair_indices=atom_pairs)

        md_distances = md.compute_distances(self.mdtraj, atom_pairs)

        my_distances = atom_pair_op( self.traj )

        np.testing.assert_allclose(md_distances, my_distances)