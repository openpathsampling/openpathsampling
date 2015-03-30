
from ala_dipeptide_example import AlanineDipeptideTrajectorySimulator
from openpathsampling.snapshot import Snapshot
from openpathsampling.ensemble import LengthEnsemble

import time

from simtk.unit import femtoseconds, picoseconds, nanometers, kelvin, dalton
from simtk.unit import Quantity

import mdtraj as md

if __name__ == "__main__":
    options = {
                'temperature' : 300.0 * kelvin,
                'collision_rate' : 1.0 / picoseconds,
                'timestep' : 2.0 * femtoseconds,
                'nframes_per_iteration' : 10,
                'n_frames_max' : 5000,
                'start_time' : time.time(),
                'fn_initial_pdb' : "../data/Alanine_solvated.pdb",
                'platform' : 'CPU',
                'solute_indices' : range(22),
                'forcefield_solute' : 'amber96.xml',
                'forcefield_solvent' : 'tip3p.xml'
               }

    simulator = AlanineDipeptideTrajectorySimulator(
                    filename="smalltraj.nc",
                    topology_file="../data/Alanine_solvated.pdb",
                    options=options,
                    mode='create'
                    )

    simulator.equilibrate(5)
    snap = Snapshot(simulator.simulation.context)

    ensemble = LengthEnsemble(10)

    traj = simulator.generate(snap, [ensemble.can_append])

    mdtrajectory = traj.md()

    mdtrajectory.save_pdb("ala_small_traj.pdb")

