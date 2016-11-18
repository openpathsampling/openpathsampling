from celery import Celery
import openpathsampling as paths
import openpathsampling.engines.openmm as peng
import simtk.openmm as mm

import simtk.unit as u


app = Celery()

@app.task
def run_steps(steps):
    template = peng.snapshot_from_pdb(
        "/Users/jan-hendrikprinz/Studium/git/openpathsampling/examples/data/Alanine_solvated.pdb")
    topology = peng.to_openmm_topology(template)

    # Generated using OpenMM Script Builder
    # http://builder.openmm.org

    forcefield = mm.app.ForceField(
        'amber96.xml',  # solute FF
        'tip3p.xml'     # solvent FF
    )

    # OpenMM System
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=mm.app.PME,
        nonbondedCutoff=1.0*u.nanometers,
        constraints=mm.app.HBonds,
        ewaldErrorTolerance=0.0005
    )

    # OpenMM Integrator
    integrator = mm.LangevinIntegrator(
        300 * u.kelvin,
        1.0 / u.picoseconds,
        2.0 * u.femtoseconds
    )
    integrator.setConstraintTolerance(0.00001)

    # Engine options
    options = {
        'n_steps_per_frame': 2,
        'n_frames_max': 5
    }

    engine = peng.Engine(
        template.topology,
        system,
        integrator,
        options=options
    )

    engine.initialize('CPU')

    traj = engine.generate(template, paths.LengthEnsemble(steps).can_append)

    return traj
