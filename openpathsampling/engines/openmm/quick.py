from . import Engine, snapshot_from_pdb, to_openmm_topology
import simtk.openmm as mm
import simtk.unit as u

# TODO: Add this to example part create some helpers for
# it or make openmmtools easier to just setup a sytem to run on
# could also be part of openmm quick system setups


def create_simple_openmm_from_pdb(pdb_file):
    template = snapshot_from_pdb(pdb_file)
    topology = to_openmm_topology(template)

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
        'n_frames_max': 1000
    }

    engine = Engine(
        template.topology,
        system,
        integrator,
        options=options
    )

    return engine, template
