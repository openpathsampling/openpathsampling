import simtk.unit as u
import numpy as np

@property
def masses_per_mole(snapshot):
    """
    Returns
    -------
    masses_per_mole : list of simtk.unit.Quantity with length n_atoms
        atomic masses (with simtk.unit attached) in units of mass/mole
    """
    try:
        simulation = snapshot.engine.simulation
    except AttributeError:
        # OpenMM snapshot not from simulation engine (setup snapshot,
        # probably)
        ops_topology = snapshot.topology
        topology = ops_topology.mdtraj.to_openmm()
        masses_per_mole = [a.element.mass for a in topology.atoms()]
    else:
        # standard case from simulation
        system = simulation.context.getSystem()
        n_particles = system.getNumParticles()
        masses_per_mole = [system.getParticleMass(i)
                  for i in range(system.getNumParticles())]
    masses_per_mole = u.Quantity(
        value=np.array([m.value_in_unit(u.dalton)
                        for m in masses_per_mole]),
        unit=u.dalton
    )
    return masses_per_mole

@property
def masses(snapshot):
    """
    Returns
    -------
    masses : list of simtk.unit.Quantity with length n_atoms
        atomic masses (with simtk.unit attached) in units of mass
    """
    masses_per_mole = snapshot.masses_per_mole
    masses = masses_per_mole / u.AVOGADRO_CONSTANT_NA
    # [m / u.AVOGADRO_CONSTANT_NA for m in masses_per_mole]
    return masses
