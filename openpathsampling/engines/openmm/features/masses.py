@property
def masses(snapshot):
    """
    Returns
    -------
    masses : list of length n_atoms
        atomic masses (with simtk.unit attached)
    """
    try:
        simulation = snapshot.engine.simulation
    except AttributeError:
        # OpenMM snapshot not from simulation engine (setup snapshot,
        # probably)
        ops_topology = snapshot.topology
        topology = ops_topology.mdtraj.to_openmm()
        masses = [a.element.mass for a in topology.atoms()]
        return masses
    else:
        # standard case from simulation
        system = simulation.context.getSystem()
        n_particles = system.getNumParticles()
        masses = [system.getParticleMass(i)
                  for i in range(system.getNumParticles())]
        return masses
