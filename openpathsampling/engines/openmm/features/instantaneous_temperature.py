import simtk.openmm as mm
import simtk.unit as u

@property
def n_degrees_of_freedom(snapshot):
    # dof calculation taken from OpenMM's StateDataReporter
    n_spatial = 3  # can we get this programmatically?
    system = snapshot.engine.simulation.system
    n_particles = system.getNumParticles()
    dofs_particles = sum([n_spatial for i in range(n_particles)
                               if system.getParticleMass(i) > 0*u.dalton])
    dofs_constaints = system.getNumConstraints()
    dofs_motion_removers = 0
    if any(type(system.getForce(i)) == mm.CMMotionRemover 
           for i in range(system.getNumForces())):
        dofs_motion_removers += 3
    dofs = dofs_particles - dofs_constaints - dofs_motion_removers
    return dofs

@property
def instantaneous_temperature(snapshot):
    # TODO: this can be generalized as a feature that works with any
    # snapshot that has features for KE (in units of kB) and n_dofs
    state = snapshot.engine.simulation.context.getState(getEnergy=True)
    # divide by Avogadro b/c OpenMM reports energy/mole
    ke_in_energy = state.getKineticEnergy() / u.AVOGADRO_CONSTANT_NA
    ke_per_kB = ke_in_energy / u.BOLTZMANN_CONSTANT_kB

    dofs = snapshot.n_degrees_of_freedom

    temperature = 2 * ke_per_kB / dofs
    return temperature

