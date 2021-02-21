import simtk.openmm as mm
import simtk.unit as u

@property
def n_degrees_of_freedom(snapshot):
    """
    Returns
    -------
    n_degrees_of_freedom: int
        number of degrees of freedom in this system (after accounting for
        constraints)
    """
    return snapshot.engine.n_degrees_of_freedom()

@property
def instantaneous_temperature(snapshot):
    """
    Returns
    -------
    instantaneous_temperature : simtk.unit.Quantity (temperature)
        instantaneous temperature from the kinetic energy of this snapshot
    """
    # TODO: this can be generalized as a feature that works with any
    # snapshot that has features for KE (in units of kB) and n_dofs

    # if no engine, error here; don't get caught in try/except below
    engine = snapshot.engine
    try:
        old_snap = engine.current_snapshot
    except Exception:  # openmm doesn't use a custom exception class yet
        # Exception: Particle positions have not been set
        old_snap = None
    engine.current_snapshot = snapshot
    state = engine.simulation.context.getState(getEnergy=True)
    # divide by Avogadro b/c OpenMM reports energy/mole
    ke_in_energy = state.getKineticEnergy() / u.AVOGADRO_CONSTANT_NA
    ke_per_kB = ke_in_energy / u.BOLTZMANN_CONSTANT_kB

    dofs = snapshot.n_degrees_of_freedom

    temperature = 2 * ke_per_kB / dofs
    if old_snap is not None:
        engine.current_snapshot = old_snap
    return temperature

