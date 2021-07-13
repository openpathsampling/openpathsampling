import pytest
import openpathsampling as paths

from openpathsampling.experimental.parallel.dask_integration import \
        SerialScheduler, DaskDistributedScheduler

from openpathsampling.integration_tools import HAS_OPENMM, HAS_MDTRAJ
from openpathsampling.tests.test_gromacs_engine import has_gmx as HAS_GMX
from openpathsampling.tests.test_helpers import data_filename

from openpathsampling.engines import toy as toys
from openpathsampling.engines import openmm as ops_openmm
from openpathsampling.engines import gromacs as ops_gmx

try:
    import dask
    import distributed
except ImportError:
    HAS_DASK = False
else:
    HAS_DASK = True

if HAS_OPENMM:
    from simtk import openmm as mm
    from simtk import unit

if HAS_MDTRAJ:
    import mdtraj as md



from openpathsampling.experimental.storage.collective_variables import \
        MDTrajFunctionCV, CoordinateFunctionCV

# NOTE: things in here are not too wrapped up in pytest (e.g., as fixtures)
# so that they can easily be imported for use in development of other
# schedulers/engines.

class SimpleEngineTest:
    def __init__(self):
        self._engine = None
        self._snapshot = None
        self._universal_volume = None

    def _make_engine(self):
        raise NotImplementedError()

    @property
    def engine(self):
        if self._engine is None:
            self._engine = self._make_engine()
        return self._engine

    @property
    def snapshot(self):
        return self.engine.current_snapshot

    @property
    def universal_volume(self):
        if self._universal_volume is None:
            self._universal_volume = self._make_volume()
        return self._universal_volume

    def cleanup(self):
        pass

class FlatEngine(SimpleEngineTest):
    def _make_engine(self):
        pes = toys.LinearSlope([0.0], 0.0)
        integ = toys.LeapfrogVerletIntegrator(dt=0.1)
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        toy_opts = {'integ': integ,
                         'n_frames_max': 1000,
                         'n_steps_per_frame': 1}
        engine = toys.Engine(options=toy_opts, topology=topology)
        snapshot = toys.Snapshot(coordinates=[[0.0]],
                                 velocities=[[0.0]])
        engine.current_snapshot = snapshot
        return engine

    def _make_volume(self):
        x = CoordinateFunctionCV(lambda s: s.xyz[0][0]).named('x')
        volume = paths.CVDefinedVolume(x, -1e6, 1e6).named('near-infinity')
        return volume

class OpenMMADEngine(SimpleEngineTest):
    def _make_engine(self):
        pdb_filename = data_filename('ala_small_traj.pdb')
        snapshot = ops_openmm.snapshot_from_pdb(pdb_filename)

        forcefield = mm.app.ForceField('amber96.xml', 'tip3p.xml')
        system = forcefield.createSystem(
            mm.app.PDBFile(pdb_filename).topology,
            nonbondedMethod=mm.app.PME,
            nonbondedCutoff=1.0*unit.nanometers,
            constraints=mm.app.HBonds,
            ewaldErrorTolerance=0.0005
        )
        integrator = mm.LangevinIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picoseconds,
            2.0 * unit.femtoseconds
        )
        integrator.setConstraintTolerance(0.00001)

        options = {
            'n_steps_per_frame': 2,
            'n_frames_max': 10
        }
        engine = ops_openmm.Engine(
            snapshot.topology,
            system,
            integrator,
            options=options
        )
        engine.current_snapshot = snapshot
        return engine

    def _make_volume(self):
        return make_ad_universal_volume(self.engine.topology)

class GromacsADEngine(SimpleEngineTest):
    def _make_engine(self):
        pass

    def _make_volume(self):
        return make_ad_universal_volume(self.engine.topology)

    def cleanup(self):
        pass

def make_ad_universal_volume(topology):
    # we use a non-periodic approach here, and include more than the whole
    # period in our resulting volume
    psi = MDTrajFunctionCV(md.compute_dihedrals, topology,
                           indices=[[6, 8, 14, 16]]).named('psi')
    volume = paths.CVDefinedVolume(psi, -4, 4).named('all')
    return volume

def engine_run_5_steps(engine, universal_volume, initial_snapshot):
    # we include the volume to test that we're actually calculating the
    # volume (loading snapshot data)
    ensemble = (paths.LengthEnsemble(5)
                & paths.AllInXEnsemble(universal_volume))
    traj = engine.generate(initial_snapshot, ensemble.can_append)
    return traj

# @pytest.mark.parametrize('engine_name', ['openmm', 'gromacs', 'toy'])
@pytest.mark.parametrize('engine_name', ['toy', 'openmm'])
@pytest.mark.parametrize('scheduler_name', ['serial', 'dask'])
def test_scheduler_engine_integration(engine_name, scheduler_name):
    eligible_engine = {
        'openmm': HAS_OPENMM & HAS_MDTRAJ,
        'gromacs': HAS_GMX & HAS_MDTRAJ,
        'toy': True,
    }[engine_name]

    eligible_scheduler = {
        'dask': HAS_DASK,
        'serial': True,
    }[scheduler_name]

    if not (eligible_engine & eligible_scheduler):
        pytest.skip()

    scheduler_builder = {
        'serial': SerialScheduler,
        'dask': DaskDistributedScheduler,
    }[scheduler_name]
    scheduler = scheduler_builder()

    model_builder = {
        'openmm': OpenMMADEngine,
        'gromacs': GromacsADEngine,
        'toy': FlatEngine
    }[engine_name]

    model = model_builder()
    engine = model.engine
    volume = model.universal_volume

    scheduler_engine_test(scheduler, engine, volume)

def scheduler_engine_test(scheduler, engine, volume):
    task = scheduler.wrap_task(engine_run_5_steps)
    future = task(engine, volume, engine.current_snapshot)
    result = scheduler.get_result(future)
    assert isinstance(result, paths.Trajectory)
    assert len(result) == 5
    for snap in result:
        assert volume(snap)
