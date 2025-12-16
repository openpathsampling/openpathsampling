import pytest
import openpathsampling as paths
import copy
import numpy as np

from openpathsampling.engines.openmm.mcengine import *

### FIXTURES FOR TEST SETUP ################################################

# TODO: move this fixture to a global conftest.py -- this should be reused
@pytest.fixture(scope='session')
def ad_vacuum_testsystem():
    testsystems = pytest.importorskip('openmmtools.testsystems')
    return testsystems.AlanineDipeptideVacuum()

@pytest.fixture(scope='session')
def thermodynamic_state(ad_vacuum_testsystem):
    states = pytest.importorskip('openmmtools.states')
    unit = pytest.importorskip('simtk.unit')
    return states.ThermodynamicState(
        system=ad_vacuum_testsystem.system,
        temperature=298.0 * unit.kelvin
    )

@pytest.fixture(scope='session')
def sampler_state(ad_vacuum_testsystem):
    states = pytest.importorskip('openmmtools.states')
    return states.SamplerState(positions=ad_vacuum_testsystem.positions)

@pytest.fixture(scope='session')
def minimized(sampler_state, thermodynamic_state):
    mcmc = pytest.importorskip('openmmtools.mcmc')
    move = mcmc.MCDisplacementMove()
    sampler = mcmc.MCMCSampler(thermodynamic_state, sampler_state, move)
    sampler.minimize()
    return sampler.sampler_state

@pytest.fixture
def snapshot(minimized, mcengine):
    return Snapshot.construct(
        coordinates=minimized.positions,
        velocities=None,
        box_vectors=None,
        engine=mcengine
    )

@pytest.fixture
def mcengine(thermodynamic_state, minimized):
    mcmc = pytest.importorskip('openmmtools.mcmc')
    move = mcmc.MCDisplacementMove()
    mcengine = OpenMMToolsMCEngine(
        thermodynamic_state,
        move,
        options={'n_steps_per_frame': 1,
                 'n_frames_max': 1000}
    ).named('mcengine')
    return mcengine

### HELPER CLASSES FOR MOCKING #############################################

if HAS_OPENMMTOOLS:
    class AlwaysMove(mcmc.MetropolizedMove):
        def __init__(self, accept):
            super(AlwaysMove, self).__init__()
            self.accept = {'accept': 1, True: 1,
                           'reject': 0, False: 0}[accept]

        def apply(self, thermodynamic_state, sampler_state):
            self.n_accepted += self.accept
            self.n_proposed += 1

        def _propose_positions(self, positions):
            pass

### TESTS FOR OPENMMTOOLS MC ENGINE ########################################

def test_snapshot_from_sampler_state(sampler_state):
    # the snapshot generated from snapshot_from_sampler_state should have
    # attributes that correspond to those of the sampler_state
    snapshot = snapshot_from_sampler_state(sampler_state, engine='foo')
    assert snapshot.velocities is None is sampler_state.velocities
    assert snapshot.coordinates is not sampler_state.positions
    assert (snapshot.coordinates == sampler_state.positions).all()


class TestOpenMMToolsMCEngine(object):
    def test_current_snapshot(self, snapshot):
        # when we set then get the current snapshot, we should get the same
        # object (identical in memory) back
        mcengine.current_snapshot = snapshot
        reloaded = mcengine.current_snapshot
        assert snapshot is reloaded

    def test_current_snapshot_set_same(self, mcengine, snapshot):
        # if we try to set the current snapshot the same snapshot we had
        # before, the snapshot remains the same object in memory
        assert mcengine._current_snapshot is None
        mcengine.current_snapshot = snapshot
        assert mcengine._current_snapshot is not None
        mcengine.current_snapshot = snapshot
        reloaded = mcengine.current_snapshot
        assert snapshot is reloaded

    def test_current_snapshot_update(self, mcengine, snapshot):
        # when we set a new snapshot (and one already exists) the
        # internal information should be updated
        unit = pytest.importorskip('simtk.unit')
        pos = snapshot.coordinates.value_in_unit(unit.nanometers)
        new_pos = np.random.random(pos.shape) * unit.nanometers
        new_snapshot = snapshot.construct(
            coordinates=new_pos,
            velocities=snapshot.velocities,
            box_vectors=snapshot.box_vectors,
            engine=snapshot.engine
        )
        mcengine.current_snapshot = snapshot
        mcengine.current_snapshot = new_snapshot
        sampler_state = mcengine._sampler.sampler_state
        assert not (sampler_state.positions == snapshot.coordinates).all()
        assert (sampler_state.positions == new_snapshot.coordinates).all()

    def test_current_snapshot_uninitialized_error(self, mcengine):
        # if we try to get the current snapshot before the engine has been
        # initialized (by setting the snapshot), we should raise an error
        with pytest.raises(EngineNotInitializedError):
            mcengine.current_snapshot

    def test_generate_next_frame(self, mcengine, minimized):
        # generate_next_frame should create a valid snapshot
        pytest.skip()

    def test_generate_next_frame_uninitialized_error(self, mcengine,
                                                     snapshot):
        # if we try to generate the next frame before setting the current
        # snapshot, we should error
        with pytest.raises(EngineNotInitializedError):
            mcengine.generate_next_frame()

    def test_generate_integration(self, mcengine, snapshot):
        # check that integration with the generate method: we should be able
        # to create a candidate trajectory for the ensemble
        ensemble = paths.LengthEnsemble(5)
        traj = mcengine.generate(snapshot, ensemble.can_append)
        assert isinstance(traj, paths.Trajectory)
        assert ensemble(traj)

    def test_serialization_cycle(self, mcengine):
        # if we serialize then deserialize, the to_dict of the deserialized
        # object should be the same as the to_dict of the original
        # breakpoint()
        serialized = mcengine.to_dict()
        deserialized = mcengine.__class__.from_dict(copy.copy(serialized))
        reserialized = deserialized.to_dict()
        assert serialized == reserialized

    def test_mdtraj_topology_error(self, mcengine):
        # if no MDTraj topology can be found, raise an error
        with pytest.raises(RuntimeError, match="no topology"):
            mcengine.mdtraj_topology

    def test_trajectory_to_mdtraj(self, snapshot, ad_vacuum_testsystem,
                                  thermodynamic_state):
        # if an MDTraj topology is associated with the engine, we
        # should be able to convert trajectories from this engine into
        # MDTraj trajectories
        md = pytest.importorskip('mdtraj')
        mdtraj_topology = ad_vacuum_testsystem.mdtraj_topology
        engine = OpenMMToolsMCEngine(
            thermodynamic_state,
            mcmc.MCDisplacementMove(),
            options={'n_steps_per_frame': 1,
                     'n_frames_max': 1000},
            topology=paths.engines.MDTrajTopology(mdtraj_topology)
        )
        new_snapshot = snapshot.copy_with_replacement(engine=engine)
        traj = paths.Trajectory([new_snapshot])
        mdt = traj.to_mdtraj()
        assert isinstance(mdt, md.Trajectory)
        assert len(mdt) == len(traj)

    @pytest.mark.parametrize('wrapper', [None, 'weighted', 'sequence'])
    @pytest.mark.parametrize('accept', ['accept', 'reject'])
    def test_get_n_accepted(self, wrapper, accept, thermodynamic_state,
                            snapshot):
        # (not API) ensure that we get the correct number of accepted steps
        # from a mover
        wrap = {
            None: lambda x: x,
            'weighted': lambda x: mcmc.WeightedMove([(x, 1)]),
            'sequence': lambda x: mcmc.SequenceMove([x]),
        }[wrapper]
        move = wrap(AlwaysMove(accept))
        results = {'accept': [1, 2], 'reject': [0, 0]}[accept]
        always_accept = OpenMMToolsMCEngine(
            move=move,
            thermodynamic_state=thermodynamic_state,
            options={'n_steps_per_frame': 1,
                     'n_frames_max': 1000},
        )
        always_accept.current_snapshot = snapshot
        assert always_accept._get_n_accepted() == 0
        for result in results:
            always_accept.generate_next_frame()
            assert always_accept._get_n_accepted() == result

    def test_generate_next_frame_accepted_move(self, thermodynamic_state,
                                               snapshot):
        # accepted moves should create new snapshots
        move = mcmc.WeightedMove([(AlwaysMove('accept'), 1)])
        always_accept = OpenMMToolsMCEngine(
            move=move,
            thermodynamic_state=thermodynamic_state,
            options={'n_steps_per_frame': 1,
                     'n_frames_max': 1000},
        )
        always_accept.current_snapshot = snapshot
        snap = always_accept.generate_next_frame()
        assert snap is not snapshot
        assert snap.__uuid__ != snapshot.__uuid__

    def test_generate_next_frame_rejected_move(self, thermodynamic_state,
                                               snapshot):
        # rejected moves should not create new snapshots
        move = mcmc.WeightedMove([(AlwaysMove('reject'), 1)])
        always_reject = OpenMMToolsMCEngine(
            move=move,
            thermodynamic_state=thermodynamic_state,
            options={'n_steps_per_frame': 1,
                     'n_frames_max': 1000},
        )
        always_reject.current_snapshot = snapshot
        snap = always_reject.generate_next_frame()
        assert snap is snapshot
