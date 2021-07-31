import pytest
import openpathsampling as paths

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
def snapshot(minimized):
    return Snapshot.construct(
        coordinates=minimized.positions,
        velocities=None,
        box_vectors=None
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
        assert snapshot.__uuid__ is reloaded.__uuid__
        # NOTE: snapshot.coordinates returns a copy to preserve immutability

    def test_current_snapshot_set_same(self, mcengine, snapshot):
        # if we try to set the current snapshot the same snapshot we had
        # before, the snapshot remains the same object in memory
        assert mcengine._current_snapshot is None
        mcengine.current_snapshot = snapshot
        assert mcengine._current_snapshot is not None
        mcengine.current_snapshot = snapshot
        reloaded = mcengine.current_snapshot
        assert snapshot is reloaded
        assert snapshot.__uuid__ is reloaded.__uuid__

    def test_current_snapshot_uninitialized_error(self, mcengine,
                                                  sampler_state):
        # if we try to get the current snapshot before the engine has been
        # initialized (by setting the snapshot), we should raise an error
        with pytest.raises(EngineNotInitializedError):
            mcengine.current_snapshot

    def test_generate_next_frame(self, mcengine, minimized):
        # generate_next_frame should create a valid snapshot
        pytest.skip()

    def test_serialization_cycle(self, mcengine):
        # if we serialize then deserialize, the to_dict of the deserialized
        # object should be the same as the to_dict of the original
        pytest.skip()

    def test_mdtraj_topology_error(self, mcengine):
        # if no MDTraj topology can be found, raise an error
        pytest.skip()

    def test_trajectory_to_mdtraj(self, ):
        # if an MDTraj topology is associated with the engine, we
        # should be able to convert trajectories from this engine into
        # MDTraj trajectories
        pytest.skip()

    def test_get_n_accepted(self):
        # (not API) ensure that we get the correct number of accepted steps
        # from a mover
        pytest.skip()

    def test_generate_next_frame_accepted_move(self):
        pytest.skip()

    def test_generate_next_frame_rejected_move(self):
        pytest.skip()
