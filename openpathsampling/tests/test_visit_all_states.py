import openpathsampling as paths

import pytest
from openpathsampling.tests.test_helpers import \
        make_1d_traj, CalvinistDynamics


from openpathsampling.visit_all_states import *

def test_default_state_progress_report():
    cv = paths.FunctionCV("x", lambda x: x.xyz[0][0])
    vol_A = paths.CVDefinedVolume(cv, 0.0, 1.0).named("A")
    vol_B = paths.CVDefinedVolume(cv, 2.0, 3.0).named("B")
    vol_C = paths.CVDefinedVolume(cv, 4.0, 5.0).named("C")
    vol_D = paths.CVDefinedVolume(cv, 6.0, 7.0).named("D")

    n_steps = 100
    found_vol = [vol_A, vol_B]
    all_vol  = [vol_A, vol_B, vol_C, vol_D]
    tstep = 0.5

    f = default_state_progress_report  # keep on one line
    assert f(n_steps, found_vol, all_vol) == \
            "Ran 100 steps. Found states [A,B]. Looking for [C,D]."
    assert f(n_steps, found_vol, all_vol, tstep) == \
            "Ran 100 steps [50.0]. Found states [A,B]. Looking for [C,D]."


class TestVisitAllStatesEnsemble(object):
    def setup(self):
        self.cv = paths.FunctionCV("x", lambda x: x.xyz[0][0])
        vol_A = paths.CVDefinedVolume(self.cv, 0.0, 1.0).named("A")
        vol_B = paths.CVDefinedVolume(self.cv, 2.0, 3.0).named("B")
        vol_C = paths.CVDefinedVolume(self.cv, 4.0, 5.0).named("C")
        vol_D = paths.CVDefinedVolume(self.cv, 6.0, 7.0).named("D")
        self.states  = [vol_A, vol_B, vol_C, vol_D]
        self.ensemble = VisitAllStatesEnsemble(self.states)
        sequence = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]
        self.state_seq = [[], [vol_A], [], [vol_B], [], [vol_C], [],
                          [vol_D], []]
        self.traj = make_1d_traj(sequence)
        self.engine = CalvinistDynamics(self.traj)

    def test_initialization(self):
        assert self.ensemble.progress == default_state_progress_report
        assert self.ensemble.report_frequency == 10
        assert self.ensemble.timestep == None

    def test_progress_indicator(self):
        assert self.ensemble._progress_indicator('default') == \
                default_state_progress_report
        assert self.ensemble._progress_indicator(None) is None

        def some_callable(n_steps, found_states, all_states, timestep):
            return "some_callable"

        assert self.ensemble._progress_indicator(some_callable) == \
                some_callable

    def test_state_for_frame(self):
        for snap, expected in zip(self.traj, self.state_seq):
            assert self.ensemble._state_for_frame(snap) == set(expected)

    def test_state_for_frame_error(self):
        vol_A = self.states[0]
        vol_A_prime = paths.CVDefinedVolume(self.cv, 0.25, 0.75).named("A'")
        snap = self.traj[1]
        assert vol_A(snap)
        assert vol_A_prime(snap)

        ensemble = VisitAllStatesEnsemble(states=[vol_A, vol_A_prime])
        with pytest.raises(RuntimeError):
            ensemble._state_for_frame(snap)

    def test_can_append(self):
        pass

    def test_strict_can_append(self):
        pass

    def test_call(self):
        pass

    def test_can_prepend(self):
        pass

    def test_strict_can_prepend(self):
        pass
