import openpathsampling as paths

import re
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

def extract_info_from_default_report(report):
    pattern = (r"Ran ([0-9]*) steps. Found states \[(.*)\]\. "
               + r"Looking for \[(.*)\]\.")
    result = re.match(pattern, report)
    match_groups = [result.group(i) for i in [1,2,3]]
    ret_val = [int(match_groups[0])]
    # extra complexity to handle the fact that "".split(',') gives [""]
    # (instead of [], as it does for "".split())
    # see https://stackoverflow.com/questions/16645083
    for state_list in match_groups[1:]:
        as_list = state_list.split(',')
        state_set = set(as_list) if as_list != [''] else {}
        ret_val.append(state_set)
    return ret_val

def test_extract_info_from_default_report():
    reports = {
        0: "Ran 0 steps. Found states []. Looking for [A,B,C,D].",
        3: "Ran 3 steps. Found states [A,B]. Looking for [C,D].",
        7: "Ran 7 steps. Found states [A,B,C,D]. Looking for []."
    }
    results = {
        0: [0, {}, {'A', 'B', 'C', 'D'}],
        3: [3, {'A', 'B'}, {'C', 'D'}],
        7: [7, {'A', 'B', 'C', 'D'}, {}]
    }
    for i in [0, 3, 7]:
        assert extract_info_from_default_report(reports[i]) == results[i]


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

        self.reports = [
            "Ran 0 steps. Found states []. Looking for [A,B,C,D].",
            "Ran 1 steps. Found states [A]. Looking for [B,C,D].",
            "Ran 2 steps. Found states [A]. Looking for [B,C,D].",
            "Ran 4 steps. Found states [A,B]. Looking for [C,D].",
            "Ran 5 steps. Found states [A,B,C]. Looking for [D].",
            "Ran 6 steps. Found states [A,B,C]. Looking for [D].",
            "Ran 7 steps. Found states [A,B,C,D]. Looking for [].",
        ]
        self.traj = make_1d_traj(sequence)
        # self.engine = CalvinistDynamics(self.traj)

    def _run_trajectory(self, can_append, trusted):
        n_frames = 1
        done = False
        my_traj = None
        while not done:
            my_traj = self.traj[:n_frames]
            done = not can_append(my_traj, trusted=trusted)
            n_frames += 1
        return my_traj

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

    @pytest.mark.parametrize('strict', [True, False],
                             ids=['strict', 'normal'])
    @pytest.mark.parametrize('trusted', [True, False],
                             ids=['trusted', 'untrusted'])
    def test_can_append(self, strict, trusted):
        can_append = {False: self.ensemble.can_append,
                      True: self.ensemble.strict_can_append}[strict]
        my_traj = self._run_trajectory(can_append, trusted)
        assert len(my_traj) == 8
        assert self.ensemble.found_states == set(self.states)

    @pytest.mark.parametrize('strict', [True, False],
                             ids=['strict', 'normal'])
    def test_can_append_trusted_incorrect(self, strict):
        # things that the trusted version can get wrong
        pytest.skip()

    @pytest.mark.parametrize('strict', [True, False],
                             ids=['strict', 'normal'])
    @pytest.mark.parametrize('trusted', [True, False],
                             ids=['trusted', 'untrusted'])
    def test_can_append_new_trajectory(self, strict, trusted):
        can_append = {False: self.ensemble.can_append,
                      True: self.ensemble.strict_can_append}[strict]
        my_traj = self._run_trajectory(can_append, trusted)
        report = self.ensemble.progress_report(my_traj)
        steps, found, looking = extract_info_from_default_report(report)
        assert steps == 7
        assert found == {'A','B','C','D'}
        assert looking == {}
        new_traj = make_1d_traj([-0.6])
        assert can_append(new_traj, trusted=trusted)
        report = self.ensemble.progress_report(new_traj)
        steps, found, looking = extract_info_from_default_report(report)
        assert steps == 0
        assert found == {}
        assert looking == {'A','B','C','D'}


    def test_call(self):
        pytest.skip()

    @pytest.mark.parametrize('strict', [True, False],
                             ids=['strict', 'normal'])
    def test_can_prepend(self, strict):
        can_prepend = {False: self.ensemble.can_prepend,
                       True: self.ensemble.strict_can_prepend}[strict]
        with pytest.raises(NotImplementedError):
            can_prepend(self.traj[0])
