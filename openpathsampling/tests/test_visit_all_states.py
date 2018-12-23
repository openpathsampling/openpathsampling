import pytest

import openpathsampling as paths
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
        cv = paths.FunctionCV("x", lambda x: x.xyz[0][0])
        vol_A = paths.CVDefinedVolume(cv, 0.0, 1.0).named("A")
        vol_B = paths.CVDefinedVolume(cv, 2.0, 3.0).named("B")
        vol_C = paths.CVDefinedVolume(cv, 4.0, 5.0).named("C")
        vol_D = paths.CVDefinedVolume(cv, 6.0, 7.0).named("D")
        states  = [vol_A, vol_B, vol_C, vol_D]
        self.ensemble = VisitAllStatesEnsemble(states)
        # TODO: add trajectory

    def test_initialization(self):
        assert self.ensemble.progress == default_state_progress_report
        assert self.ensemble.report_frequency == 10
        assert self.ensemble.timestep == None

    def test_progress_indicator(self):
        assert self.ensemble._progress_indicator('default') == \
                default_state_progress_report
        assert self.ensemble._progress_indicator(None) is None

        def some_callable(n_steps, found_states, all_states, timestep):
            return 0.0

        assert self.ensemble._progress_indicator(some_callable) == \
                some_callable

    def test_state_for_frame(self):
        pass

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
