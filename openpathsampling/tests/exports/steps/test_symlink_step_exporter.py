import pytest
import os
import pathlib
import openpathsampling as paths
from unittest.mock import Mock

from openpathsampling.exports.steps.symlink_step_exporter import (
    SymLinkStepExporter,
    _DEFAULT_TRIAL_PATTERN, _DEFAULT_ACTIVE_PATTERN,
    export_steps,
)

from openpathsampling.tests.analysis.utils.mock_movers import (
    MockForwardShooting, MockRepex,
    MockPathReversal, run_moves,
)

from openpathsampling.tests.test_helpers import make_1d_traj

@pytest.fixture
def shooting_step(default_unidirectional_tis):
    scheme = default_unidirectional_tis.scheme
    traj = default_unidirectional_tis.make_tis_trajectory(5)
    init_conds = scheme.initial_conditions_from_trajectories(traj)
    ensemble = scheme.network.sampling_ensembles[0]

    # Create a shooting move that will be accepted
    partial_traj = default_unidirectional_tis.make_trajectory(-1, 2).reversed
    move = MockForwardShooting(
        shooting_index=2,
        partial_traj=partial_traj,
        scheme=scheme,
        ensemble=ensemble,
        accepted=True
    )

    steplist = list(run_moves(init_conds, [move]))
    assert len(steplist) == 1
    step = steplist[0]
    assert ensemble(step.change.trials[0])  # acceptance is allowed
    return step

@pytest.fixture
def repex_step(default_unidirectional_tis):
    scheme = default_unidirectional_tis.scheme
    t1 = default_unidirectional_tis.make_tis_trajectory(4)
    t2 = default_unidirectional_tis.make_tis_trajectory(10)
    init_conds = scheme.initial_conditions_from_trajectories([t1, t2])

    e1, e2 = scheme.network.sampling_ensembles[:2]
    ensembles = [e1, e2]

    move = MockRepex(scheme, ensembles)

    steplist = list(run_moves(init_conds, [move]))
    assert len(steplist) == 1
    step = steplist[0]
    return step

@pytest.fixture
def pathreversal_step(default_unidirectional_tis):
    scheme = default_unidirectional_tis.scheme
    traj = default_unidirectional_tis.make_tis_trajectory(6)
    init_conds = scheme.initial_conditions_from_trajectories(traj)

    move = MockPathReversal(scheme)

    steplist = list(run_moves(init_conds, [move]))
    assert len(steplist) == 1
    step = steplist[0]
    return step

@pytest.fixture
def all_steps(shooting_step, repex_step, pathreversal_step):
    return [shooting_step, repex_step, pathreversal_step]

@pytest.fixture
def unnamed_ensemble():
    cv = paths.FunctionCV("x", lambda s: s.xyz[0][0])
    ensemble = paths.CVDefinedVolume(cv, -1.0, 1.0)
    return ensemble

@pytest.fixture
def named_ensemble():
    cv = paths.FunctionCV("x", lambda s: s.xyz[0][0])
    ensemble = paths.CVDefinedVolume(cv, -1.0, 1.0).named("EnsembleName")
    return ensemble


class TestSymLinkStepExporter:
    def setup_method(self):
        """Setup method to create exporter instance for tests."""
        self.mock_writer = Mock()
        self.mock_writer.ext = "dat"
        def mock_write_func(trajectory, filename):
            pathlib.Path(filename).parent.mkdir(parents=True, exist_ok=True)
            pathlib.Path(filename).touch()
        self.mock_writer.side_effect = mock_write_func
        self.exporter = SymLinkStepExporter(writer=self.mock_writer)

    @pytest.mark.parametrize("source", ["name", "uuid"])
    def test_get_ensemble_id(self, source, request):
        fixture = {'name': 'named_ensemble',
                   'uuid': 'unnamed_ensemble'}[source]
        ensemble = request.getfixturevalue(fixture)
        expected = {'name': 'EnsembleName',
                    'uuid': str(ensemble.__uuid__)}[source]

        sample = Mock()
        sample.ensemble = ensemble
        assert self.exporter._get_ensemble_id(sample) == expected

    @pytest.mark.parametrize("source", ["obj", "most_common"])
    def test_get_writer(self, source):
        if source == "obj":
            mock_writer = Mock()
            mock_writer.ext = "dat"
            exporter = SymLinkStepExporter(writer=mock_writer)
            sample = Mock()
            sample.trajectory = []

            result = exporter._get_writer(sample)
            assert result is mock_writer

        elif source == "most_common":
            exporter = SymLinkStepExporter(writer=None)

            traj1 = make_1d_traj([1.0])
            traj2 = make_1d_traj([2.0, 3.0])

            engine1 = traj1[0].engine
            engine2 = traj2[0].engine
            assert engine1 is not engine2

            combined_traj = traj1 + traj2

            sample = Mock()
            sample.trajectory = combined_traj

            result = exporter._get_writer(sample)
            default_writer = engine2._default_trajectory_writer()
            assert result.__class__ is default_writer.__class__

    def test_substitution_dict(self, shooting_step):
        exporter = SymLinkStepExporter(writer=self.mock_writer)

        step = shooting_step
        sample = step.change.trials[0]

        subs = exporter._substitution_dict(step, sample)

        assert subs["step"] is step
        assert subs["sample"] is sample
        assert subs["ensemble_id"] == exporter._get_ensemble_id(sample)
        assert subs["ext"] == self.mock_writer.ext

    def test_export_trial_sample(self, shooting_step, tmp_path):
        step = shooting_step
        sample = step.change.trials[0]

        exporter = SymLinkStepExporter(
            base_dir=tmp_path, writer=self.mock_writer
        )
        assert tmp_path.exists()
        assert tmp_path.is_dir()
        assert len(list(tmp_path.iterdir())) == 0

        original_cwd = os.getcwd()
        os.chdir(tmp_path)

        try:
            exporter.export_trial_sample(step, sample)

            subs_dict = exporter._substitution_dict(step, sample)
            raw_data_path = exporter.raw_data_pattern.format(**subs_dict)
            trial_path = exporter.trial_pattern.format(**subs_dict)
            assert pathlib.Path(raw_data_path).exists()

            assert pathlib.Path(trial_path).exists()
            assert pathlib.Path(trial_path).is_symlink()

            assert pathlib.Path(raw_data_path).samefile(trial_path)

        finally:
            os.chdir(original_cwd)

    def test_export_active_sample(self, shooting_step, tmp_path):
        step = shooting_step
        sample = step.active[0]

        exporter = SymLinkStepExporter(
            base_dir=tmp_path, writer=self.mock_writer
        )

        original_cwd = os.getcwd()
        os.chdir(tmp_path)

        try:
            exporter.export_active_sample(step, sample)

            subs_dict = exporter._substitution_dict(step, sample)
            raw_data_path = exporter.raw_data_pattern.format(**subs_dict)
            assert pathlib.Path(raw_data_path).exists()

            active_path = exporter.active_pattern.format(**subs_dict)
            assert pathlib.Path(active_path).exists()
            assert pathlib.Path(active_path).is_symlink()

            assert pathlib.Path(raw_data_path).samefile(active_path)

        finally:
            os.chdir(original_cwd)

    def test_export_raw_sample(self, shooting_step, tmp_path):
        step = shooting_step
        sample = step.active[0]

        exporter = SymLinkStepExporter(
            base_dir=tmp_path, writer=self.mock_writer
        )

        original_cwd = os.getcwd()
        os.chdir(tmp_path)

        try:
            exporter.export_raw_sample(step, sample)

            subs_dict = exporter._substitution_dict(step, sample)
            raw_data_path = exporter.raw_data_pattern.format(**subs_dict)
            assert pathlib.Path(raw_data_path).exists()

            assert pathlib.Path(raw_data_path).is_file()
            assert not pathlib.Path(raw_data_path).is_symlink()

        finally:
            os.chdir(original_cwd)

    @pytest.mark.parametrize("step_type", ["shooting", "repex",
                                           "pathreversal"])
    def test_export_step(self, step_type, request, tmp_path):
        expected_trials_by_type = {
            "shooting": 1,
            "repex": 2,
            "pathreversal": 1,
        }
        expected_trials = expected_trials_by_type[step_type]

        step = request.getfixturevalue(f"{step_type}_step")

        unique_trajectories = {
            s.trajectory for s in step.change.trials + list(step.active)
        }

        actual_trials = len(step.change.trials)
        expected_active = len(step.active)

        exporter = SymLinkStepExporter(
            base_dir=tmp_path, writer=self.mock_writer
        )

        original_cwd = os.getcwd()
        os.chdir(tmp_path)

        try:
            assert actual_trials == expected_trials

            exporter.export_step(step)

            raw_dir = pathlib.Path("raw_data")
            assert raw_dir.exists()
            raw_files = list(raw_dir.glob(f"*.{self.mock_writer.ext}"))
            assert len(raw_files) == len(unique_trajectories)

            trials_dir = pathlib.Path(f"{step.mccycle}/trials")
            assert trials_dir.exists()
            trial_files = list(trials_dir.glob(f"*.{self.mock_writer.ext}"))
            assert len(trial_files) == actual_trials

            for trial_file in trial_files:
                assert trial_file.is_symlink()

            active_dir = pathlib.Path(f"{step.mccycle}/active")
            assert active_dir.exists()
            active_files = list(
                active_dir.glob(f"*.{self.mock_writer.ext}")
            )

            assert len(active_files) == expected_active

            for active_file in active_files:
                assert active_file.is_symlink()

        finally:
            os.chdir(original_cwd)


def test_export_steps(all_steps, tmp_path):
    mock_writer = Mock()
    mock_writer.ext = "db"
    def mock_write_func(trajectory, filename):
        pathlib.Path(filename).parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(filename).touch()
    mock_writer.side_effect = mock_write_func

    original_cwd = os.getcwd()
    os.chdir(tmp_path)

    try:
        export_steps(all_steps, writer=mock_writer)

        assert pathlib.Path("raw_data").exists()
        raw_files = list(pathlib.Path("raw_data").glob(f"*.{mock_writer.ext}"))
        assert len(raw_files) > 0

        for step in all_steps:
            step_dir = pathlib.Path(str(step.mccycle))
            assert step_dir.exists()

            trials_dir = step_dir / "trials"
            assert trials_dir.exists()

            active_dir = step_dir / "active"
            assert active_dir.exists()

        writer_db2 = Mock()
        writer_db2.ext = "db2"
        writer_db2.side_effect = mock_write_func

        export_steps(all_steps, writer=writer_db2, export_trials=False)

        for step in all_steps:
            trials_dir = pathlib.Path(str(step.mccycle)) / "trials"
            assert trials_dir.exists()
            trial_files = list(trials_dir.glob(f"*.{writer_db2.ext}"))
            assert len(trial_files) == 0

        writer_db3 = Mock()
        writer_db3.ext = "db3"
        writer_db3.side_effect = mock_write_func

        export_steps(all_steps, writer=writer_db3, export_active=False)

        for step in all_steps:
            active_dir = pathlib.Path(str(step.mccycle)) / "active"
            assert active_dir.exists()
            active_files = list(active_dir.glob(f"*.{writer_db3.ext}"))
            assert len(active_files) == 0

    finally:
        os.chdir(original_cwd)
