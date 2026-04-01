_DEFAULT_RAW_DATA_PATTERN = "raw_data/{sample.trajectory.__uuid__}.{ext}"
_DEFAULT_TRIAL_PATTERN = "{step.mccycle}/trials/{ensemble_id}.{ext}"
_DEFAULT_ACTIVE_PATTERN = "{step.mccycle}/active/{ensemble_id}.{ext}"

import os
import collections
import pathlib

class SymLinkStepExporter:
    """Export steps as raw data and symlink to the raw data.

    In the patterns used for raw data, trials, and active data, the
    following substitutions are available:

    - step: The step being exported. This can be used with, e.g.,
      step.mccycle to get the Monte Carlo cycle number.
    - sample: The sample being exported.
    - ensemble_id: The ensemble ID for the sample. This is the name of the
      ensemble if it is set, or the UUID of the ensemble if the name is not
      set. This behavior can be customized by subclassing this class and
      overriding the :meth:`._get_ensemble_id` method.
    - ext: The extension for the raw data files.

    These parameters can be further customized by subclassing this class and
    overriding the :meth:`._substitution_dict` method.

    Parameters
    ----------
    writer : Optional[Callable]
       The TrajectoryWriter to use for exporting the raw data. If None, the
       default writer for the most common engine in the step will be used.
    raw_data_pattern : str
       File pattern for the raw data files.
    trial_pattern : str
       File pattern for the trial data symlinks.
    active_pattern : str
       File pattern for the active data symlinks.
    """
    def __init__(
        self,
        writer=None,
        *,
        base_dir='.',
        # TODO: add in a base using storage handlers, default to cwd
        raw_data_pattern=_DEFAULT_RAW_DATA_PATTERN,
        trial_pattern=_DEFAULT_TRIAL_PATTERN,
        active_pattern=_DEFAULT_ACTIVE_PATTERN,
    ):
        self.writer = writer
        self.base_dir = pathlib.Path(base_dir)
        self.raw_data_pattern = raw_data_pattern
        self.trial_pattern = trial_pattern
        self.active_pattern = active_pattern

    def _get_ensemble_id(self, sample):
        """Get the ensemble ID (used in file names) for a sample.
        """
        if sample.ensemble.is_named:
            ensemble_id = sample.ensemble.name
        else:
            ensemble_id = str(sample.ensemble.__uuid__)

        return ensemble_id

    def _get_writer(self, sample):
        """Get the TrajectoryWriter to use for a sample.

        If ``self.writer`` is set, it will be used. Otherwise, the most
        common engine in the sample's trajectory will be used.
        """
        if self.writer is not None:
            writer = self.writer
        else:
            engines = collections.Counter([s.engine for s in
                                          sample.trajectory])
            engine = engines.most_common(1)[0][0]
            writer = engine._default_trajectory_writer()
        return writer

    def _substitution_dict(self, step, sample):
        writer = self._get_writer(sample)
        ensemble_id = self._get_ensemble_id(sample)
        return {
            "step": step,
            "sample": sample,
            "ensemble_id": ensemble_id,
            "ext": writer.ext,
        }

    def _export_sample_symlink(self, pattern, step, sample):
        if pattern is None:
            return

        subs_dict = self._substitution_dict(step, sample)
        path = pattern.format(**subs_dict)
        raw_data_path = self.raw_data_pattern.format(**subs_dict)
        if not pathlib.Path(raw_data_path).exists():
            self.export_raw_sample(step, sample)

        pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
        symlink_path = pathlib.Path(path)
        target_path = pathlib.Path(raw_data_path)
        relative_target = os.path.relpath(target_path, symlink_path.parent)

        if not symlink_path.exists():
            os.symlink(relative_target, path)

    def export_trial_sample(self, step, sample):
        """Export a symlink to the raw data for a trial sample.

        Parameters
        ----------
        step : Step
           The step containing the sample.
       sample : Sample
           The trial sample to export.
        """
        self._export_sample_symlink(self.trial_pattern, step, sample)

    def export_active_sample(self, step, sample):
        """Export a symlink to the raw data for an active sample.

        Parameters
        ----------
        step : Step
           The step containing the sample.
        sample : Sample
           The active sample to export.
        """
        self._export_sample_symlink(self.active_pattern, step, sample)

    def export_raw_sample(self, step, sample):
        """Export the raw data for a sample.

        Parameters
        ----------
        step : Step
           The step containing the sample.
        sample : Sample
           The sample to export.
        """
        subs_dict = self._substitution_dict(step, sample)
        raw_data_path = self.raw_data_pattern.format(**subs_dict)
        if os.path.exists(raw_data_path):
            return

        # ensure parent directory exists
        pathlib.Path(raw_data_path).parent.mkdir(parents=True, exist_ok=True)
        writer = self._get_writer(sample)
        writer(sample.trajectory, raw_data_path)

    def export_step(self, step):
        """Export a step.

        Parameters
        ----------
        step : Step
           The step to export.
        """
        for sample in step.change.trials:
            self.export_raw_sample(step, sample)
            self.export_trial_sample(step, sample)

        for sample in step.active:
            self.export_active_sample(step, sample)


def export_steps(steps, writer=None, *, export_trials=True,
                 export_active=True):
    trial_pattern = _DEFAULT_TRIAL_PATTERN if export_trials else None
    active_pattern = _DEFAULT_ACTIVE_PATTERN if export_active else None
    exporter = SymLinkStepExporter(
        writer=writer,
        trial_pattern=trial_pattern,
        active_pattern=active_pattern,
    )
    for step in steps:
        exporter.export_step(step)
