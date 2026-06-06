from abc import ABC, abstractmethod


class StepExporter(ABC):
    """Abstract interface for exporting simulation steps.

    Implementations define how raw trajectory data and trial/active
    references are represented. The default orchestration assumes
    single-process use and does not provide locking around duplicate exports.
    """

    @abstractmethod
    def export_raw_sample(self, step, sample):
        """Export the raw trajectory data for ``sample``."""
        raise NotImplementedError()

    @abstractmethod
    def export_trial_sample(self, step, sample):
        """Export the trial reference/data for ``sample``."""
        raise NotImplementedError()

    @abstractmethod
    def export_active_sample(self, step, sample):
        """Export the active reference/data for ``sample``."""
        raise NotImplementedError()

    def export_step(self, step):
        """Export all trial and active samples from ``step``."""
        for sample in step.change.trials:
            self.export_raw_sample(step, sample)
            self.export_trial_sample(step, sample)

        for sample in step.active:
            self.export_active_sample(step, sample)
