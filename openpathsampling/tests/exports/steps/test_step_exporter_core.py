from openpathsampling.exports.steps.core import StepExporter


class FakeStepExporter(StepExporter):
    def __init__(self):
        self.calls = []

    def export_raw_sample(self, step, sample):
        self.calls.append(("raw", sample))

    def export_trial_sample(self, step, sample):
        self.calls.append(("trial", sample))

    def export_active_sample(self, step, sample):
        self.calls.append(("active", sample))


def test_step_exporter_export_step():
    class Change:
        trials = ["trial-1", "trial-2"]

    class Step:
        change = Change()
        active = ["active-1"]

    exporter = FakeStepExporter()
    exporter.export_step(Step())

    assert exporter.calls == [
        ("raw", "trial-1"),
        ("trial", "trial-1"),
        ("raw", "trial-2"),
        ("trial", "trial-2"),
        ("active", "active-1"),
    ]
