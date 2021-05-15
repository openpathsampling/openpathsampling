import openpathsampling as paths
import enum

class Stages(enum.Enum):
    UNKNOWN = 0
    STEP_FILTER = 1
    SAMPLE_SELECTOR = 2
    SAMPLE_FILTER = 3
    POSTPROCESSING = 4

class _Filter(object):
    """Abstract base class for filters.
    """
    STAGE = Stages.UNKNOWN
    def condition(self, inp):
        raise NotImplementedError()

    def __and__(self, other):
        # need fancier logic here to combine multiple stages
        return AndFilterCombination(self, other)

    def __or__(self, other):
        return OrFilterCombination(self, other)

    def __invert__(self):
        return NegationFilter(self)

    def __call__(self, inputs):
        return filter(self.condition, inputs)


class NegationFilter(_Filter):
    @property
    def STAGE(self):
        return self.filter.STAGE

    def __init__(self, filt):
        self.filter = filt

    def condition(self, inp):
        return not self.filter(inp)


class AndFilterCombination(_Filter):
    @property
    def STAGE(self):
        return self.filter1.STAGE

    def __init__(self, filter1, filter2):
        if filter1.STAGE != filter2.STAGE:
            raise TypeError(f"Cannot combine filters {filter1} and {filter2}")
        self.filter1 = filter1
        self.filter2 = filter2

    def condition(self, inp):
        return self.filter1.condition(inp) and self.filter2.condition(inp)


class OrFilterCombination(_Filter):
    @property
    def STAGE(self):
        return self.filter1.STAGE

    def __init__(self, filter1, filter2):
        if filter1.STAGE != filter2.STAGE:
            raise TypeError(f"Cannot combine filters {filter1} and {filter2}")
        self.filter1 = filter1
        self.filter2 = filter2

    def condition(self, inp):
        return self.filter1.condition(inp) or self.filter2.condition(inp)


class CanonicalMover(_Filter):
    STAGE = Stages.STEP_FILTER
    def __init__(self, mover):
        if isinstance(mover, paths.PathMover):
            self.mover = mover
            self._condition = self.equality_check
        elif isinstance(mover, str):
            self.mover = getattr(paths, mover)
            self._condition = self.isinstance_check
        elif issubclass(mover, paths.PathMover):
            self.mover = mover
            self._condition = self.isinstance_check
        else:
            raise ValueError(f"{mover} does not appear to be a path mover")

    def equality_check(self, mover):
        return mover == self.mover

    def isinstance_check(self, mover):
        return isinstance(mover, self.mover)

    def condition(self, step):
        return self._condition(step.change.canonical.mover)

    def __str__(self):
        return f"(canonical mover is {self.mover}"


class StepFilter(_Filter):
    STAGE = Stages.STEP_FILTER
    def __init__(self, condition):
        self.condition = condition

    def condition(self, step):
        return self.condition(step)

RejectedSteps = StepFilter(lambda step: not step.change.accepted)
AcceptedSteps = StepFilter(lambda step: step.change.accepted)

class SampleSelector(_Filter):
    STAGE = Stages.SAMPLE_SELECTOR
    def __init__(self, selector):
        self.selector = selector

    def condition(self, step):
        return self.selector(step)

    def __call__(self, step):
        return self.selector(step)

ActiveSamples = SampleSelector(lambda step: step.active.samples)
TrialSamples = SampleSelector(lambda step: step.change.canonical.trials)


class Ensemble(_Filter):
    STAGE = Stages.SAMPLE_FILTER
    def __init__(self, ensemble):
        self.ensemble = ensemble

    def condition(self, sample):
        return sample.ensemble == self.ensemble


class Replica(_Filter):
    STAGE = Stages.SAMPLE_FILTER
    def __init__(self, replica):
        self.replica = replica

    def condition(self, sample):
        return sample.replica == self.replica

class PostProcess(_Filter):
    STAGE = Stages.POSTPROCESSING
    def __init__(self, method):
        self.method = method

    def condition(self, samples):
        for item in self.method(samples):
            yield item

def _flatten(samples):
    for sample in samples:
        yield sample

def _list_per_step(samples):
    return [list(samples)]

Flatten = PostProcess(_flatten)
ListPerStep = PostProcess(_list_per_step)

class SampleFilter(object):
    def __init__(self, step_filter=None, sample_selector=None,
                 sample_filter=None, postprocess_filter=None):
        self.step_filter = step_filter
        self.sample_selector = sample_selector
        self.sample_filter = sample_filter
        if postprocess_filter is None:
            postprocess_filter = ListPerStep
        self.postprocess_filter = postprocess_filter

    def __call__(self, steps, flatten=None):
        postprocess = None
        if flatten is True:
            postprocess = Flatten
        elif flatten is False:
            postprocess = ListPerStep
        elif flatten is None:
            postprocess = self.postprocess_filter

        if postprocess is None:
            raise RuntimeError()

        for step in self.step_filter(steps):
            samples = self.sample_selector(step)
            filtered = self.sample_filter(samples)
            finalized = postprocess(filtered)
            for result in finalized:
                yield result

# once we can mix stages:
# TPSActiveSamples = ActiveSamples & Flatten
# RejectedTrials = RejectedSteps & TrialSamples
