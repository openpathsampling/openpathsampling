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

    def __repr__(self):
        return "~" + repr(self.filter)


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

    def __repr__(self):
        return "(" + repr(self.filter1) + " & " repr(self.filter2) + ")"


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

    def __repr__(self):
        return "(" + repr(self.filter1) + " | " repr(self.filter2) + ")"


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


class StepFilterCondition(_Filter):
    STAGE = Stages.STEP_FILTER
    def __init__(self, condition, name):
        self._condition = condition
        self.name = name

    def condition(self, step):
        return self._condition(step)

    def __repr__(self):
        return self.name

RejectedSteps = StepFilterCondition(lambda step: not step.change.accepted,
                                    name="RejectedSteps")
AcceptedSteps = StepFilterCondition(lambda step: step.change.accepted,
                                    name="AcceptedSteps")
AllSteps = StepFilterCondition(lambda step: True,
                               name="AllSteps")


def _pass_thru(items):
    for item in items:
        yield item

class SampleFilterCondition(_Filter):
    STAGE = Stages.SAMPLE_FILTER
    def __init__(self, condition, name):
        self.condition = condition
        self.name = name

    def __repr__(self):
        return self.name

AllSamples = SampleFilterCondition(lambda sample: True, name="AllSamples")

class Ensemble(_Filter):
    STAGE = Stages.SAMPLE_FILTER
    def __init__(self, ensemble):
        self.ensemble = ensemble

    def condition(self, sample):
        return sample.ensemble == self.ensemble

    def __repr__(self):
        return f"Ensemble(<Ensemble name='{self.ensemble.name}'>)"


class Replica(_Filter):
    STAGE = Stages.SAMPLE_FILTER
    def __init__(self, replica):
        self.replica = replica

    def condition(self, sample):
        return sample.replica == self.replica

    def __repr__(self):
        return f"Replica({self.replica})"

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

class ExtractorFilter(object):
    def __init__(self, step_filter, extractor, sample_filter, hook=None):
        if step_filter is None:
            step_filter = AllSteps
        if sample_filter is None:
            sample_filter = AllSamples
        if hook is None:
            hook = _pass_thru

        self.step_filter = step_filter
        self.extractor = extractor
        self.sample_filter = sample_filter
        self.hook = hook

    def _get_postprocess(self, flatten):
        postprocess = {True: Flatten,
                       False: ListPerStep,
                       None: ListPerStep}[flatten]
        return postprocess

    def __call__(self, steps, flatten=None):
        raise NotImplementedError()


class SampleFilter(ExtractorFilter):
    def __call__(self, steps, flatten=None):
        postprocess = self._get_postprocess(flatten)
        for step in self.step_filter(steps):
            samples = self.extractor(step)
            filtered = self.sample_filter(samples)
            processed = self.hook(filtered)
            # TODO: self.extractor_hook here?
            finalized = postprocess(filtered)
            for result in finalized:
                yield result


class Extractor(object):
    STAGE = Stages.SAMPLE_SELECTOR
    hook = _pass_thru
    def __init__(self, extractor, name, extract_filter):
        self.extractor = extractor
        self.name = name
        self.extract_filter = extract_filter

    def __call__(self, step):
        return self.extractor(step)

    def __matmul__(self, other):
        if len(other) != 2:
            # TODO: better error
            raise ValueError("Requires a tuple of 2 items")
        # TODO: test that the inputs are the correct type
        step_filter, sample_filter = other
        return self.using(step_filter, sample_filter)

    def using(self, step_filter=None, sample_filter=None):
        return self.extract_filter(step_filter=step_filter,
                                   extractor=self,
                                   sample_filter=sample_filter)

    def __repr__(self):
        return self.name

class SampleExtractor(Extractor):
    def __init__(self, extractor, name):
        super(SampleExtractor, self).__init__(extractor=extractor,
                                              name=name,
                                              extract_filter=SampleFilter)

class ActiveEnsembles(Extractor):
    def __init__(self, ensemble):
        super().__init__(extractor=lambda step: step.active_ensembles,
                         name=(f"<Ensemble named={ensemble.name} "
                               f"uuid={ensemble.__uuid__}>"))
        self.ensemble = ensemble

    def hook(self, samples):
        return [paths.SampleSet(samples)[self.ensemble]]




ActiveSamples = SampleExtractor(lambda step: step.active.samples,
                               name="ActiveSamples")
TrialSamples = SampleExtractor(lambda step: step.change.canonical.trials,
                              name="TrialSamples")




def _get_shooting_point(step):
    details = step.change.canonical.details
    try:
        shooting_pt = details.shooting_snapshot
    except AttributeError:
        shooting_pt = None
    return shooting_pt

ShootingSteps = StepFilterCondition(
   lambda step: _get_shooting_point(step) is not None,
   name="ShootingSteps"
)

ShootingPoints = Extractor(_get_shooting_point, name="ShootingPoints",
                           extract_filter=None)

class _ShootingPointCondition(_Filter):
   STAGE = Stages.SAMPLE_FILTER
   def __init__(self, shooting_point):
       super().__init__()
       self.shooting_point = shooting_point

   def condition(self, sample):
       return self.shooting_point in sample.trajectory
        

# TODO: can this ExtractFilter just be generalized?
class ShootingPointFilter(ExtractorFilter):
   PayloadFilter = _ShootingPointCondition
   def __call__(self, steps, flatten=None):
       postprocess = self._get_postprocess(flatten)
       for step in self.step_filter(steps):
           payload = self.extractor(step)
           samples = TrialSamples(step)
           sample_filter = self.sample_filter & self.PayloadFilter(payload)
           filtered = self.sample_filter(samples)


# once we can mix stages:
# TPSActiveSamples = ActiveSamples & Flatten
# RejectedTrials = RejectedSteps & TrialSamples
