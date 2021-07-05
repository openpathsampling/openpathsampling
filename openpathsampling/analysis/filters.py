import openpathsampling as paths
import enum

def _default_repr(obj):
    return obj.name if obj.name is not None \
            else super(obj.__class__, obj).__repr__()

class _Filter(object):
    """Abstract base class for filters.
    """
    FILTER_TYPE = None
    def condition(self, inp):
        raise NotImplementedError()

    def __and__(self, other):
        # need fancier logic here to combine multiple stages
        return AndFilterCombination(self, other)

    def __or__(self, other):
        return OrFilterCombination(self, other)

    def __invert__(self):
        return NegationFilter(self)

    def __sub__(self, other):
        return self & ~other

    def __xor__(self, other):
        return (self & ~other) | (~self & other)

    def __call__(self, inputs):
        return filter(self.condition, inputs)


class NegationFilter(_Filter):
    @property
    def FILTER_TYPE(self):
        return self.filter.FILTER_TYPE

    def __init__(self, filt):
        self.filter = filt

    def condition(self, inp):
        return not self.filter(inp)

    def __repr__(self):
        return "~" + repr(self.filter)

class _FilterCombination(_Filter):
    @property
    def FILTER_TYPE(self):
        retrun self.filter1.FILTER_TYPE

    def __init__(self, filter1, filter2):
        if filter1.FILTER_TYPE is not filter2.FILTER_TYPE:
            raise TypeError(f"Cannot combine filters {filter1} and {filter2}")
        self.filter1 = filter1
        self.filter2 = filter2


class AndFilterCombination(_FilterCombination):
    def condition(self, inp):
        return self.filter1.condition(inp) and self.filter2.condition(inp)

    def __repr__(self):
        return "(" + repr(self.filter1) + " & " + repr(self.filter2) + ")"


class OrFilterCombination(_FilterCombination):
    def condition(self, inp):
        return self.filter1.condition(inp) or self.filter2.condition(inp)

    def __repr__(self):
        return "(" + repr(self.filter1) + " | " + repr(self.filter2) + ")"


class GenericFilter(_Filter):
    def __init__(self, condition, name):
        self._condition = condition
        self.name = name

    def condition(self, inp):
        return self._condition(inp)

    def __repr__(self):
        return _default_repr(self)

### OPS STEP FILTERS #######################################################

class CanonicalMover(_Filter):
    FILTER_TYPE = paths.MCStep
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
        return f"(canonical mover is {self.mover})"


class TrialReplica(_Filter):
    FILTER_TYPE = paths.MCStep
    def __init__(self, replica):
        self.replica = replica

    def condition(self, step):
        replicas = [s.replica for s in step.change.canonical.trials]
        return self.replica in replicas


class TrialEnsemble(_Filter):
    FILTER_TYPE = paths.MCStep
    def __init__(self, ensemble):
        self.ensemble = ensemble

    def condition(self, step):
        ensembles = [s.ensemble for s in step.change.canonical.trials]
        return self.ensemble in ensembles


class StepFilter(GenericFilter):
    FILTER_TYPE = paths.MCStep

RejectedSteps = StepFilter(lambda step: not step.change.accepted,
                           name="RejectedSteps")
AcceptedSteps = StepFilter(lambda step: step.change.accepted,
                           name="AcceptedSteps")
AllSteps = StepFilter(lambda step: True,
                      name="AllSteps")


### SAMPLE FILTERS #########################################################

class SampleFilterCondition(GenericFilter):
    FILTER_TYPE = paths.Sample

AllSamples = SampleFilterCondition(lambda sample: True, name="AllSamples")

class Ensemble(_Filter):
    FILTER_TYPE = paths.Sample
    def __init__(self, ensemble):
        self.ensemble = ensemble

    def condition(self, sample):
        return sample.ensemble == self.ensemble

    def __repr__(self):
        return f"Ensemble(<Ensemble name='{self.ensemble.name}'>)"


class Replica(_Filter):
    FILTER_TYPE = paths.Sample
    def __init__(self, replica):
        self.replica = replica

    def condition(self, sample):
        return sample.replica == self.replica

    def __repr__(self):
        return f"Replica({self.replica})"


class _NetworkEnsemble(_Filter):
    FILTER_TYPE = paths.Sample
    def __init__(self, network):
        super().__init__()
        self.network = network
        self.ensembles = set(self._get_ensembles())

    def _get_ensembles(self):
        raise NotImplementedError()

    def condition(self, sample):
        return sample.ensemble in self.ensembles

    def __repr__(self):
        return f"{self.__class__.__name__}({self.network})"


class SamplingEnsemble(_NetworkEnsemble):
    def _get_ensembles(self):
        return self.network.sampling_ensembles


class MinusEnsemble(_NetworkEnsemble):
    def _get_ensembles(self):
        return self.network.minus_ensembles


class MSOuterEnsemble(_NetworkEnsemble):
    def _get_ensembles(self):
        return self.network.ms_outers



### EXTRACTOR-FILTERS ######################################################

NOT_EXTRACTED = object()  # flag for errors on extraction

def _flatten(items):
    for item in items:
        if item is not NOT_EXTRACTED:
            yield item

def _list_per_step(items):
    return [list(items)]

class ExtractorFilter(object):
    _default_secondary = None
    def __init__(self, step_filter, extractor):
        if step_filter is None:
            step_filter = AllSteps
        # if postprocessor is None:
            # postprocessor = _default_postprocessor

        self.step_filter = step_filter
        self.extractor = extractor
        # self.postprocessor = postprocessor
        self.secondary_filter = self._default_secondary

    def _get_finalizer(self, flatten):
        if flatten is None:
            flatten = True  # default behavior is True
        finalizer = {True: _flatten,
                     False: _list_per_step}[flatten]
        return finalizer

    def with_filter(self, secondary_filter):
        if self.secondary_filter != self._default_secondary:
            raise RuntimeError("A secondary filter has already been set.")
        self.secondary_filter = secondary_filter
        return self

    def __call__(self, steps, flatten=None):
        finalizer = self._get_finalizer(flatten)
        for step in self.step_filter(steps):
            extracted = self.extractor(step)

            if not self.extractor.returns_iterable:
                extracted = [extracted]

            if self.secondary_filter is not None:
                extracted = self.secondary_filter(extracted)
            # extracted = self.postprocessor(extracted)
            finalized = finalizer(extracted)
            for result in finalized:
                yield result

# For the time being, leaving all this combination-related code commented
# out. It's enough typing to keep it around in case we decide we want it,
# but I have the feeling right now that it would be overly complicated from
# a user standpoint.
#
#    def __and__(self, other):
#        return _AndExtractorFilter(self, other)
#
#    def __or__(self, other):
#        return _OrExtractorFilter(self, other)
#
#    def __invert__(self):
#        return _NegationExtractorFilter(self)
#
#    def _test_compatibility(self, other):
#        value_err_str = ("Can not combine ExtractorFilters that do not "
#                         "have the same {attr}: {self_attr} != "
#                         "{other_attr}")
#        if self.extractor != other.extractor:
#            raise ValueError(value_err_str.format(
#                attr='extractor',
#                self_attr=self.extractor,
#                other_attr=other.extractor
#            ))
#        if self.secondary_filter != other.secondary_filter:
#            raise ValueError(value_err_str.format(
#                attr='secondary filter',
#                self_attr=self.secondary_filter,
#                other_attr=other.secondary_filter
#            ))
#
#class _AndExtractorFilter(ExtractorFilter):
#    def __init__(self, ef1, ef2):
#        ef1._test_compatibility(ef2)
#        step_filter = ef1.step_filter & ef2.step_filter
#        self._default_secondary = ef1._default_secondary
#        super(_AndExtractorFilter, self).__init__(step_filter,
#                                                  ef1.extractor)
#
#class _OrExtractorFilter(ExtractorFilter):
#    def __init__(self, ef1, ef2):
#        ef1._test_compatibility(ef2)
#        step_filter = ef1.step_filter | ef2.step_filter
#        self._default_secondary = ef1._default_secondary
#        super(_OrExtractorFilter, self).__init__(step_filter,
#                                                  ef1.extractor)
#
#class _NegationExtractorFilter(ExtractorFilter):
#    def __init__(self, ef):
#        self._default_secondary = ef._default_secondary
#        super(_NegationExtractorFilter, self).__init__(~ef.step_filter,
#                                                       ef.extractor)
#        self.with_filter(ef.secondary_filter)


class SampleExtractorFilter(ExtractorFilter):
    _default_secondary = AllSamples


### EXTRACTORS #############################################################

class Extractor(object):
    def __init__(self, extractor, name=None,
                 returns_iterable=False,
                 extractor_filter=ExtractorFilter):
        self.extractor = extractor
        self.name = name
        self.returns_iterable = returns_iterable
        self.extractor_filter = extractor_filter

    def __call__(self, step):
        return self.extractor(step)

    def using(self, step_filter=None):
        return self.extractor_filter(step_filter=step_filter,
                                     extractor=self)

    def __repr__(self):
        return _default_repr(self)


class CanonicalDetails(Extractor):
    def __init__(self, detail_name, returns_iterable=False):
        super(CanonicalDetails, self).__init__(
            extractor=lambda step: getattr(step.change.canonical.details,
                                           detail_name,
                                           None),
            name="CanonicalDetails({detail_name})".format(detail_name),
            returns_iterable=returns_iterable
        )


class SampleExtractor(Extractor):
    def __init__(self, extractor, name):
        super(SampleExtractor, self).__init__(
            extractor=extractor,
            name=name,
            returns_iterable=True,
            extractor_filter=SampleExtractorFilter
        )

ActiveSamples = SampleExtractor(lambda step: step.active.samples,
                                name="ActiveSamples")
TrialSamples = SampleExtractor(lambda step: step.change.canonical.trials,
                               name="TrialSamples")


def _get_shooting_point(step):
    details = step.change.canonical.details
    try:
        shooting_pt = details.shooting_snapshot
    except AttributeError:
        shooting_pt = NOT_EXTRACTED
    return shooting_pt

ShootingSteps = StepFilter(
   lambda step: _get_shooting_point(step) is not NOT_EXTRACTED,
   name="ShootingSteps"
)

ShootingPoints = Extractor(_get_shooting_point, name="ShootingPoints")

def _get_modified_shooting_point(step):
    if not ShootingSteps.condition(step):
        return NOT_EXTRACTED
    details = step.change.canonical.details
    try:
        shooting_pt = ...
    except AttributeError:
        # return the original shooting point for one-way
        shooting_pt = details.shooting_snapshot
    return shooting_pt

ModifiedShootingPoints = Extractor(_get_modified_shooting_point,
                                   name="ModifiedShootingPoints")
