import openpathsampling as paths
import enum
import inspect

def _default_repr(obj):  # pragma: no cover
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
        return not self.filter.condition(inp)

    def __repr__(self):
        return "~" + repr(self.filter)

class _FilterCombination(_Filter):
    @property
    def FILTER_TYPE(self):
        return self.filter1.FILTER_TYPE

    def __init__(self, filter1, filter2):
        if filter1.FILTER_TYPE is not filter2.FILTER_TYPE:
            raise TypeError("Cannot combine filters " + str(filter1)
                            + " and " + str(filter2))
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
    """A generic version of the user-facing code for filters.

    Filters for specific types are subclasses that set the ``FILTER_TYPE``
    class variable.

    Parameters
    ----------
    condition : Callable[Any] -> bool
        a function that takes an instance of the input type of this filter
        and returns True if that instance should be kept, or False if it
        should be filtered out
    name : str
        name string, primarily for use in error reporting
    """
    def __init__(self, condition, name):
        self._condition = condition
        self.name = name

    def condition(self, inp):
        return self._condition(inp)

    def __repr__(self):
        return _default_repr(self)

### OPS STEP FILTERS #######################################################

class canonical_mover(_Filter):
    FILTER_TYPE = paths.MCStep
    def __init__(self, mover):
        if isinstance(mover, paths.PathMover):
            self.mover = mover
            self._condition = self.equality_check
        elif isinstance(mover, str):
            self.mover = getattr(paths, mover)
            self._condition = self.isinstance_check
        elif inspect.isclass(mover) and issubclass(mover, paths.PathMover):
            self.mover = mover
            self._condition = self.isinstance_check
        elif isinstance(mover, paths.MCStep):
            raise TypeError("canonical_mover filter initialized with a "
                            "step: Did you mean to use the extractor, "
                            "canonical_movers? (note the final 's')")
        else:
            raise TypeError(str(mover) + " does not appear to be a path "
                            "mover")

    def equality_check(self, mover):
        return mover == self.mover

    def isinstance_check(self, mover):
        return isinstance(mover, self.mover)

    def condition(self, step):
        return self._condition(step.change.canonical.mover)

    def __str__(self):
        return "(canonical mover is " + str(self.mover) + ")"


class trial_replica(_Filter):
    FILTER_TYPE = paths.MCStep
    def __init__(self, replica):
        self.replica = replica

    def condition(self, step):
        replicas = [s.replica for s in step.change.canonical.trials]
        return self.replica in replicas


class trial_ensemble(_Filter):
    FILTER_TYPE = paths.MCStep
    def __init__(self, ensemble):
        self.ensemble = ensemble

    def condition(self, step):
        ensembles = [s.ensemble for s in step.change.canonical.trials]
        return self.ensemble in ensembles


class StepFilter(GenericFilter):
    FILTER_TYPE = paths.MCStep

rejected_steps = StepFilter(lambda step: not step.change.accepted,
                            name="rejected_steps")
accepted_steps = StepFilter(lambda step: step.change.accepted,
                            name="accepted_steps")
all_steps = StepFilter(lambda step: True,
                       name="all_steps")


### SAMPLE FILTERS #########################################################

class SampleFilterCondition(GenericFilter):
    FILTER_TYPE = paths.Sample

all_samples = SampleFilterCondition(lambda sample: True, name="all_samples")

class ensemble(_Filter):
    FILTER_TYPE = paths.Sample
    def __init__(self, ensemble):
        self.ensemble = ensemble

    def condition(self, sample):
        return sample.ensemble == self.ensemble

    def __repr__(self):
        return "Ensemble(<Ensemble name='" + self.ensemble.name + "'>)"


class replica(_Filter):
    FILTER_TYPE = paths.Sample
    def __init__(self, replica):
        self.replica = replica

    def condition(self, sample):
        return sample.replica == self.replica

    def __repr__(self):
        return "Replica(" + str(self.replica) + ")"


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
        return str(self.__class__.__name__) + "(" + repr(self.network) + ")"


class sampling_ensemble(_NetworkEnsemble):
    def _get_ensembles(self):
        return self.network.sampling_ensembles


class minus_ensemble(_NetworkEnsemble):
    def _get_ensembles(self):
        return self.network.minus_ensembles


class ms_outer_ensemble(_NetworkEnsemble):
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
            step_filter = all_steps
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

class SampleExtractorFilter(ExtractorFilter):
    _default_secondary = all_samples


### EXTRACTORS #############################################################

class Extractor(object):
    """Extract information from a step.

    Parameters
    ----------
    extractor : Callable[[MCStep], Any]
        function that extracts the desired information from an MCStep
    name : Optional[str]
        name of this extractor (used in error messages)
    returns_iterable : bool
        whether the return object is an iterable
    extractor_filter : class
        the class to use to create the extractor filter from the
        meth:`.using` method.
    """
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
        """Create an :class:`.ExtractorFilter` with the given filter

        Parameters
        ----------
        step_filter : :class:`.StepFilter`
            the filter to use for in the resulting extractor-filter

        Returns
        -------
        ExtractorFilter :
            object that filters steps and extracts the desired information
        """
        return self.extractor_filter(step_filter=step_filter,
                                     extractor=self)

    def __repr__(self):
        return _default_repr(self)


class canonical_details(Extractor):
    def __init__(self, detail_name, returns_iterable=False):
        super(CanonicalDetails, self).__init__(
            extractor=lambda step: getattr(step.change.canonical.details,
                                           detail_name,
                                           None),
            name="canonical_details({detail_name})".format(detail_name),
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

active_samples = SampleExtractor(lambda step: step.active.samples,
                                 name="active_samples")
"""Extractor for active samples

Active samples count toward statistical sampling, and appear repeatedly to
increase their (Monte Carlo) weight.
"""

trial_samples = SampleExtractor(lambda step: step.change.canonical.trials,
                                name="trial_samples")
"""Extractor for trial samples"""


def _get_shooting_point(step):
    details = step.change.canonical.details
    try:
        shooting_pt = details.shooting_snapshot
    except AttributeError:
        shooting_pt = NOT_EXTRACTED
    return shooting_pt

shooting_steps = StepFilter(
   lambda step: _get_shooting_point(step) is not NOT_EXTRACTED,
   name="shooting_steps"
)
"""
Filter to select steps 
"""

shooting_points = Extractor(_get_shooting_point, name="shooting_points")

def _get_modified_shooting_point(step):
    if not ShootingSteps.condition(step):
        return NOT_EXTRACTED
    details = step.change.canonical.details
    try:
        shooting_pt = details.modified_shooting_snapshot
    except AttributeError:
        # return the original shooting point for one-way
        shooting_pt = details.shooting_snapshot
    return shooting_pt

modified_shooting_points = Extractor(_get_modified_shooting_point,
                                     name="modified_shooting_points")

canonical_movers = Extractor(lambda step: step.change.canonical.mover,
                             name="canonical_movers")
