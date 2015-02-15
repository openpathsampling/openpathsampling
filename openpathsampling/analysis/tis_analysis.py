from histogram import Histogram
"""
Experimental analysis module.

The idea here is to simplify the vast majority of common analysis routines.
Interestingly, the process should also simplify a lot of calculation
preparation.

Goal: RETIS for a simple A->B transition (one direction) boils down to

>>> # things that would be hypothetically already set up
>>> import openpathsampling as paths
>>> engine = ??? something that sets up the MD engine
>>> storage = ??? something that sets up storage
>>> globalstate0 = ??? something that sets up initial trajectories
>>> orderparameter = paths.OP_Function("lambda", some_function)
>>>
>>> # from here, this is real code
>>> stateA = paths.LambdaVolume(orderparameter, min=-infinity, max=0.0)
>>> stateB = paths.LambdaVolume(orderparameter, min=1.0, max=infinity)
>>> interfaces = paths.VolumeSet(orderparameter, min=-infinity, max=[0.0, 0.1, 0.2])
>>> transitionAB = paths.RETISTransition(stateA, stateB, orderparameter, interfaces, storage)
>>> retis_calc = PathSampling(
>>>     storage=storage,
>>>     engine=engine,
>>>     root_mover=transitionAB.default_movers(engine),
>>>     globalstate=globalstate0
>>> )
>>> retis_calc.run(nsteps=10000)
>>> tcp = transitionAB.total_crossing_probability()
>>> flow = transitionAB.replica_flow()
>>> rate = transitionAB.rate()

Note that once the total crossing probability has been calculated once, it
does not need to be recalculated as part of the rate. (Or, if it were
calculated as part of the rate, it would be already available on its own.)
In the order listed above, the time for the rate calculation is almost
entirely in determining the flux from the information in the minus mover.
"""

def pathlengths(samples):
    for sample in samples:
        yield len(sample.trajectory)

def max_lambdas(sample, orderparameter):
    for sample in samples:
        yield max([orderparameter(frame) for frame in sample.trajectory])


class Histogrammer(object):
    """
    Basically a dictionary to track what each histogram should be making.
    """
    def __init__(self, f, f_args=None, hist_args=None):
        self.f = f
        self.f_args = f_args
        self._hist_args = hist_args
        self.empty_hist = Histogram(**self._hist_args)

    @property
    def hist_args(self):
        return self._hist_args

    @hist_args.setter
    def hist_args(self, val):
        self._hist_args = val
        self.empty_hist = Histogram(**self._hist_args)


class Transition(object):
    """
    Describes (in general) a transition between two states.
    """
    def __init__(self, stateA, stateB, storage=None):
        self.movers = {}
        self.stateA = stateA
        self.stateB = stateB
        self.storage = storage
        pass


class TPSTransition(object):
    """
    Transition using TPS ensembles
    """
    def __init__(self, stateA, stateB, storage=None):
        super(TPSTransition, self).__init__(stateA, stateB, storage)
        self.movers['shooting'] = []
        self.movers['shifting'] = []
        self.movers['pathreversal'] = []
        #self.ensembles = [paths.TPSEnsemble(stateA, stateB)]


class TISTransition(Transition):
    """
    Transition using TIS ensembles.

    The additional information from the TIS ensembles allows us to set up
    all the analysis (assuming we built these are proper TIS ensembles,
    which we DO in the intitialization!)
    """
    def __init__(self, stateA, stateB, orderparameter, interfaces, storage=None):
        super(TISTransition, self).__init__(stateA, stateB, storage)
        # NOTE: making these into dictionaries like this will make it easy
        # to combine them in order to make a PathSampling calculation object
        self.movers['shooting'] = []
        self.movers['pathreversal'] = []

        self.total_crossing_probability_method="wham" 

        # TODO: eventually I'll generalize this to include the function to
        # be called, possibly some parameters ... can't this go to a 
        self.ensemble_histogram_info = {
            'pathlength_histogram' : Histogrammer(
                f=sample_pathlength,
                f_args=None,
                hist_args=None
            )
        }
        self.histograms = {}

        # caches for the results of our calculations
        self._flux = None
        self._rate = None
        pass

    # path movers
    @property
    def shooting_movers(self):
        return self.movers['shooting']

    @property
    def pathreversal_movers(self):
        return self.movers['pathreversal']

    # parameters for different types of output
    def ensemble_statistics(self, ensemble, data, weights, force=False):
        """Calculate stats for a given ensemble: path length, crossing prob

        In general we do all of these at once because the extra cost of
        running through the samples twice is worse than doing the extra
        calculations.
        """
        # figure out which histograms need to updated for this ensemble
        run_it = []
        if not force:
            # figure out which need to be rerun
            pass
        else:
            run_it = self._ensemble_histograms.keys()

        for hist in run_it:
            hist_info = self.ensemble_histogram_info[hist]
            if hist not in self.histograms:
                self.histograms[hist] = {}
            self.histograms[hist][ensemble] = Histogram(**(hist_info.hist_args))
            self.histograms[hist][ensemble].histogram(
                hist_info.f(data, **hist_info.f_args), weights
            )

        pass

    def pathlength_histogram(self, ensemble):
        # check existence and correctness of self.histograms[pl][ens]
        hist = self.histograms['pathlength'][ensemble]
        return hist.normalized()

    def crossing_probability(self, ensemble):
        # check existence and correctness of self.histograms[cp][ens]
        hist = self.histograms['crossing_probability'][ensemble]
        return hist.reverse_cumulative()

    def total_crossing_probability(self, method="wham", force=False):
        """Return the total crossing probability using `method`"""
        if method == "wham":
            cp = {}
            for ens in self.ensembles:
                cp[ens] = self.crossing_probability(ens)
            wham = WHAM()
            wham.initial_histograms = cp
            wham.clean_leading_ones()
            tcp = wham.wham_bam_histogram()
        elif method == "mbar":
            pass

        return tcp

    def rate(self, flux=None, flux_error=None, force=False):
        """Calculate the rate for this transition.

        For TIS transitions, this requires the result of an external
        calculation of the flux. 
        """
        if flux is not None:
            self._flux = flux

        if self._flux is None:
            raise ValueError("No flux available to TISTransition. Cannot calculate rate")
        
        tcp = self.total_crossing_probability(force=force)
        pass

    def default_movers(self, engine):
        """Create reasonable default movers for a `PathSampling` calculation"""
        pass



class RETISTransition(TISTransition):
    """Transition class for RETIS."""
    def __init__(self, stateA, stateB, interfaces, storage=None):
        super(RETISTransition, self).__init__(stateA, stateB, interfaces, storage)
        self.movers['repex'] = []
        self.movers['minus'] = []

    @property
    def repex_movers(self):
        pass

    @property
    def minus_movers(self):
        pass

    @property
    def replica_flow(self):
        # this should check to build the replica exchange network. If the
        # number of neighbors at any station is more than 2, we can't do
        # "normal" replica flow -- instead produce a network graph. Or,
        # actually, ALWAYS produce a network graph (although this will be a
        # feature to implement later)
        pass

    @property
    def minus_move_flux(self):
        pass

    @property
    def multiple_set_minus_switching(self):
        pass

    @property
    def rate(self, flux=None, flux_error=None, force=False):
        tcp = self.total_crossing_probability()
        if flux is None:
            (flux, flux_error) = self.minus_move_flux()

        pass

    def default_movers(self, engine):
        """Create reasonable default movers for a `PathSampling` calculation
        
        Extends `TISTransition.default_movers`.
        """
        pass

class RETISBirectionalSetup(object):
    def __init__(self, A_to_B, B_to_A):
        pass
