"""
Experimental analysis module.

The idea here is to simplify the vast majority of common analysis routines.
Interestingly, the process should also simplify a lot of calculation
preparation.

Goal: RETIS for a simple A->B transition (one direction) boils down to

>>> import openpathsampling as paths
>>> engine = ??? something that sets up the MD engine
>>> storage = ??? something that sets up storage
>>> globalstate0 = ??? something that sets up 
>>> orderparameter = paths.OP_Function("lambda", some_function)
>>>
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

class Transition(object):
    """
    Describes (in general) a transition between two states.
    """
    def __init__(self, stateA, stateB, storage=None):
        self.movers = {}
        pass

    def load_samples(self, storage=None):
        pass

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

        self._calcd_crossprob_params = None # check if this changed
        self._calcd_pathlen_params = None # check if this changed 

        # caches for the results of our calculations
        self._individual_crossing_probabilities = None
        self._total_crossing_probability = None
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

    @property
    def crossing_probability_parameters(self):
        """Dictionary of parameters for crossing probabilities"""
        pass

    # analysis results: note that these are cached, but can be overridden by
    # using `force=True` (or by changing the nature of histogram parameters)
    def individual_crossing_probabilities(self, force=False):
        """Return the crossing probability for each interface."""
        params_match = (self.crossing_probability_parameters == self._calcd_crossprob_params)
        if force == True or params_match == False:
            # redo the calculation
            pass

        return self._individual_crossing_probabilities

    def total_crossing_probability(self, method="wham", force=False):
        """Return the total crossing probability using `method`"""
        pass

    def pathlength_histograms(self, force=False):
        """Return histogram of the path length for each ensemble"""
        pass

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
            flux = self.minus_move_flux()

        pass

    def default_movers(self, engine):
        """Create reasonable default movers for a `PathSampling` calculation
        
        Extends `TISTransition.default_movers`.
        """
        pass

class RETISBirectionalSetup(object):
    def __init__(self, A_to_B, B_to_A):
        pass
