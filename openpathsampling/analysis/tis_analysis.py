"""
Experimental analysis module.

The idea here is to simplify the vast majority of common analysis routines.
Interestingly, the process should also simplify a lot of calculation
preparation.

Goal: RETIS for a simple A->B transition (one direction) boils down to

>>> # define orderparameter; stateA; stateB; engine; storage
>>> interfaces = VolumeSet(orderparameter, min=-infinity, max=[0.0, 0.1, 0.2])
>>> transition = RETISTransition(stateA, stateB, interfaces, storage)
>>> retis_calc = transition.calculation(engine, globalstate0)
>>> retis_calc.run(nsteps=10000)
>>> tcp = retis_calc.total_crossing_probability()
>>> flow = retis_calc.replica_flow()
>>> rate = retis_calc.rate()
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
    def __init__(self, stateA, stateB, interfaces, storage=None):
        super(TISTransition, self).__init__(stateA, stateB, storage)
        # NOTE: making these into dictionaries like this will make it easy
        # to combine them in order to make a PathSampling calculation object
        self.movers['shooting'] = []
        self.movers['pathreversal'] = []

        self._calcd_crossprob_params = None # check if this changed
        self._calcd_pathlen_params = None # check if this changed 

        # these get set when we run the calculation
        self._individual_crossing_probabilities = None
        self._total_crossing_probability = None
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

    # properties which define the output format of analyzed data
    @property
    def individual_crossing_probabilities(self, force=False):
        """Return the crossing probability for each interface."""
        params_match = (self.crossing_probability_parameters == self._calcd_crossprob_params)
        if force == True or params_match == False:
            # redo the calculation
            pass

        return self._individual_crossing_probabilities

    @property
    def total_crossing_probability(self, method="wham", force=False):
        """Return the total crossing probability using `method`"""
        pass

    @property
    def pathlength_histograms(self), force=False):
        """Return histogram of the path length for each ensemble"""
        pass

    @property
    def rate(self, flux, flux_error=None, force=False):
        pass



class RETISTransition(TISTransition):
    """
    Transition class for RETIS
    """
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
        pass

