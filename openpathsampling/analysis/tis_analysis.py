from histogram import Histogram
import openpathsampling as paths
from openpathsampling.util.todict import restores_as_full_object

"""
Experimental analysis module.

The idea here is to simplify the vast majority of common analysis routines.
Interestingly, the process should also simplify a lot of pathsimulator
preparation.

Goal: RETIS for a simple A->B transition (one direction) boils down to

>>> # things that would be hypothetically already set up
>>> import openpathsampling as paths
>>> engine = ??? something that sets up the MD engine
>>> storage = ??? something that sets up storage
>>> globalstate0 = ??? something that sets up initial trajectories
>>> collectivevariable = paths.CV_Function("lambda", some_function)
>>>
>>> # from here, this is real code
>>> stateA = paths.LambdaVolume(collectivevariable, min=-infinity, max=0.0)
>>> stateB = paths.LambdaVolume(collectivevariable, min=1.0, max=infinity)
>>> interfaces = paths.VolumeSet(collectivevariable, min=-infinity, max=[0.0, 0.1, 0.2])
>>> transitionAB = paths.RETISTransition(stateA, stateB, collectivevariable, interfaces, storage)
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
In the order listed above, the time for the rate pathsimulator is almost
entirely in determining the flux from the information in the minus mover.
"""


def pathlength(sample):
    return len(sample.trajectory)

def max_lambdas(sample, collectivevariable):
    return max([collectivevariable(frame) for frame in sample.trajectory])


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
    def __init__(self, stateA, stateB):
        self.movers = {}
        self.stateA = stateA
        self.stateB = stateB

        self._mover_acceptance = {}
        pass

    def calculate_mover_acceptance(self, samples):
        for sample in samples:
            pass
        pass


class TPSTransition(Transition):
    """
    Transition using TPS ensembles
    """
    def __init__(self, stateA, stateB):
        super(TPSTransition, self).__init__(stateA, stateB)
        self.movers['shooting'] = []
        self.movers['shifting'] = []
        self.movers['pathreversal'] = []
        #self.ensembles = [paths.TPSEnsemble(stateA, stateB)]


@restores_as_full_object
class TISTransition(Transition):
    """
    Transition using TIS ensembles.

    The additional information from the TIS ensembles allows us to set up
    all the analysis (assuming we built these are proper TIS ensembles,
    which we DO in the intitialization!)

    Parameters
    ----------
    stateA : Volume
        Volume for the state from which the transition begins
    stateB : Volume
        Volume for the state in which the transition ends
    interfaces : list of Volume
        Volumes for the interfaces
    collectivevariable : CollectiveVariable
        order parameter to be used in the analysis (does not need to be the
        parameter which defines the interfaces, although it usually is)
    name : string
        name for the transition

    """
    
    def __init__(self, stateA, stateB, interfaces, collectivevariable=None, name=None):
        super(TISTransition, self).__init__(stateA, stateB)
        # NOTE: making these into dictionaries like this will make it easy
        # to combine them in order to make a PathSampling pathsimulator object

        self.stateA = stateA
        self.stateB = stateB
        self.interfaces = interfaces
        self.name = name

        # If we reload from a storage file, we want to use the
        # ensembles/movers from the file, not the automatically generated
        # ones here

        # build ensembles if we don't already have them
        if not hasattr(self, "ensembles"):
            self.build_ensembles(self.stateA, self.stateB, self.interfaces)

        # build movers if we don't already have them
        if self.movers == {}:
            self.build_movers()

        self.collectivevariable = collectivevariable
        self.default_collectivevariable = self.collectivevariable

        self.total_crossing_probability_method="wham" 
        self.histograms = {}
        self._ensemble_histograms = {}
        # caches for the results of our calculations
        self._flux = None
        self._rate = None

        # TODO: eventually I'll generalize this to include the function to
        # be called, possibly some parameters ... can't this go to a 
        self.ensemble_histogram_info = {
            'max_lambda' : Histogrammer(
                f=max_lambdas,
                f_args={'collectivevariable' : self.collectivevariable},
                hist_args={}
            ),
            'pathlength' : Histogrammer(
                f=pathlength,
                f_args={},
                hist_args={}
            )
        }


    def to_dict(self):
        ret_dict = {
            'stateA' : self.stateA,
            'stateB' : self.stateB,
            'collectivevariable' : self.collectivevariable,
            'interfaces' : self.interfaces,
            'name' : self.name,
            'movers' : self.movers,
            'ensembles' : self.ensembles
        }
        return ret_dict

    @staticmethod
    def from_dict(adict):
        mytrans = TISTransition(
            stateA=adict['stateA'],
            stateB=adict['stateB'],
            interfaces=adict['interfaces'],
            collectivevariable=adict['collectivevariable'],
            name=adict['name']
        )
        mytrans.movers = adict['movers']
        mytrans.ensembles = adict['ensembles']
        return mytrans

    def build_ensembles(self, stateA, stateB, interfaces):
        self.ensembles = paths.EnsembleFactory.TISEnsembleSet(
            stateA, stateB, self.interfaces
        )
        for ensemble in self.ensembles:
            ensemble.name = "I'face "+str(self.ensembles.index(ensemble))

    def build_movers(self):
        self.movers['shooting'] = paths.PathMoverFactory.OneWayShootingSet(
            paths.UniformSelector(), self.ensembles
        )
        self.movers['pathreversal'] = paths.PathReversalSet(self.ensembles)

    # parameters for different types of output
    def ensemble_statistics(self, ensemble, samples, weights=None, force=False):
        """Calculate stats for a given ensemble: path length, crossing prob

        In general we do all of these at once because the extra cost of
        running through the samples twice is worse than doing the extra
        calculations.

        Parameters
        ----------
        ensemble: Ensemble
        samples : iterator over samples

        """
        # figure out which histograms need to updated for this ensemble
        run_it = []
        if not force:
            # TODO figure out which need to be rerun
            pass
        else:
            run_it = self.ensemble_histogram_info.keys()

        buflen = 10
        for hist in run_it:
            in_ens_samples = (s for s in samples if s.ensemble == ensemble)
            hist_info = self.ensemble_histogram_info[hist]
            if hist not in self.histograms.keys():
                self.histograms[hist] = {}
            self.histograms[hist][ensemble] = Histogram(**(hist_info.hist_args))
            hist_data = []
            for sample in in_ens_samples:
                hist_data.append(hist_info.f(sample, **hist_info.f_args))
            self.histograms[hist][ensemble].histogram(hist_data, weights)
            self.histograms[hist][ensemble].name = (hist + " " + self.name
                                                    + " " + ensemble.name)


    def all_statistics(self, samples, weights=None, force=False):
        # TODO: speed this up by just running over all samples once and
        # dealing them out to the appropriate histograms
        for ens in self.ensembles:
            self.ensemble_statistics(ens, samples, weights, force)

    def pathlength_histogram(self, ensemble):
        # check existence and correctness of self.histograms[pl][ens]
        if "pathlength" not in self.histograms:
            self.histograms['pathlength'] = {}
        hist = self.histograms['pathlength'][ensemble]
        return hist.normalized()

    def crossing_probability(self, ensemble, nblocks=1):
        # check existence and correctness of self.histograms[cp][ens]
        hist = self.histograms['crossing_probability'][ensemble]
        return hist.reverse_cumulative()

    def total_crossing_probability(self, method="wham", force=False, nblocks=1):
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

        self.tcp = tcp
        return tcp

    def rate(self, flux=None, flux_error=None, force=False):
        """Calculate the rate for this transition.

        For TIS transitions, this requires the result of an external
        pathsimulator of the flux.
        """
        if flux is not None:
            self._flux = flux

        if self._flux is None:
            raise ValueError("No flux available to TISTransition. Cannot calculate rate")
        
        tcp = self.total_crossing_probability(force=force)
        pass

    def default_movers(self, engine):
        """Create reasonable default movers for a `PathSampling` pathsimulator"""
        shoot_sel = paths.RandomChoiceMover(
            movers=self.movers['shooting'],
            name="ShootingChooser"
        )
        pathrev_sel = paths.RandomChoiceMover(
            movers=self.movers['pathreversal'],
            name="ReversalChooser"
        )
        root_mover = paths.RandomChoiceMover(
            movers=[shoot_sel, pathrev_sel], 
            weights=[1.0, 0.5],
            name="RootMover"
        )
        return root_mover


class RETISTransition(TISTransition):
    """Transition class for RETIS."""
    def __init__(self, stateA, stateB, interfaces, collectivevariable=None, name=None):
        super(RETISTransition, self).__init__(stateA, stateB, interfaces,
                                              collectivevariable, name)

        self.minus_ensemble = paths.MinusInterfaceEnsemble(
            state_vol=stateA, 
            innermost_vol=interfaces[0]
        )

        try:
            self.movers['repex']
        except KeyError:
            self.movers['repex'] = paths.NeighborEnsembleReplicaExchange(self.ensembles)
        try:
            self.movers['minus']
        except KeyError:
            self.movers['minus'] = paths.MinusMover(self.minus_ensemble, self.ensembles[0])


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
        # TODO: move this to network
        pass

    @property
    def rate(self, flux=None, flux_error=None, force=False):
        tcp = self.total_crossing_probability()
        if flux is None:
            (flux, flux_error) = self.minus_move_flux()

        pass

    def default_movers(self, engine):
        """Create reasonable default movers for a `PathSampling` pathsimulator
        
        Extends `TISTransition.default_movers`.
        """
        repex_sel = paths.RandomChoiceMover(
            movers=self.movers['repex'],
            name="ReplicaExchange"
        )
        tis_root_mover = super(RETISTransition, self).default_movers(engine)
        movers = tis_root_mover.movers + [repex_sel, self.movers['minus']]
        weights = tis_root_mover.weights + [0.5, 0.2 / len(self.ensembles)]
        root_mover = paths.RandomChoiceMover(
            movers=movers,
            weights=weights,
            name="RootMover"
        )
        return root_mover

