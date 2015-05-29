from histogram import Histogram, histograms_to_pandas_dataframe
from wham import WHAM
import numpy as np
from lookup_function import LookupFunction
import openpathsampling as paths
from openpathsampling.todict import OPSNamed

import logging
logger = logging.getLogger(__name__)


def pathlength(sample):
    return len(sample.trajectory)

def max_lambdas(sample, orderparameter):
    return max(orderparameter(sample.trajectory))

def sampleset_sample_generator(storage):
    for step in storage.steps:
        sset = step.active # take the sampleset after the move
        for sample in sset:
            yield sample

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

class Transition(OPSNamed):
    """
    Describes (in general) a transition between two states.
    """
    def __init__(self, stateA, stateB):
        super(Transition, self).__init__()
        self.movers = {}
        self.stateA = stateA
        self.stateB = stateB

    @property
    def all_movers(self):
        """
        All the path movers for this transition.
        """
        all_movers = []
        for movetype in self.movers.keys():
            all_movers += self.movers[movetype]
        return all_movers

    # TODO: @dwhs can this function be removed?
    # def _assign_sample(self, sample):
    #     if sample.ensemble in self.all_ensembles:
    #         try:
    #             self._samples_by_ensemble[sample.ensemble].append(sample)
    #         except KeyError:
    #             self._samples_by_ensemble[sample.ensemble] = [sample]
    #         try:
    #             self._samples_by_id[sample.replica].append(sample)
    #         except KeyError:
    #             self._samples_by_id[sample.replica] = [sample]

    # TODO: can this be removed? Only used in _assign_sample
    # @property
    # def all_ensembles(self):
    #     return self.ensembles


    def to_dict(self):
        return {
            'stateA' : self.stateA,
            'stateB' : self.stateB,
            'movers' : self.movers
        }

    @classmethod
    def from_dict(cls, dct):
        return Transition(
            stateA=dct['stateA'],
            stateB=dct['stateB']
        )

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
    orderparameter : CollectiveVariable
        order parameter to be used in the analysis (does not need to be the
        parameter which defines the interfaces, although it usually is)
    name : string
        name for the transition

    """
    
    def __init__(self, stateA, stateB, interfaces, orderparameter=None, name=None):
        super(TISTransition, self).__init__(stateA, stateB)
        # NOTE: making these into dictionaries like this will make it easy
        # to combine them in order to make a PathSampling PathSimulator object


        self.stateA = stateA
        self.stateB = stateB
        self.interfaces = interfaces
        self.name = name

        # If we reload from a storage file, we want to use the
        # ensembles/movers from the file, not the automatically generated
        # ones here

        # build ensembles if we don't already have them
        self.orderparameter = orderparameter
        if not hasattr(self, "ensembles"):
            self.build_ensembles(self.stateA, self.stateB, 
                                 self.interfaces, self.orderparameter)

        # build movers if we don't already have them
        if self.movers == {}:
            self.build_movers()

        self.default_orderparameter = self.orderparameter

        self.total_crossing_probability_method="wham" 
        self.histograms = {}
        # caches for the results of our calculation
        self._flux = None
        self._rate = None

        self.hist_args = {} # shortcut to ensemble_histogram_info[].hist_args
        self.ensemble_histogram_info = {
            'max_lambda' : Histogrammer(
                f=max_lambdas,
                f_args={'orderparameter' : self.orderparameter},
                hist_args={}
            ),
            'pathlength' : Histogrammer(
                f=pathlength,
                f_args={},
                hist_args={}
            )
        }

    # TODO: replace with copy.copy()
    def copy(self, with_results=True):
        copy = self.from_dict(self.to_dict())
        copy.default_orderparameter = self.default_orderparameter
        copy.total_crossing_probability_method = self.total_crossing_probability_method
        copy.hist_args = self.hist_args
        copy.ensemble_histogram_info = self.ensemble_histogram_info
        copy.histograms = self.histograms
        copy._flux = self._flux
        copy._rate = self._rate
        if hasattr(self, "tcp"):
            copy.tcp = self.tcp
        if hasattr(self, "ctp"):
            copy.ctp = self.ctp
        return copy


    def __str__(self):
        mystr = str(self.__class__.__name__) + ": " + str(self.name) + "\n"
        mystr += (str(self.stateA.name) + " -> " + str(self.stateA.name) 
                  + " or " + str(self.stateB.name) + "\n")
        for iface in self.interfaces:
            mystr += "Interface: " + str(iface.name) + "\n"
        return mystr


    def to_dict(self):
        ret_dict = {
            'stateA' : self.stateA,
            'stateB' : self.stateB,
            'orderparameter' : self.orderparameter,
            'interfaces' : self.interfaces,
            'name' : self.name,
            'movers' : self.movers,
            'ensembles' : self.ensembles
        }
        return ret_dict

    @classmethod
    def from_dict(cls, dct):
        mytrans = paths.TISTransition(
            stateA=dct['stateA'],
            stateB=dct['stateB'],
            interfaces=dct['interfaces'],
            orderparameter=dct['orderparameter'],
            name=dct['name']
        )
        mytrans.movers = dct['movers']
        mytrans.ensembles = dct['ensembles']
        return mytrans

    def build_ensembles(self, stateA, stateB, interfaces, orderparameter):
        self.ensembles = paths.EnsembleFactory.TISEnsembleSet(
            stateA, stateB, self.interfaces, orderparameter
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

        for hist in run_it:
            hist_info = self.ensemble_histogram_info[hist]
            if hist_info.hist_args == {} and self.hist_args[hist] != {}:
                hist_info.hist_args = self.hist_args[hist]

            if hist not in self.histograms.keys():
                self.histograms[hist] = {}
            self.histograms[hist][ensemble] = Histogram(**(hist_info.hist_args))

        in_ens_samples = (s for s in samples if s.ensemble is ensemble)
        hist_data = {}
        buflen = -1
        sample_buf = []
        for sample in in_ens_samples:
            for hist in run_it:
                hist_info = self.ensemble_histogram_info[hist]
                hist_data_sample = hist_info.f(sample, **hist_info.f_args)
                try:
                    hist_data[hist].append(hist_data_sample)
                except KeyError:
                    hist_data[hist] = [hist_data_sample]


        for hist in run_it:
            self.histograms[hist][ensemble].histogram(hist_data[hist], weights)
            self.histograms[hist][ensemble].name = (hist + " " + self.name
                                                    + " " + ensemble.name)


    def all_statistics(self, storage, weights=None, force=False):
        """
        Run all statistics for all ensembles.
        """
        # TODO: speed this up by just running over all samples once and
        # dealing them out to the appropriate histograms
        for ens in self.ensembles:
            samples = sampleset_sample_generator(storage)
            self.ensemble_statistics(ens, samples, weights, force)

    def pathlength_histogram(self, ensemble):
        """
        Return the pathlength histogram for the given ensemble.
        """
        # check existence and correctness of self.histograms[pl][ens]
        if "pathlength" not in self.histograms:
            self.histograms['pathlength'] = {}
        hist = self.histograms['pathlength'][ensemble]
        return hist.normalized()

    def crossing_probability(self, ensemble, nblocks=1):
        """
        Return the crossing probability for the given ensemble.
        """
        # check existence and correctness of self.histograms[cp][ens]
        hist = self.histograms['crossing_probability'][ensemble]
        return hist.reverse_cumulative()

    def total_crossing_probability(self, storage=None, method="wham", force=False):
        """Return the total crossing probability using `method`
        
        Parameters
        ----------
        storage : storage
            cycles to be analyzed
        method : "wham" (later: or "mbar" or "tram")
            approach to use to combine the histograms
        force : bool (False)
            if true, cached results are overwritten
        """

        tcp = dict()

        if method == "wham":
            run_ensembles = False
            for ens in self.ensembles:
                try:
                    hist = self.histograms['max_lambda'][ens]
                except KeyError:
                    run_ensembles = True
            if run_ensembles or force:
                if storage is None:
                    raise RuntimeError("Unable to build histograms without storage source")
                self.all_statistics(storage, force=True)
                         
            df = histograms_to_pandas_dataframe(
                self.histograms['max_lambda'].values(),
                fcn="reverse_cumulative"
            ).sort(axis=1)
            wham = WHAM()
            wham.load_from_dataframe(df)
            wham.clean_leading_ones()
            tcp = wham.wham_bam_histogram()
        elif method == "mbar":
            pass
        else:
            raise ValueError("Only supported methods are 'wham' and 'mbar'. " + \
                             "Whereas 'mbar' is not yet implemented!")

        self.tcp = LookupFunction(tcp.keys(), tcp.values())
        return self.tcp

    def conditional_transition_probability(self, storage, ensemble, force=False):
        """
        This transition's conditional transition probability for a given
        ensemble.

        The conditional transition probability for an ensemble is the
        probability that a path in that ensemble makes the transition from
        state A to state B.

        Parameters
        ----------
        storage : storage
            cycles to analyze
        ensemble : Ensemble
            which ensemble to calculate the CTP for
        force : bool (False)
            if true, cached results are overwritten
        """
        samples = sampleset_sample_generator(storage)
        n_acc = 0
        n_try = 0
        for samp in samples:
            if samp.ensemble is ensemble:
                if self.stateB(samp.trajectory[-1]):
                    n_acc += 1
                n_try += 1
        ctp = float(n_acc)/n_try
        logger.info("CTP: " + str(n_acc) + "/" + str(n_try) + "=" + str(ctp) + "\n")
        try:
            self.ctp[ensemble] = ctp
        except AttributeError:
            self.ctp = { ensemble : ctp }

        return ctp

    def rate(self, storage, flux=None, outer_ensemble=None,
             outer_lambda=None, error=None, force=False):
        """Calculate the rate for this transition.

        For TIS transitions, this requires the result of an external
        calculation of the flux.

        Parameters
        ==========
        storage : openpathsampling.storage.Storage
        flux : float
        outer_ensemble : openpathsampling.TISEnsemble
        error : list(3) or None
        """
        # get the flux
        if flux is not None:
            self._flux = flux

        if self._flux is None:
            raise ValueError("No flux available to TISTransition. Cannot calculate rate")
        
        flux = self._flux

        # get the total crossing probability
        if not force and hasattr(self, 'tcp'):
            tcp = self.tcp
        else:
            tcp = self.total_crossing_probability(storage=storage, force=force)

        # get the conditional transition probability
        if outer_ensemble is None:
            outer_ensemble = self.ensembles[-1]
        outer_cross_prob = self.histograms['max_lambda'][outer_ensemble]
        if outer_lambda is None:
            lambda_bin = -1
            outer_cp_vals = outer_cross_prob.reverse_cumulative().values()
            while (outer_cp_vals[lambda_bin+1] == 1.0):
                lambda_bin += 1
            outer_lambda = outer_cross_prob.bins[lambda_bin]

        ctp = self.conditional_transition_probability(storage,
                                                      outer_ensemble,
                                                      force=force)
        outer_tcp = tcp(outer_lambda)
        #print flux, outer_tcp, ctp
        self._rate = flux*outer_tcp*ctp
        return self._rate

    def default_movers(self, engine):
        """Create reasonable default movers for a `PathSampling` pathsimulator"""
        shoot_sel = paths.RandomChoiceMover(
            movers=self.movers['shooting']
        )
        shoot_sel.name = "ShootingChooser"
        pathrev_sel = paths.RandomChoiceMover(
            movers=self.movers['pathreversal']
        )
        pathrev_sel.name = "ReversalChooser"
        root_mover = paths.RandomChoiceMover(
            movers=[shoot_sel, pathrev_sel], 
            weights=[1.0, 0.5]
        )
        root_mover.name = "RootMover"
        return root_mover

class RETISTransition(TISTransition):
    """Transition class for RETIS."""
    def __init__(self, stateA, stateB, interfaces, orderparameter=None, name=None):
        super(RETISTransition, self).__init__(stateA, stateB, interfaces, orderparameter, name)

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
            self.movers['minus'] = [paths.MinusMover(self.minus_ensemble, self.ensembles[0])]


    def to_dict(self):
        ret_dict = {
            'stateA' : self.stateA,
            'stateB' : self.stateB,
            'orderparameter' : self.orderparameter,
            'interfaces' : self.interfaces,
            'name' : self.name,
            'movers' : self.movers,
            'ensembles' : self.ensembles,
            'minus_ensemble' : self.minus_ensemble
        }
        return ret_dict

    @classmethod
    def from_dict(cls, dct):
        mytrans = RETISTransition(
            stateA=dct['stateA'],
            stateB=dct['stateB'],
            interfaces=dct['interfaces'],
            orderparameter=dct['orderparameter'],
            name=dct['name']
        )
        mytrans.minus_ensemble = dct['minus_ensemble']
        mytrans.movers = dct['movers']
        mytrans.ensembles = dct['ensembles']
        return mytrans

    @property
    def all_ensembles(self):
        return self.ensembles + [self.minus_ensemble]


    def minus_move_flux(self, storage, force=False):
        """
        Calculate the flux based on the minus ensemble trajectories.
        """
        if not force and self._flux != None:
            return self._flux

        self.minus_count_sides = { "in" : [], "out" : [] }
        minus_moves = (d.change for d in storage.steps
                       if self.movers['minus'][0] in
                       d.change and d.change.accepted)
        for move in minus_moves:
            minus_samp = [s for s in move.results
                          if s.ensemble is self.minus_ensemble][0]
            minus_trajectory = minus_samp.trajectory
            minus_summ = minus_sides_summary(minus_trajectory,
                                             self.minus_ensemble)
            for key in self.minus_count_sides.keys():
                self.minus_count_sides[key].extend(minus_summ[key])
       
        for key in self.minus_count_sides.keys():
            if len(self.minus_count_sides[key]) == 0:
                logger.warn("No instances of "+str(key)+" for minus move.")

        t_in_avg = np.array(self.minus_count_sides['in']).mean()
        t_out_avg = np.array(self.minus_count_sides['out']).mean()
        self._flux = 1.0 / (t_in_avg + t_out_avg)
        return self._flux


    def rate(self, storage, flux=None, outer_ensemble=None,
             outer_lambda=None, error=None, force=False):
        if flux is None:
            flux = self.minus_move_flux(storage)

        return super(RETISTransition, self).rate(
            storage=storage, 
            flux=flux, 
            outer_ensemble=outer_ensemble,
            outer_lambda=outer_lambda,
            error=error,
            force=force
        )

    def populate_minus_ensemble(self, partial_traj, minus_replica_id, engine):
        """
        Generate a sample for the minus ensemble by extending `partial_traj`

        Parameters
        ----------
        partial_traj : Trajectory
            trajectory to extend
        minus_replica_id : integer or string
            replica ID for this sample
        engine : DynamicsEngine
            engine to use for MD extension
        """
        last_frame = partial_traj[-1]
        if not self.minus_ensemble._segment_ensemble(partial_traj):
            raise RuntimeError(
                "Invalid input trajectory for minus extension. (Not A-to-A?)"
            )
        extension = engine.generate(last_frame,
                                    [self.minus_ensemble.can_append])
        first_minus = paths.Trajectory(partial_traj + extension[1:])
        minus_samp = paths.Sample(
            replica=minus_replica_id,
            trajectory=first_minus,
            ensemble=self.minus_ensemble
        )
        return minus_samp

    def default_movers(self, engine):
        """Create reasonable default movers for a `PathSampling` pathsimulator
        
        Extends `TISTransition.default_movers`.
        """
        repex_sel = paths.RandomChoiceMover(
            movers=self.movers['repex']
        )
        repex_sel.name = "ReplicaExchange"
        tis_root_mover = super(RETISTransition, self).default_movers(engine)
        minus = self.movers['minus']
        movers = tis_root_mover.movers + [repex_sel] + minus
        weights = tis_root_mover.weights + [0.5, 0.2 / len(self.ensembles)]
        root_mover = paths.RandomChoiceMover(
            movers=movers,
            weights=weights
        )
        root_mover.name = "RootMover"
        return root_mover


def minus_sides_summary(trajectory, minus_ensemble):
    # note: while this could be refactored so vol_dict is external, I don't
    # think this hurts speed very much, and it it really useful for testing
    minus_state = minus_ensemble.state_vol
    minus_innermost = minus_ensemble.innermost_vol
    minus_interstitial = minus_innermost & ~minus_state
    vol_dict = { 
        "A" : minus_state,
        "X" : ~minus_innermost,
        "I" : minus_interstitial
    }
    summary = trajectory.summarize_by_volumes(vol_dict)
    # this is a per-trajectory loop
    count_sides = {"in" : [], "out" : []}
    side=None
    local_count = 0
    # strip off the beginning and ending in A
    for (label, count) in summary[1:-1]:
        if label == "X" and side != "out":
            if side == "in":
                count_sides["in"].append(local_count)
            side = "out"
            local_count = 0
        elif label == "A" and side != "in":
            if side == "out":
                count_sides["out"].append(local_count)
            side = "in"
            local_count = 0
        local_count += count
    return count_sides
