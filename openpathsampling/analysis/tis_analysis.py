from histogram import Histogram, histograms_to_pandas_dataframe
from wham import WHAM
import numpy as np
from lookup_function import LookupFunction
import openpathsampling as paths
from openpathsampling.todict import OPSNamed
import sys

import inspect

import time 

"""
Experimental analysis module.

The idea here is to simplify the vast majority of common analysis routines.
Interestingly, the process should also simplify a lot of simulation
preparation.

Goal: RETIS for a simple A->B transition (one direction) boils down to

>>> # things that would be hypothetically already set up
>>> import openpathsampling as paths
>>> engine = ??? something that sets up the MD engine
>>> storage = ??? something that sets up storage
>>> globalstate0 = ??? something that sets up initial trajectories
>>> orderparameter = paths.CV_Function("lambda", some_function)
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

import logging
logger = logging.getLogger(__name__)


def pathlength(sample):
    return len(sample.trajectory)

def max_lambdas(sample, orderparameter):
    return max([orderparameter(frame) for frame in sample.trajectory])

def sampleset_sample_generator(storage):
    for sset in storage.samplesets:
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

        self._mover_acceptance = {}

    @property
    def all_movers(self):
        all_movers = []
        for movetype in self.movers.keys():
            all_movers += self.movers[movetype]
        return all_movers

    def _assign_sample(self, sample):
        if sample.ensemble in self.all_ensembles:
            try:
                self._samples_by_ensemble[sample.ensemble].append(sample)
            except KeyError:
                self._samples_by_ensemble[sample.ensemble] = [sample]
            try: 
                self._samples_by_id[sample.replica].append(sample)
            except KeyError:
                self._samples_by_id[sample.replica] = [sample]

    def initialize_with_storage(self, storage, force=False):
        if self._samples_by_ensemble == { }:
            force=True
        if force:
            for sample in sampleset_sample_generator(storage):
                self._assign_sample(sample)

    def _move_summary_line(self, move_name, n_accepted, n_trials,
                           n_total_trials, indentation):
        line = ("* "*indentation + str(move_name) +
                " ran " + str(float(n_trials)/n_total_trials*100) + 
                "% of the cycles with acceptance " + str(n_accepted) + "/" + 
                str(n_trials) + " (" + str(float(n_accepted) / n_trials) + 
                ") \n")
        return line

    def move_acceptance(self, storage):
        for step in storage.steps:
            delta = step.change
            for m in delta:
                acc = 1 if m.accepted else 0
                key = (m.mover, str(delta.key(m)))
                try:
                    self._mover_acceptance[key][0] += acc
                    self._mover_acceptance[key][1] += 1
                except KeyError:
                    self._mover_acceptance[key] = [acc, 1]

    def move_summary(self, storage, movers=None, output=sys.stdout, depth=0):
        """
        Provides a summary of the movers in `storage` based on this transition.

        The summary includes the number of moves attempted and the
        acceptance rate. In some cases, extra lines are printed for each of
        the submoves.

        Parameters
        ----------
        storage : Storage
            The storage object
        movers : None or string or list of PathMover
            If None, provides a short summary of the keys in self.mover. If
            a string, provides a short summary using that string as a key in
            the `movers` dict. If a mover or list of movers, provides
            summary for each of those movers.
        output : file
            file to direct output
        depth : integer or None
            depth of submovers to show: if integer, shows that many
            submovers for each move; if None, shows all submovers
        """
        my_movers = { }
        if movers is None:
            movers = self.movers.keys()
        if type(movers) is str:
            movers = self.movers[movers]
        for key in movers:
            try:
                my_movers[key] = self.movers[key]
            except KeyError:
                my_movers[key] = [key]


        stats = { } 
        for groupname in my_movers.keys():
            stats[groupname] = [0, 0]

        if self._mover_acceptance == { }:
            self.move_acceptance(storage)

        tot_trials = len(storage.steps)
        for groupname in my_movers.keys():
            group = my_movers[groupname]
            for mover in group:
                key_iter = (k for k in self._mover_acceptance.keys()
                            if k[0] == mover)

                for k in key_iter:
                    stats[groupname][0] += self._mover_acceptance[k][0]
                    stats[groupname][1] += self._mover_acceptance[k][1]

        for groupname in my_movers.keys():
            line = self._move_summary_line(
                move_name=groupname, 
                n_accepted=stats[groupname][0],
                n_trials=stats[groupname][1], 
                n_total_trials=tot_trials,
                indentation=0
            )
            output.write(line)

    @property
    def all_movers(self):
        all_movers = []
        for movetype in self.movers.keys():
            all_movers += self.movers[movetype]
        return all_movers

    @property
    def all_ensembles(self):
        return self.ensembles


    def calculate_mover_acceptance(self, storage, movers=None):
        # regularlize format of movers argument:
        if movers is None:
            movers = self.all_movers
        try:
            nmovers = len(movers)
        except TypeError:
            nmovers = 1
            movers = [movers]

        for movechange in storage.movepath:
            pass
        pass


    def to_dict(self):
        return {
            'stateA' : self.stateA,
            'stateB' : self.stateB,
            'movers' : self.movers
        }

    @staticmethod
    def from_dict(dct):
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
        self._ensemble_histograms = {}
        # caches for the results of our calculation
        self._flux = None
        self._rate = None

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

    @staticmethod
    def from_dict(dct):
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
            if hist not in self.histograms.keys():
                self.histograms[hist] = {}
            self.histograms[hist][ensemble] = Histogram(**(hist_info.hist_args))

        in_ens_samples = (s for s in samples if s.ensemble == ensemble)
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
        # TODO: speed this up by just running over all samples once and
        # dealing them out to the appropriate histograms
        for ens in self.ensembles:
            samples = sampleset_sample_generator(storage)
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

    def total_crossing_probability(self, method="wham", storage=None, force=False, nblocks=1):
        """Return the total crossing probability using `method`"""
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

        self.tcp = LookupFunction(tcp.keys(), tcp.values())
        return self.tcp

    def conditional_transition_probability(self, storage, ensemble, force=False):
        """
        This transition's conditional transition probability for a given
        ensemble.

        The conditional transition probability for an ensemble is the
        probability that a path in that ensemble makes the transition from
        state A to state B.
        """
        samples = sampleset_sample_generator(storage)
        n_acc = 0
        n_try = 0
        for samp in samples:
            if samp.ensemble == ensemble:
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
            while (outer_cross_prob.reverse_cumulative()[lambda_bin+1] == 1.0):
                lambda_bin += 1
            outer_lambda = outer_cross_prob.bins[lambda_bin]

        ctp = self.conditional_transition_probability(storage,
                                                      outer_ensemble,
                                                      force=force)
        outer_tcp = tcp(outer_lambda)
        #print flux, outer_tcp, ctp
        return flux*outer_tcp*ctp

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

    @staticmethod
    def from_dict(dct):
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


    @property
    def replica_flow(self, bottom=None, top=None):
        if bottom is None:
            bottom = self.minus_ensemble
        if top is None:
            top = self.ensembles[-1]

        updown = {}
        round_trips = {}
        for ens in [self.minus_ensemble] + self.ensembles:
            updown[ens] = 0
            nvisits[ens] = 0
            up_trips[ens] = []
            down_trips[ens] = []
            round_trips[ens] = []

        # loop over moves which include a repex

        # run stats to get replica flow and stats on up_trips, down_trips,
        # round_trips
        pass

    def minus_move_flux(self, storage, force=False):
        if not force and self._flux != None:
            return self._flux

        self.minus_count_sides = { "in" : [], "out" : [] }
        minus_moves = (d.change for d in storage.steps
                       if self.movers['minus'][0] in
                       d.change and d.change.accepted)
        for move in minus_moves:
            minus_samp = [s for s in move.results
                          if s.ensemble==self.minus_ensemble][0]
            minus_trajectory = minus_samp.trajectory
            minus_summ = minus_sides_summary(minus_trajectory,
                                             self.minus_ensemble)
            for key in self.minus_count_sides.keys():
                self.minus_count_sides[key].extend(minus_summ[key])
       
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


# TODO: move this to trajectory.summarize_volumes(label_dict)?
def summarize_trajectory_volumes(trajectory, label_dict):
    """Summarize trajectory based on number of continuous frames in volumes.

    This uses a dictionary of disjoint volumes: the volumes must be disjoint
    so that every frame can be mapped to one volume. If the frame maps to
    none of the given volumes, it returns the label None.

    Parameters
    ----------
    trajectory : Trajectory
        input trajectory
    label_dict : dict
        dictionary with labels for keys and volumes for values

    Returns
    -------
    list of tuple
        format is (label, number_of_frames)
    """
    last_vol = None
    count = 0
    segment_labels = []
    for frame in trajectory:
        in_state = []
        for key in label_dict.keys():
            vol = label_dict[key]
            if vol(frame):
                in_state.append(key)
        if len(in_state) > 1:
            raise RuntimeError("Volumes given to summarize_trajectory not disjoint")
        if len(in_state) == 0:
            current_vol = None
        else:
            current_vol = in_state[0]
        
        if last_vol == current_vol:
            count += 1
        else:
            if count > 0:
                segment_labels.append( (last_vol, count) )
            last_vol = current_vol
            count = 1
    segment_labels.append( (last_vol, count) )
    return segment_labels


# TODO: test this with manufactured summaries that have interstitials and
# multiple excursions
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
    summary = summarize_trajectory_volumes(trajectory, vol_dict)
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
