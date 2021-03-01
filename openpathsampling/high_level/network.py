import logging
import itertools

import pandas as pd

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject

# from functools import reduce  # not built-in for py3

logger = logging.getLogger(__name__)


def _default_state_name(state):
    return state.name if state.is_named else str(state)


def _name_unnamed_states(unnamed_states, all_names):
    name_index = 0
    for state in unnamed_states:
        while index_to_string(name_index) in all_names:
            name_index += 1
        state = state.named(index_to_string(name_index))
        name_index += 1


def _or_bar_namer(volumes):
    return "|".join([v.name for v in volumes])


# TODO: this should be moved into a general tools module
def listify(obj):
    try:
        _ = iter(obj)
    except TypeError:
        obj = [obj]
    return obj


def index_to_string(index):
    n_underscore = index // 26
    letter_value = index % 26
    mystr = "_"*n_underscore + chr(65 + letter_value)
    return mystr


# TODO: this will be removed when we start using the analysis.tis methods
# for the network.rate_matrix
def _set_hist_args(transition, hist_args):
    for histname in hist_args.keys():
        trans_hist = transition.ensemble_histogram_info[histname]
        if trans_hist.hist_args == {}:
            trans_hist.hist_args = hist_args[histname]


class TransitionNetwork(StorableNamedObject):
    """
    Subclasses of TransitionNetwork are the main way to set up calculations

    Attributes
    ----------
    sampling_ensembles
    all_ensembles
    sampling_transitions
    """
    def __init__(self):
        super(TransitionNetwork, self).__init__()
        # self.transitions = {}
        # self.special_ensembles = {}

    @property
    def sampling_ensembles(self):
        """
        Ensembles from the sampling transitions, excluding special ensembles.
        """
        return sum([t.ensembles for t in self.sampling_transitions], [])

    @property
    def analysis_ensembles(self):
        """
        Ensembles from the analysis transitions, excluding special ensembles.
        """
        return sum([t.ensembles for t in self.transitions.values()], [])

    @property
    def all_ensembles(self):
        """
        All ensembles in the sampling transitions, including special
        ensembles.
        """
        all_ens = self.sampling_ensembles
        for special_dict in self.special_ensembles.values():
            all_ens.extend(list(special_dict.keys()))
        return all_ens

    @property
    def sampling_transitions(self):
        """The transitions used in sampling"""
        try:
            return self._sampling_transitions
        except AttributeError:
            return None


class GeneralizedTPSNetwork(TransitionNetwork):
    """General class for TPS-based method.

    The main differences between fixed-length and flexible-length TPS is a
    small change in the ensemble. In implementation, this means that they
    use different transition classes, and that they have slightly different
    function signatures (fixed-length requires a length argument).

    To simplify this, and to make the docstrings readable, we make each
    class into a simple subclass of this GeneralizedTPSNetwork, which acts
    as an abstract class that manages most of the relevant code.

    Parameters
    ----------
    initial_states : list of :class:`.Volume`
        acceptable initial states
    final_states : list of :class:`.Volume`
        acceptable final states
    allow_self_transitions : bool
        whether self-transitions (A->A) are allowed; default is False

    Attributes
    ----------
    TransitionType : :class:`paths.Transition`
        Type of transition used here. Sets, for example, fixed or flexible
        pathlengths.
    """
    TransitionType = NotImplemented

    def __init__(self, initial_states, final_states,
                 allow_self_transitions=False, **kwargs):
        # **kwargs gets passed to the transition
        super(GeneralizedTPSNetwork, self).__init__()
        self.initial_states = listify(initial_states)
        self.final_states = listify(final_states)
        self.special_ensembles = {None: {}}

        all_initial = paths.join_volumes(self.initial_states, _or_bar_namer)

        if set(self.initial_states) == set(self.final_states):
            all_final = all_initial  # so we don't create 2 objs for it
        else:
            all_final = paths.join_volumes(self.final_states, _or_bar_namer)

        self._sampling_transitions, self.transitions = \
            self._build_transitions(self.initial_states, self.final_states,
                                    allow_self_transitions, **kwargs)

    def _build_transitions(self, initial_states, final_states,
                           allow_self_transitions, **kwargs):
        sampling_transitions = self._build_sampling_transitions(
            initial_states, final_states, allow_self_transitions, **kwargs
        )
        transitions = self._build_analysis_transitions(
            initial_states, final_states, allow_self_transitions, **kwargs
        )
        return sampling_transitions, transitions

    def _sampling_transitions_from_pairs(self, state_pairs, **kwargs):
        initial, final = state_pairs[0]
        sampling_transition = self.TransitionType(initial, final, **kwargs)
        for initial, final in state_pairs[1:]:
            sampling_transition.add_transition(initial, final)
        return [sampling_transition]

    def _build_sampling_transitions(self, initial_states, final_states,
                                    allow_self_transitions, **kwargs):
        if allow_self_transitions:
            initial_to_joined_final = {
                initial: paths.join_volumes(final_states, _or_bar_namer)
                for initial in initial_states
            }
        else:
            initial_to_joined_final = {
                initial: paths.join_volumes([final for final in final_states
                                             if initial != final],
                                            _or_bar_namer)
                for initial in initial_states
            }
        sampling_transitions = self._sampling_transitions_from_pairs(
            state_pairs=list(initial_to_joined_final.items()),
            **kwargs
        )
        return sampling_transitions

    def _build_analysis_transitions(self, initial_states, final_states,
                                    allow_self_transitions, **kwargs):
        transitions = {
            (initial, final): self.TransitionType(initial, final, **kwargs)
            for (initial, final) in itertools.product(initial_states,
                                                      final_states)
            if initial != final
        }
        return transitions

    def to_dict(self):
        ret_dict = {
            'transitions': self.transitions,
            'x_sampling_transitions': self._sampling_transitions,
            'special_ensembles': self.special_ensembles
        }
        try:
            ret_dict['initial_states'] = self.initial_states
            ret_dict['final_states'] = self.final_states
        except AttributeError:  # pragma: no cover
            # DEPRECATED: remove for 2.0
            from openpathsampling.deprecations import \
                    SAVE_RELOAD_OLD_TPS_NETWORK
            SAVE_RELOAD_OLD_TPS_NETWORK.warn()
            pass  # backward compatibility
        return ret_dict

    @property
    def all_states(self):
        """list of all initial and final states"""
        return list(set(self.initial_states + self.final_states))

    @classmethod
    def from_dict(cls, dct):
        network = cls.__new__(cls)
        super(GeneralizedTPSNetwork, network).__init__()
        network._sampling_transitions = dct['x_sampling_transitions']
        network.transitions = dct['transitions']
        try:
            network.initial_states = dct['initial_states']
            network.final_states = dct['final_states']
        except KeyError:  # pragma: no cover
            # DEPRECATED: remove for 2.0
            pass  # backward compatibility
        try:
            network.special_ensembles = dct['special_ensembles']
        except KeyError:  # pragma: no cover
            # DEPRECATED: remove for 2.0
            network.special_ensembles = {None: {}}
            # default behavior for backward compatibility
        return network

    @classmethod
    def from_state_pairs(cls, state_pairs, **kwargs):
        # TODO: redo this to use the new _sampling_transitions_from_pairs
        # method
        sampling = []
        transitions = {}
        initial_states = []
        final_states = []
        for (initial, final) in state_pairs:
            initial_states += [initial]
            final_states += [final]
            if len(sampling) == 1:
                sampling[0].add_transition(initial, final)
            elif len(sampling) == 0:
                sampling = [cls.TransitionType(initial, final, **kwargs)]
            else:
                raise RuntimeError("More than one sampling transition for TPS?")

            transitions[(initial, final)] = cls.TransitionType(initial,
                                                               final,
                                                               **kwargs)

        dict_result = {
            'x_sampling_transitions': sampling,
            'transitions': transitions
        }
        dict_result.update(kwargs)
        network = cls.from_dict(dict_result)
        network.initial_states = initial_states
        network.final_states = final_states
        return network

    @classmethod
    def from_states_all_to_all(cls, states, allow_self_transitions=False,
                               **kwargs):
        return cls(states, states,
                   allow_self_transitions=allow_self_transitions, **kwargs)


class TPSNetwork(GeneralizedTPSNetwork):
    """
    Class for flexible pathlength TPS networks (2-state or multiple state).
    """
    TransitionType = paths.TPSTransition

    # we implement these functions entirely to fix the signature (super's
    # version allow arbitrary kwargs) so the documentation can read them.
    def __init__(self, initial_states, final_states,
                 allow_self_transitions=False):
        super(TPSNetwork, self).__init__(initial_states, final_states,
                                         allow_self_transitions)

    @classmethod
    def from_state_pairs(cls, state_pairs, allow_self_transitions=False):
        return super(TPSNetwork, cls).from_state_pairs(state_pairs)

    @classmethod
    def from_states_all_to_all(cls, states, allow_self_transitions=False):
        return super(TPSNetwork, cls).from_states_all_to_all(
            states, allow_self_transitions
        )


class FixedLengthTPSNetwork(GeneralizedTPSNetwork):
    """
    Class for fixed pathlength TPS networks (2-states or multiple states).
    """
    TransitionType = paths.FixedLengthTPSTransition

    # as with TPSNetwork, we don't really need to add these functions.
    # However, without them, we need to explicitly name `length` as
    # length=value in these functions. This frees us of that, and gives us
    # clearer documentation.
    def __init__(self, initial_states, final_states, length,
                 allow_self_transitions=False):
        super(FixedLengthTPSNetwork, self).__init__(
            initial_states, final_states,
            allow_self_transitions=allow_self_transitions, length=length
        )

    @classmethod
    def from_state_pairs(cls, state_pairs, length):
        return super(FixedLengthTPSNetwork, cls).from_state_pairs(
            state_pairs, length=length
        )

    @classmethod
    def from_states_all_to_all(cls, states, length,
                               allow_self_transitions=False):
        return super(FixedLengthTPSNetwork, cls).from_states_all_to_all(
            states=states,
            allow_self_transitions=allow_self_transitions,
            length=length
        )


class TISNetwork(TransitionNetwork):
    # NOTE: this is an abstract class with several properties used by many
    # TIS-based networks
    # TODO: most of the analysis stuff should end up in here; the bigger
    # differences are in setup, not analysis
    def __init__(self, trans_info, ms_outers):
        self.trans_info = trans_info
        try:
            ms_outers = list(ms_outers)
        except TypeError:
            if ms_outers is not None:
                ms_outers = [ms_outers]
        self.ms_outer_objects = ms_outers
        self._sampling_to_analysis = None
        self._analysis_to_sampling = None
        self._sampling_ensemble_for = None
        super(TISNetwork, self).__init__()

    @property
    def sampling_to_analysis(self):
        """dict mapping sampling transitions to analysis transitions"""
        if self._sampling_to_analysis is None:
            self._sampling_to_analysis = {
                sampling_t: [t for t in self.transitions.values()
                             if sampling_t.interfaces == t.interfaces]
                for sampling_t in self.sampling_transitions
            }
        return self._sampling_to_analysis

    @property
    def analysis_to_sampling(self):
        """dict mapping analysis transitions to sampling transitions"""
        # in current examples, the result list here is always length 1, but
        # perhaps future methods will use multiple sampling transitions
        # (different order parameters?) to describe one physical transition
        if self._analysis_to_sampling is None:
            self._analysis_to_sampling = {
                t: [sampling_t for sampling_t in self.sampling_to_analysis
                    if t in self.sampling_to_analysis[sampling_t]]
                for t in self.transitions.values()
            }
        return self._analysis_to_sampling

    @property
    def sampling_ensemble_for(self):
        """dict mapping ensembles (incl. sampling) to sampling ensemble"""
        if self._sampling_ensemble_for is None:
            self._sampling_ensemble_for = {ens: ens
                                           for ens in self.sampling_ensembles}
            for ens in self.analysis_ensembles:
                analysis_transitions = [t for t in self.transitions.values()
                                        if ens in t.ensembles]
                analysis_trans = analysis_transitions[0]  # could use any
                ens_idx = analysis_trans.ensembles.index(ens)
                sampling_trans = self.analysis_to_sampling[analysis_trans]
                assert len(sampling_trans) == 1  # this only works in this case
                sampling_ens = sampling_trans[0].ensembles[ens_idx]
                self._sampling_ensemble_for[ens] = sampling_ens
        return self._sampling_ensemble_for

    def set_fluxes(self, flux_dictionary):
        """
        Parameters
        ----------
        flux_dictionary : dict of 2-tuple to float
            keys are in the form (state, interface), and values are the
            associated flux

        Raises
        ------
        KeyError
            If the flux for one of the transitions isn't in the dictionary.
        """
        # for now, if you don't have all the fluxes needed, it raises a
        # KeyError
        for trans in self.transitions.values():
            trans._flux = flux_dictionary[(trans.stateA, trans.interfaces[0])]

    @property
    def minus_ensembles(self):
        return list(self.special_ensembles['minus'].keys())

    @property
    def ms_outers(self):
        return list(self.special_ensembles['ms_outer'].keys())

    def add_ms_outer_interface(self, ms_outer, transitions, forbidden=None):
        relevant = ms_outer.relevant_transitions(transitions)
        ensemble = ms_outer.make_ensemble(relevant, forbidden)
        # TODO: this should use defaultdict, I think
        dct = {ensemble: relevant}
        try:
            self.special_ensembles['ms_outer'].update(dct)
        except KeyError:
            self.special_ensembles['ms_outer'] = dct

    @property
    def all_states(self):
        return list(set(self.initial_states + self.final_states))

    def get_state(self, snapshot):
        """
        Find which core state a snapshot is in, if any

        Parameters
        ----------
        snapshot : `openpathsampling.engines.BaseSnapshot`
            the snapshot to be tested

        Returns
        -------
        `openpathsampling.Volume`
            the volume object defining the state
        """
        for state in self.all_states:
            if state(snapshot):
                return state

        return None


class MSTISNetwork(TISNetwork):
    """
    Multiple state transition interface sampling network.

    The way this works is that it sees two effective sets of transitions.
    First, there are sampling transitions. These are based on ensembles
    which go to any final state. Second, there are analysis transitions.
    These are based on ensembles which go to a specific final state.

    Sampling is done using the sampling transitions. Sampling transitions
    are stored in the `from_state[state]` dictionary. For MSTIS, the flux
    and total crossing probabilities are independent of the final state, and
    so the analysis calculates them in the sampling transitions, and copies
    the results into the analysis transitions. This way flux and total
    crossing probably are only calculated once per interface set.

    The conditional transition probability depends on the final state, so it
    (and the rate) are calculated using the analysis transitions. The
    analysis transitions are obtained using `.transition[(stateA, stateB)]`.
    """
    def to_dict(self):
        ret_dict = {
            'from_state': self.from_state,
            'states': self.states,
            'special_ensembles': self.special_ensembles,
            'trans_info': self.trans_info,
            'ms_outer_objects': self.ms_outer_objects
        }
        return ret_dict

    @classmethod
    def from_dict(cls, dct):
        network = cls.__new__(cls)

        # replace automatically created attributes with stored ones
        network.from_state = dct['from_state']
        network.special_ensembles = dct['special_ensembles']
        network.states = dct['states']
        network.__init__(
            trans_info=dct['trans_info'],
            ms_outers=dct['ms_outer_objects']
        )
        return network

    def __init__(self, trans_info, ms_outers=None):
        """
        Creates MSTISNetwork, including interfaces.

        Parameters
        ----------
        trans_info : list of tuple
            Details of each state-based ensemble set. 2-tuple in the order
            (state, interface_set) where state is a Volume, and
            interface_set is an InterfaceSet (with associated
            CollectiveVariable)
        ms_outers : MSOuterTISInterface or list of MSOuterTISInterface
            mutliple state outer interfaces for this network
        """
        super(MSTISNetwork, self).__init__(trans_info, ms_outers)
        # build sampling transitions
        states, interfaces = zip(*trans_info)
        self.states = states
        if not hasattr(self, "from_state"):
            self.special_ensembles = {}
            self.from_state = {}
            self._build_fromstate_transitions(trans_info)
            if self.ms_outer_objects is not None:
                for ms_outer in self.ms_outer_objects:
                    all_transitions = list(self.from_state.values())
                    self.add_ms_outer_interface(ms_outer, all_transitions)

        self._sampling_transitions = list(self.from_state.values())

        # by default, we set assign these values to all ensembles
        self.hist_args = {}

        self.transitions = self._build_analysis_transitions()

    @property
    def all_states(self):
        return self.states

    def _build_transitions(self, trans_info, ms_outers, special_ensembles):
        sampling_ensembles = self._build_sampling_ensembles(trans_info)


        return sampling_transitions, transitions, special_ensembles

    @staticmethod
    def _build_analysis_transition_for_sampling(sampling_transition,
                                                all_states):
        local_transitions = {}
        state_A = sampling_transition.stateA
        other_states = set(all_states) - set([state_A])
        str_A = _default_state_name(state_A)
        for state_B in other_states:
            str_B = _default_state_name(state_B)
            trans = paths.TISTransition(
                stateA=state_A,
                stateB=state_B,
                interfaces=sampling_transition.interfaces,
                name=str_A + "->" + str_B,
                orderparameter=sampling_transition.orderparameter
            )
            # override created stuff
            trans.ensembles = sampling_transition.ensembles
            for i in range(len(trans.ensembles)):
                trans.ensembles[i].named(trans.name + "[" + str(i) + "]")

            trans.minus_ensemble = sampling_transition.minus_ensemble
            local_transitions[(state_A, state_B)] = trans
        return local_transitions

    def _build_analysis_transitions(self):
        # set up analysis transitions (not to be saved)
        transitions = {}
        for from_A in self.from_state.values():
            local_transitions = self._build_analysis_transition_for_sampling(
                sampling_transition=from_A,
                all_states=self.all_states
            )
            transitions.update(local_transitions)

        return transitions

    @staticmethod
    def build_one_state_sampling_transition(state, interfaces, all_states):
        other_states = list(set(all_states) - set([state]))
        union_others = paths.join_volumes(
            volume_list=other_states,
            name="all states except " + str(state.name)
        )
        this_trans = paths.TISTransition(
            stateA=state,
            stateB=union_others,
            interfaces=interfaces,
            name="Out " + state.name,
            orderparameter=interfaces.cv
        )
        return this_trans

    def _build_fromstate_transitions(self, trans_info):
        """
        Builds the sampling transitions (the self.from_state dictionary).

        This also sets self.states (list of states volumes), self.outers
        (list of interface volumes making the MS-outer interface), and
        self.outer_ensembles (list of TISEnsembles associated with the
        self.outers interfaces). Additionally, it gives default names
        volumes, interfaces, and transitions.

        Parameters
        ----------
        trans_info : list of 2-tuples
            See description in __init__.

        """
        states, interfaces = zip(*trans_info)
        orderparams = [iface_set.cv for iface_set in interfaces]

        # NAMING STATES (give default names)
        all_states = paths.join_volumes(states).named("all states")
        all_names = list(set([s.name for s in states]))
        unnamed_states = [s for s in states if not s.is_named]
        _name_unnamed_states(unnamed_states, all_names)

        # BUILDING ENSEMBLES
        self.states = states
        for (state, ifaces) in trans_info:
            this_trans = self.build_one_state_sampling_transition(
                state=state,
                interfaces=ifaces,
                all_states=states
            )
            # op = ifaces.cv
            # state_index = states.index(state)
            # other_states = states[:state_index]+states[state_index+1:]
            # other_states = list(set(states) - set([state]))
            # union_others = paths.join_volumes(
                # volume_list=other_states,
                # name="all states except " + str(state.name)
            # )
            # union_others = paths.volume.join_volumes(other_states)
            # union_others.named("all states except " + str(state.name))
            # out_others = paths.AllOutXEnsemble(union_others)

            # this_trans = paths.TISTransition(
                # stateA=state,
                # stateB=union_others,
                # interfaces=ifaces,
                # name="Out " + state.name,
                # orderparameter=op
            # )

            self.from_state[state] = this_trans

            this_minus = self.from_state[state].minus_ensemble #& out_others
            this_inner = self.from_state[state].ensembles[0]
            # TODO: this should use defaultdict, I think
            try:
                self.special_ensembles['minus'][this_minus] = [this_trans]
            except KeyError:
                self.special_ensembles['minus'] = {this_minus : [this_trans]}

    def __str__(self):
        mystr = "Multiple State TIS Network:\n"
        for state in self.from_state.keys():
            mystr += str(self.from_state[state])
        return mystr

    def rate_matrix(self, steps, force=False):
        """
        Calculate the matrix of all rates.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            steps to be analyzed
        force : bool (False)
            if True, cached results are overwritten

        Returns
        -------
        pandas.DataFrame
            Rates from row_label to column_label. Diagonal is NaN.
        """
        # for each transition in from_state:
        # 1. Calculate the flux and the TCP
        names = [s.name for s in self.states]
        self._rate_matrix = pd.DataFrame(columns=names, index=names)
        for stateA in self.from_state.keys():
            transition = self.from_state[stateA]
            # set up the hist_args if necessary
            _set_hist_args(transition, self.hist_args)
            # for histname in self.hist_args.keys():
                # trans_hist = transition.ensemble_histogram_info[histname]
                # if trans_hist.hist_args == {}:
                    # trans_hist.hist_args = self.hist_args[histname]

            transition.total_crossing_probability(steps=steps,
                                                  force=force)
            transition._minus_move_flux(steps=steps, force=force)
            for stateB in self.from_state.keys():
                if stateA != stateB:
                    analysis_trans = self.transitions[(stateA, stateB)]
                    analysis_trans.copy_analysis_from(transition)

        for trans in self.transitions.values():
            rate = trans.rate(steps)
            # self._rate_matrix.set_value(trans.stateA.name,
                                        # trans.stateB.name,
                                        # rate)
            name_A = trans.stateA.name
            name_B = trans.stateB.name
            self._rate_matrix.at[name_A, name_B] = rate
            #print trans.stateA.name, trans.stateB.name,
            #print rate

        return self._rate_matrix


class MISTISNetwork(TISNetwork):
    """
    Multiple interface set TIS network.

    Input is given as a list of 4-tuples. Each 4-tuple represents a
    transition, and is in the order::

        (initial_state, interfaces, order_parameter, final_states)

    This will create the `input_transitions` objects.

    Attributes
    ----------
    input_transitions : list of TISTransition
        the transitions given as input
    sampling_transitions : list of TISTransition
        the transitions used in sampling
    transitions : list of TISTransition
        the transitions used in analysis

    Note
    ----
        The distinction between the three types of transitions in the object
        are a bit subtle, but important. The `input_transitions` are, of
        course, the transitions given in the input. These are A->B
        transitions, but would allow any other state. The
        `sampling_transitions` are what are used in sampling. These are
        A->any transitions if strict sampling is off, or "A->B & not_others"
        if strict sampling is on. Finally, the regular `transitions` are the
        transitions that are used for analysis (use the sampling ensembles
        for the interfaces, but also A->B).

    Parameters
    ----------
    trans_info : list of tuple
        Details of each interface set. 3-tuple in the order (initial_state,
        interfaces, final_state) where initial_state and final_state are
        Volumes, and interfaces is an InterfaceSet
    ms_outers : MSOuterTISInterface or list of MSOuterTISInterface
        mutliple state outer interfaces for this network
    strict_sampling : bool
        whether the final state from the tuple is the *only* allowed final
        state in the sampling; default False
    """
    # NOTE: input_transitions are in addition to the sampling_transitions
    # and the transitions (analysis transitions)
    def __init__(self, trans_info, ms_outers=None, strict_sampling=False):
        super(MISTISNetwork, self).__init__(trans_info, ms_outers)
        self.strict_sampling = strict_sampling
        states_A, interfaces, states_B = zip(*trans_info)
        orderparams = [iface_set.cv for iface_set in interfaces]
        self.initial_states = list(set(states_A))
        self.final_states = list(set(states_B))
        list_all_states = list(set(self.initial_states + self.final_states))

        # name states
        all_state_names = list(set([s.name for s in list_all_states]))
        unnamed_states = [s for s in list_all_states if not s.is_named]
        name_index = 0
        for state in unnamed_states:
            while index_to_string(name_index) in all_state_names:
                name_index += 1
            state.named(index_to_string(name_index))
            name_index += 1

        if not hasattr(self, "input_transitions"):
            self.input_transitions = {
                (stateA, stateB):
                paths.TISTransition(stateA, stateB, interface, interface.cv,
                                    name=stateA.name + "->" + stateB.name,
                                    name_suffix=" (input)")
                for (stateA, interface, stateB) in self.trans_info
            }

        if not hasattr(self, 'x_sampling_transitions'):
            self.special_ensembles = {}
            self._build_sampling_transitions(self.input_transitions.values())
            if self.ms_outer_objects is not None:
                for ms_outer in self.ms_outer_objects:
                    all_transitions = self.x_sampling_transitions
                    if not self.strict_sampling:
                        self.add_ms_outer_interface(ms_outer, all_transitions)
                    else:
                        relevant = ms_outer.relevant_transitions(all_transitions)
                        allowed = set(sum([[t.stateA, t.stateB]
                                           for t in relevant], []))
                        forbidden = set(list_all_states) - allowed
                        self.add_ms_outer_interface(ms_outer,
                                                    all_transitions,
                                                    forbidden)

        self._sampling_transitions = self.x_sampling_transitions

        # by default, we set assign these values to all ensembles
        self.hist_args = {}

        self._build_analysis_transitions()

    def to_dict(self):
        ret_dict = {
            'special_ensembles': self.special_ensembles,
            'transition_pairs': self.transition_pairs,
            'x_sampling_transitions': self.x_sampling_transitions,
            'transition_to_sampling': self.transition_to_sampling,
            'input_transitions': self.input_transitions,
            'trans_info': self.trans_info,
            'strict_sampling': self.strict_sampling,
            'ms_outer_objects': self.ms_outer_objects
        }
        return ret_dict

    @staticmethod
    def from_dict(dct):
        network = MISTISNetwork.__new__(MISTISNetwork)
        network.special_ensembles = dct['special_ensembles']
        network.transition_pairs = dct['transition_pairs']
        network.transition_to_sampling = dct['transition_to_sampling']
        network.input_transitions = dct['input_transitions']
        network.x_sampling_transitions = dct['x_sampling_transitions']
        network.__init__(trans_info=dct['trans_info'],
                         ms_outers=dct['ms_outer_objects'],
                         strict_sampling=dct['strict_sampling'])
        return network

    def _build_transition_pairs(self, transitions):
        # identify transition pairs
        transition_pair_set_dict = {}
        for initial in self.initial_states:
            for t1 in [t for t in transitions if t.stateA == initial]:
                t_reverse = [
                    t for t in transitions
                    if t.stateA == t1.stateB and t.stateB == t1.stateA
                ]
                if len(t_reverse) == 1:
                    key = frozenset([t1.stateA, t1.stateB])
                    new_v = [t1, t_reverse[0]]
                    if key not in transition_pair_set_dict.keys():
                        transition_pair_set_dict[key] = new_v
                elif len(t_reverse) > 1:  # pragma: no cover
                    raise RuntimeError("More than one reverse transition")
                # if len(t_reverse) is 0, we just pass

        transition_pairs = list(transition_pair_set_dict.values())
        return transition_pairs

    def _build_sampling_transitions(self, transitions):
        transitions = list(transitions)  # input may be iterator
        # TODO: I don't think transition pairs are used (see comment below;
        # I think that was the previous use case -- as input to all_in_pairs
        # However, existing files store this, so we won't actually remove it
        # yet.
        self.transition_pairs = self._build_transition_pairs(transitions)
        # this seems to no longer be used; I think it was necessary when the
        # MSOuter interface was done implicitly, instead of explicitly. Then
        # we turn the outermost to MS if and only if it was paired with the
        # reverse transition
        # if len(self.transition_pairs) > 0:
            # all_in_pairs = reduce(list.__add__, map(lambda x: list(x),
                                                    # self.transition_pairs))
        # else:
            # all_in_pairs = []

        # build sampling transitions
        all_states = paths.join_volumes(self.initial_states + self.final_states)
        all_states_set = set(self.initial_states + self.final_states)
        self.transition_to_sampling = {}
        for transition in transitions:
            stateA = transition.stateA
            stateB = transition.stateB
            if self.strict_sampling:
                final_state = stateB
                other_states = paths.join_volumes(all_states_set -
                                                  set([stateA, stateB]))
                ensemble_to_intersect = paths.AllOutXEnsemble(other_states)
            else:
                final_state = paths.join_volumes(all_states_set)
                ensemble_to_intersect = paths.FullEnsemble()

            sample_trans = paths.TISTransition(
                stateA=stateA,
                stateB=final_state,
                interfaces=transition.interfaces,
                name=stateA.name + "->" + stateB.name,
                orderparameter=transition.orderparameter
            )

            new_ensembles = [e & ensemble_to_intersect
                             for e in sample_trans.ensembles]
            if self.strict_sampling:
                for (old, new) in zip(new_ensembles, sample_trans.ensembles):
                    old.name = new.name + " strict"
            sample_trans.ensembles = new_ensembles
            sample_trans.named("Sampling " + str(stateA) + "->" + str(stateB))
            self.transition_to_sampling[transition] = sample_trans

        self.x_sampling_transitions = \
                list(self.transition_to_sampling.values())

        self._build_sampling_minus_ensembles()

    def _build_sampling_minus_ensembles(self):
        # combining the minus interfaces
        all_states_set = set(self.initial_states + self.final_states)
        for initial in self.initial_states:
            innermosts = []
            # trans_from_initial: list of transition from initial
            trans_from_initial = [
                t for t in self.x_sampling_transitions
                if t.stateA == initial
            ]
            for t1 in trans_from_initial:
                innermosts.append(t1.interfaces[0])

            forbidden = list(all_states_set - {initial})
            minus = paths.MinusInterfaceEnsemble(
                state_vol=initial,
                innermost_vols=innermosts,
                forbidden=forbidden
            ).named(initial.name + " MIS minus")
            try:
                self.special_ensembles['minus'][minus] = trans_from_initial
            except KeyError:
                self.special_ensembles['minus'] = {minus: trans_from_initial}

    def _build_analysis_transitions(self):
        self.transitions = {}
        for trans in self.input_transitions.values():
            sample_trans = self.transition_to_sampling[trans]
            stateA = trans.stateA
            stateB = trans.stateB
            analysis_trans = paths.TISTransition(
                stateA=stateA,
                stateB=stateB,
                interfaces=sample_trans.interfaces,
                orderparameter=sample_trans.orderparameter
            )
            analysis_trans.ensembles = sample_trans.ensembles
            analysis_trans.named(trans.name)
            #analysis_trans.special_ensembles = sample_trans.special_ensembles
            self.transitions[(stateA, stateB)] = analysis_trans

    def rate_matrix(self, steps, force=False):
        initial_names = [s.name for s in self.initial_states]
        final_names = [s.name for s in self.final_states]
        self._rate_matrix = pd.DataFrame(columns=final_names,
                                         index=initial_names)
        for trans in self.transitions.values():
            # set up the hist_args if necessary
            _set_hist_args(trans, self.hist_args)
            # for histname in self.hist_args.keys():
                # trans_hist = trans.ensemble_histogram_info[histname]
                # if trans_hist.hist_args == {}:
                    # trans_hist.hist_args = self.hist_args[histname]
            tcp = trans.total_crossing_probability(steps=steps,
                                                   force=force)
            if trans._flux is None:
                logger.warning("No flux for transition " + str(trans.name)
                               + ": Rate will be NaN")
                trans._flux = float("nan")
                # we give NaN so we can calculate the condition transition
                # probability automatically

            rate = trans.rate(steps)
            # self._rate_matrix.set_value(trans.stateA.name,
                                        # trans.stateB.name,
                                        # rate)
            name_A = trans.stateA.name
            name_B = trans.stateB.name
            self._rate_matrix.at[name_A, name_B] = rate

        return self._rate_matrix
