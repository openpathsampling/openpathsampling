import openpathsampling as paths
from openpathsampling.todict import OPSNamed

from openpathsampling.analysis.move_strategy import levels as strategy_levels
import openpathsampling.analysis.move_strategy as strategies

import sys


class MoveScheme(OPSNamed):
    def __init__(self, network):
        self.movers = {}
        self.movers = network.movers # TODO: legacy
        self.network = network
        self._mover_acceptance = {} # used in analysis
        self.strategies = {}
        self.balance_partners = {}
        self.mover_weights = {}
        self.ensemble_weights = {}
        self.choice_probability = {}
        self.root_mover = None

    def append(self, strategies, levels=None):
        """
        Adds new strategies to this scheme, organized by `level`.

        Parameters
        ----------
        strategies : MoveStrategy or list of MoveStrategy
            strategies to add to this scheme
        levels : integer or list of integer or None
            levels to associate with each strategy. If None, strategy.level.
        """
        # first we clean up the input: strategies is a list of MoveStrategy;
        # levels is a list of integers
        try:
            strategies = list(strategies)
        except TypeError:
            strategies = [strategies]

        if levels is not None:
            try:
                levels = list(levels)
            except TypeError:
                levels = [levels]*len(strategies)
        else:
            levels = []
            for strat in strategies:
                levels.append(strat.level)
        
        # now we put everything into appropriate dictionaries
        for strat, lev in zip(strategies, levels):
            try:
                self.strategies[lev].append(strat)
            except KeyError:
                self.strategies[lev] = [strat]

    # TODO: it might be nice to have a way to "lock" this once it has been
    # saved. That would prevent a (stupid) user from trying to rebuild a
    # custom-modified tree
    def build_move_decision_tree(self):
        for lev in sorted(self.strategies.keys()):
            for strat in self.strategies[lev]:
                self.apply_strategy(strat)

    # TODO: should I make this a property? make root_mover into
    # _move_decision_tree? allow the user to directly set it?
    def move_decision_tree(self, rebuild=False):
        if self.root_mover is None:
            rebuild=True
        if rebuild:
            self.build_move_decision_tree()
        return self.root_mover

    def apply_strategy(self, strategy):
        movers = strategy.make_movers(self)
        group = strategy.group
        if strategy_levels.level_type(strategy.level) == strategy_levels.GLOBAL:
            # shortcut out for the global-level stuff
            self.root_mover = movers
        elif strategy.replace_signatures:
            self.movers[group] = movers
        elif strategy.replace_movers:
            try:
                n_existing = len(self.movers[group])
            except KeyError:
                # if the group doesn't exist, set it to these movers
                self.movers[group] = movers
            else:
                # Note that the following means that if the list of new
                # movers includes two movers with the same sig, the second
                # will overwrite the first. This is desired behavior. On the
                # other hand, if the list of old movers in the group already
                # has two movers with the same signature, then both should
                # be overwritten.
                existing_sigs = {}
                for i in range(n_existing):
                    key = self.movers[group][i].ensemble_signature
                    try:
                        existing_sigs[key].append(i)
                    except KeyError:
                        existing_sigs[key] = [i]

                # For each mover, if its signature exists in the existing
                # movers, replace the existing. Otherwise, append it to the
                # list.
                for mover in movers:
                    m_sig = mover.ensemble_signature
                    if m_sig in existing_sigs.keys():
                        for idx in existing_sigs[m_sig]:
                            self.movers[group][idx] = mover
                    else:
                        self.movers[group].append(mover)
        else:
            try:
                self.movers[group].extend(movers)
            except KeyError:
                self.movers[group] = movers

    def ensembles_for_move_tree(self, root=None):
        if root is None:
            root = self.root_mover
        movers = root.map_pre_order(lambda x : x)
        mover_ensemble_dict = {}
        for m in movers:
            input_sig = m.input_ensembles
            output_sig = m.output_ensembles
            for ens in input_sig + output_sig:
                mover_ensemble_dict[ens] = 1
        mover_ensembles = mover_ensemble_dict.keys()
        return mover_ensembles

    def find_hidden_ensembles(self, root=None):
        unhidden_ensembles = set(self.network.all_ensembles)
        mover_ensembles = set(self.ensembles_for_move_tree(root))
        hidden_ensembles = mover_ensembles - unhidden_ensembles
        return hidden_ensembles

    def find_unused_ensembles(self, root=None):
        unhidden_ensembles = set(self.network.all_ensembles)
        mover_ensembles = set(self.ensembles_for_move_tree(root))
        unused_ensembles = unhidden_ensembles - mover_ensembles
        return unused_ensembles

    def check_for_root(self, fcn_name):
        if self.root_mover is None:
            warnstr = ("Can't use {fcn_name} before building the move " +
                       "decision tree").format(fcn_name=fcn_name)
            raise RuntimeWarning(warn_str)

    def build_balance_partners(self):
        self.check_for_root("build_balance_partners")
        for groupname in self.movers.keys():
            group = self.movers[groupname]
            self.balance_partners[groupname] = {}
            for mover in group:
                partner_sig_set = (set(mover.output_ensembles), 
                                   set(mover.input_ensembles))
                partners = [m for m in group 
                            if m.ensemble_signature_set==partner_sig_set]
                self.balance_partners[groupname][mover] = partners

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
                            if k[0] is mover)

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


class DefaultScheme(MoveScheme):
    """
    Just a MoveScheme with the full set of default strategies: nearest
    neighbor repex, uniform selection one-way shooting, minus move, and
    path reversals, all structured as choose move type then choose specific
    move.
    """
    def __init__(self, network):
        super(DefaultScheme, self).__init__(network)
        self.append(strategies.NearestNeighborRepExStrategy())
        self.append(strategies.OneWayShootingStrategy())
        self.append(strategies.PathReversalStrategy())
        self.append(strategies.DefaultStrategy())
        self.append(strategies.MinusMoveStrategy())

        msouters = self.network.special_ensembles['ms_outer']
        for ms in msouters.keys():
            self.append(strategies.OneWayShootingStrategy(
                ensembles=[ms],
                group="ms_outer_shooting"
            ))
            self.append(strategies.PathReversalStrategy(
                ensembles=[ms],
                replace=False
            ))
            ms_neighbors = [t.ensembles[-1] for t in msouters[ms]]
            pairs = [[ms, neighb] for neighb in ms_neighbors]
            self.append(strategies.SelectedPairsRepExStrategy(
                ensembles=pairs
            ))

