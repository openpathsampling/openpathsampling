import openpathsampling as paths
from openpathsampling.todict import OPSNamed

import sys

class MoveScheme(OPSNamed):
    def __init__(self, network):
        self.movers = network.movers
        self.network = network
        self._mover_acceptance = {} # used in analysis
        self.strategies = {}

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

    def move_decision_tree(self):
        for lev in sorted(self.strategies.keys()):
            for strat in self.strategies[lev]:
                strat.network = self.network
                strat.make_scheme(self)
        # return self.???


    def include_movers(self, movers, group, replace):
        if replace:
            try:
                n_existing = len(self.movers[group])
            except KeyError:
                # if the group doesn't exist, set it to these movers
                self.movers[group] = movers
            else:
                # Note that the following means that if the list of new
                # movers includes two movers with the same sig, the second
                # will overwrite the first. This is desired behavior.
                existing_sigs = {self.movers[group][i].ensemble_signature : i
                                 for i in range(n_existing)}
                # For each mover, if its signature exists in the existing
                # movers, replace the existing. Otherwise, append it to the
                # list.
                for mover in movers:
                    m_sig = mover.ensemble_signature
                    if m_sig in existing_sigs.keys():
                        self.movers[group][existing_sigs[m_sig]] = mover
                    else:
                        self.movers[group].append(mover)
        else:
            try:
                self.movers[group].extend(movers)
            except KeyError:
                self.movers[group] = movers

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

    

