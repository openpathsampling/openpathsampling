import openpathsampling as paths
import openpathsampling.netcdfplus as netcdfplus

class MSOuterTISInterface(netcdfplus.StorableNamedObject):
    """
    Object to manage multiple-state outer interface in TIS.

    The MS outer interface needs to know what the interface volume is for
    each interface set, and needs to know which interface sets it is
    associated with.

    Can also be initialized from a dictionary of interface sets to lambda
    values (assuming those interface sets have a .new_interface method). See
    :meth:`.from_lambdas`.

    This also includes a convenience to create the appropriate MS outer
    ensemble.

    Much of this currently depends on the idea that there is a one-to-one
    mapping from :class:`.Transition` to :class:`.InterfaceSet` (at least,
    among transitions used for sampling). So we use InterfaceSets, which the
    user must interact with anyway, in place of Transitions, which the use
    need not interact with.

    Parameters
    ----------
    interface_sets : list of :class:`.InterfaceSet`
        Interface sets associated with the volumes which will be used as
        interfaces in this MS-outer. The specific volume for the MS-outer is
        not included.
    volumes : list of :class:`.Volume`
        Volumes for the outermost interface (to be made part of the MS-outer
        ensemble) for each InterfaceSet in `interface_sets`. Must be the
        same length as `interface_sets`.
    lambdas : list of float or list of int
        Values of the CV associated with each volume. Optional, but
        recommended to facilitate other setup and analysis code.
    """
    def __init__(self, interface_sets, volumes, lambdas=None):
        super(MSOuterTISInterface, self).__init__()
        self.volumes = volumes
        self.interface_sets = interface_sets
        if lambdas is None:
            lambdas = [None]*len(volumes)
        self.lambdas = lambdas
        self._interface_set_to_volume = {
            i_set: vol for (i_set, vol) in zip(interface_sets, volumes)
        }
        self._interface_set_to_lambda = {
            i_set: lmbda for (i_set, lmbda) in zip(interface_sets, lambdas)
        }

    def volume_for_interface_set(self, interface_set):
        """
        Outer volume in this MS outer interface for given interface set

        Parameters
        ----------
        interface_set : :class:`.InterfaceSet`
            the interface set in question

        Returns
        -------
        :class:`.Volume`
            the outer volume associated with the input interface set
            according to this MS outer interface
        """
        return self._interface_set_to_volume[interface_set]

    def lambda_for_interface_set(self, interface_set):
        """
        CV value for outer volume associated with given interface set

        Parameters
        ----------
        interface_set : :class:`.InterfaceSet`
            the interface set in question

        Returns
        -------
        float or int
            the CV value (lambda) associated with the input interface set
            according to this MS outer interface

        """
        return self._interface_set_to_lambda[interface_set]

    @staticmethod
    def from_lambdas(interface_sets_lambdas):
        """
        Create an MSOuterTISInterface from desired lambda values.

        Parameters
        ----------
        interface_sets_lambdas : dict of {InterfaceSet: float or int}
            Dictionary mapping interface sets to the value of the collective
            variable associated with the volume for the MS outer interface.
            Requires that the interface sets have a working
            `.new_interface(lambda_i)` function.

        Returns
        -------
        MSOuterTISInterface :
            the desired object with correct volumes and lambdas
        """
        interface_sets = list(interface_sets_lambdas.keys())
        lambdas = list(interface_sets_lambdas.values())
        volumes = [iface_set.new_interface(interface_sets_lambdas[iface_set])
                   for iface_set in interface_sets]
        return MSOuterTISInterface(interface_sets, volumes, lambdas)

    def relevant_transitions(self, transitions):
        """
        Identify transitions which are relevant to this MS outer interface.

        Useful for networks, which need to know which transitions a given MS
        outer interface links.

        Parameters
        ----------
        transitions : list of :class:`.TISTransition`
            possible transitions of relevance

        Returns
        -------
        list of :class:`.TISTransitions`
            the transitions which have the same interface sets as this MS
            outer interface

        """
        return [t for t in transitions if t.interfaces in self.interface_sets]

    def make_ensemble(self, transitions, forbidden=None):
        """
        Create the ensemble for this MS outer interface.

        Parameters
        ----------
        transitions : list of :class:`.TISTransition`
            possible transitions of relevance
        forbidden : list of :class:`.Volume` or None
            (optional) volumes to disallow from the ensemble (e.g.,
            other states that should cause the trajectory to stop)

        Returns
        -------
        :class:`.Ensemble`
            the union of the TISEnsembles for each volume of the MS outer
            interface
        """
        if forbidden is None:
            ensemble_to_intersect = paths.FullEnsemble()
        else:
            try:
                _ = len(forbidden)
            except TypeError:
                forbidden = [forbidden]
            forbidden_vol = paths.join_volumes(forbidden)
            ensemble_to_intersect = paths.AllOutXEnsemble(forbidden_vol)

        # TODO: maybe we should crash if given transitions aren't relevant?
        relevant_transitions = self.relevant_transitions(transitions)

        outer_ensembles = []
        for trans in relevant_transitions:
            initial = trans.stateA
            final = trans.stateB
            volume = self.volume_for_interface_set(trans.interfaces)
            # TODO: move following to a logger.debug
            #print initial.name, final.name,\
                    #self.lambda_for_interface_set(trans.interfaces)
            outer_ensembles.append(
                ensemble_to_intersect & paths.TISEnsemble(initial, final,
                                                          volume)
            )

        return paths.join_ensembles(outer_ensembles)
