import openpathsampling as paths
import openpathsampling.netcdfplus as netcdfplus

class MSOuterTISInterface(netcdfplus.StorableNamedObject):
    def __init__(self, interface_sets, volumes, lambdas=None):
        super(MSOuterTISInterface, self).__init__()
        self.volumes = volumes
        self.interface_sets = interface_sets
        self.lambdas = lambdas
        self._interface_set_to_volume = {
            i_set: vol for (i_set, vol) in zip(interface_sets, volumes)
        }
        self._interface_set_to_lambda = {
            i_set: lmbda for (i_set, lmbda) in zip(interface_sets, lambdas)
        }

    def volume_for_interface_set(self, interface_set):
        return self._interface_set_to_volume[interface_set]

    def lambda_for_interface_set(self, interface_set):
        return self._interface_set_to_lambda[interface_set]

    @staticmethod
    def from_lambdas(interface_sets_lambdas):
        interface_sets = interface_sets_lambdas.keys()
        lambdas = interface_sets_lambdas.values()
        volumes = [iface_set.new_interface(interface_sets_lambdas[iface_set])
                   for iface_set in interface_sets]
        return MSOuterTISInterface(interface_sets, volumes, lambdas)

    def make_ensemble(self, network, forbidden=None):
        if forbidden is None:
            ensemble_to_intersect = paths.FullEnsemble()
        else:
            ensemble_to_intersect = paths.AllOutXEnsemble(forbidden)
        relevant_transitions = [t for t in network.sampling_transitions
                                if t.interfaces == self.interface_sets]

        outer_ensembles = []
        for trans in relevant_transitions:
            initial = trans.stateA
            final = trans.stateB
            volume = self.volume_for_interface_set(trans.interfaces)
            outer_ensembles.append(
                ensemble_to_intersect & paths.TISEnsemble(initial, final,
                                                          interface)
            )

        return paths.join_ensembles(outer_ensembles)


            


