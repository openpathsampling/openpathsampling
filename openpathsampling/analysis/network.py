import openpathsampling as paths

class TransitionNetwork(object):
    def __init__(self):
        pass

    @property
    def all_ensembles(self):
        pass

#    def replica_exchange_matrix(self):

class TISNetwork(TransitionNetwork):
    # TODO: most of the analysis stuff should end up in here; the bigger
    # differences are in setup, not analysis
    def __init__(self):
        # this should check to build the replica exchange network. If the
        # number of neighbors at any station is more than 2, we can't do
        # "normal" replica flow -- instead produce a network graph. Or,
        # actually, ALWAYS produce a network graph (although this will be a
        # feature to implement later)
        pass

    def from_transitions(self, transitions, interfaces=None):
        # this will have to be disabled until I can do something
        # better with it
        pass

def join_ms_outer(outers):
    pass

def join_mis_minus(minuses):
    pass

# hack to used a named tuple for data here.. rename from DepartureInfo?
# see: http://stackoverflow.com/questions/11351032/
class DepartureInfo(collections.namedtuple('DepartureInfo',
                                           ['state', 'interfaces', 'name',
                                            'orderparameter'])):
    def __new__(cls, state, interfaces, name, orderparameter=None):
        return super(DepartureInfo, cls).__new__(cls, state, interfaces,
                                                 name, orderparameter)


class MSTISNetwork(TISNetwork):
    def __init__(self, trans_info):
        states, interfaces, names, orderparams = zip(*trans_info)
        self.from_state = {}
        msouters = []
        for state, ifaces, op, name in trans_info:
            state_index = states.index(state)
            other_states = states[:state_index]+states[state_index+1:]
            union_others = other_states[0]
            for other in other_states[1:]:
                union_others = union_others | other

            self.from_state[state] = RETISTransition(
                stateA=state, 
                stateB=union_others,
                interfaces=ifaces[:-1],
                name="Out "+name,
                orderparameter=op
            )
            msouters.append(ifaces[-1])

        pass

#    def disallow(self, stateA, stateB):

    def default_movers(self):
        pass



class MISTISNetwork(TISNetwork):
    def __init__(self, transitions):
        pass

#    def multiple_set_minus_switching(self):

    def default_movers(self):
        pass

