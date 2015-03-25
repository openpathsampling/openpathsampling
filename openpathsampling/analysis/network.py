import openpathsampling as paths

class TransitionNetwork(object):
    def __init__(self):
        pass

    def replica_flow_map(self):
        pass

    def replica_flow(self, bottom, top):
        pass


class TISNetwork(TransitionNetwork):
    def __init__(self):
        pass

    def from_transitions(self, transitions, interfaces=None):
        # TODO: this will have to be disabled until I can do something
        # better with it
        pass

    def default_mover(self):
        pass


class MSTISNetwork(TISNetwork):
    def __init__(self, trans_info):
        states, interfaces, orderparams, names = zip(*trans_info)
        self.out_state = {}
        msouters = []
        for state, ifaces, op, name in trans_info:
            state_index = states.index(state)
            other_states = states[:state_index]+states[state_index+1:]
            union_others = other_states[0]
            for other in other_states[1:]:
                union_others = union_others | other

            self.out_state[state] = RETISTransition(
                stateA=state, 
                stateB=union_others,
                interfaces=ifaces[:-1],
                name="Out "+name,
                orderparameter=op
            )
            msouters.append(ifaces[-1])

        pass

    def disallow(self, stateA, stateB):

        pass

    def default_mover(self):
        pass



class MISTISNetwork(TISNetwork):
    def __init__(self):
        pass

    def default_mover(self):
        pass

