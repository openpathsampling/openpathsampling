from openpathsampling.high_level.network import FixedLengthTPSNetwork
from openpathsampling.high_level.transition import FixedLengthTPSTransition
import openpathsampling as paths

class NewTransition(FixedLengthTPSTransition):
    """Fixed length TPS transition accepting any frame in the final state

    See also
    --------
    NewNetwork

    References
    ----------

    C. Dellago, P.G. Bolhuis, and D. Chandler. J. Chem. Phys. 110, 6617
    (1999). http://dx.doi.org/10.1063/1.478569
    """
    def _tps_ensemble(self, stateA, stateB):
        return paths.SequentialEnsemble([
            paths.LengthEnsemble(1) & paths.AllInXEnsemble(stateA),
            paths.LengthEnsemble(self.length - 1) \
                & paths.PartInXEnsemble(stateB)
        ])

class NewNetwork(FixedLengthTPSNetwork):
    """Network for fixed-length TPS accepting any frame in the final state

    ???

    Parameters
    ----------
    intial_states : (list of) :class:`.Volume`
        acceptable initial states
    final_states : (list of) :class:`.Volume`
        acceptable final states
    length : int
        length of paths in the path ensemble, in frames
    allow_self_transitions : bool
        whether self-transitions (A->A) are allowed; default is False. For
        this network, A->B->A transitions are *always* allowed.

    References
    ----------

    C. Dellago, P.G. Bolhuis, and D. Chandler. J. Chem. Phys. 110, 6617
    (1999). http://dx.doi.org/10.1063/1.478569
    """
    TransitionType = NewTransition
