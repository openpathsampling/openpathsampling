from openpathsampling.high_level.network import FixedLengthTPSNetwork
from openpathsampling.high_level.transition import FixedLengthTPSTransition
import openpathsampling as paths

class PartInBFixedLengthTPSTransition(FixedLengthTPSTransition):
    """Fixed length TPS transition accepting any frame in the final state.

    Transition that builds an ensemble used to facilitate the rate
    calculation in fixed-length TPS. [1]_ Details in
    :class:`.PartInBFixedLengthTPSNetwork`.

    See also
    --------
    PartInBFixedLengthTPSNetwork

    References
    ----------

    .. [1] C. Dellago, P.G. Bolhuis, and D. Chandler. J. Chem. Phys. 110,
           6617 (1999). http://dx.doi.org/10.1063/1.478569
    """
    def _tps_ensemble(self, stateA, stateB):
        return paths.SequentialEnsemble([
            paths.LengthEnsemble(1) & paths.AllInXEnsemble(stateA),
            paths.LengthEnsemble(self.length - 1) \
                & paths.PartInXEnsemble(stateB)
        ])

class PartInBFixedLengthTPSNetwork(FixedLengthTPSNetwork):
    """Network for fixed-length TPS accepting any frame in the final state

    This network samples a single path ensemble where the paths must begin
    in an initial state, run for a fixed total number of frames, and must
    have at least one frame in a final state. This  was used to assist in
    the flux part of the TPS rate calculation. [1]_ This version is
    generalized to multiple states.

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

    .. [1] C. Dellago, P.G. Bolhuis, and D. Chandler. J. Chem. Phys. 110,
           6617 (1999). http://dx.doi.org/10.1063/1.478569
    """
    TransitionType = PartInBFixedLengthTPSTransition
