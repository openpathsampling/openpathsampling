from analysis.move_scheme import (
    MoveScheme, DefaultScheme, LockedMoveScheme, OneWayShootingMoveScheme
)

from analysis.tis_analysis import (
    TISTransition, Transition, TPSTransition, FixedLengthTPSTransition
)

from analysis.network import (
    MSTISNetwork, TransitionNetwork, MISTISNetwork, TPSNetwork,
    FixedLengthTPSNetwork
)

from analysis.replica_network import (
    ReplicaNetwork, trace_ensembles_for_replica,
    trace_replicas_for_ensemble, condense_repeats,
    ReplicaNetworkGraph
)

from analysis.shooting_point_analysis import (
    ShootingPointAnalysis, SnapshotByCoordinateDict
)

from analysis.single_trajectory_analysis import (
    SingleTrajectoryAnalysis,
    TrajectorySegmentContainer
)

from collectivevariable import CV_Function, CV_MDTraj_Function, CV_MSMB_Featurizer, \
    CV_Volume, CollectiveVariable

from ensemble import (
    Ensemble, EnsembleCombination, EnsembleFactory, EntersXEnsemble,
    EmptyEnsemble, ExitsXEnsemble, FullEnsemble, PartInXEnsemble,
    AllInXEnsemble, AllOutXEnsemble, WrappedEnsemble,
    SuffixTrajectoryEnsemble, PrefixTrajectoryEnsemble,
    PartOutXEnsemble, LengthEnsemble, NegatedEnsemble,
    ReversedTrajectoryEnsemble, SequentialEnsemble, VolumeEnsemble,
    SequentialEnsemble, IntersectionEnsemble, UnionEnsemble,
    SymmetricDifferenceEnsemble, RelativeComplementEnsemble,
    SingleFrameEnsemble, MinusInterfaceEnsemble, TISEnsemble,
    OptionalEnsemble, join_ensembles
)

from live_visualization import LiveVisualization

from pathmovechange import (
    EmptyPathMoveChange, ConditionalSequentialPathMoveChange,
    PathMoveChange, PartialAcceptanceSequentialPathMoveChange,
    RandomChoicePathMoveChange, SamplePathMoveChange,
    SequentialPathMoveChange, KeepLastSamplePathMoveChange,
    FilterSamplesPathMoveChange,
    PathSimulatorPathMoveChange, AcceptedSamplePathMoveChange,
    RejectedSamplePathMoveChange, SubPathMoveChange,
    FilterByEnsemblePathMoveChange
)

from pathmover import Details, MoveDetails, SampleDetails

from pathmover import (
    RandomChoiceMover, PathMover, ConditionalSequentialMover,
    PartialAcceptanceSequentialMover, BackwardShootMover, ForwardShootMover,
    BackwardExtendMover, ForwardExtendMover, MinusMover,
    SingleReplicaMinusMover, PathMoverFactory, PathReversalMover,
    ReplicaExchangeMover, EnsembleHopMover, ReplicaIDChangeMover,
    SequentialMover, ConditionalMover,
    PathSimulatorMover, PathReversalSet, NeighborEnsembleReplicaExchange,
    SampleMover, StateSwapMover, FinalSubtrajectorySelectMover, EngineMover,
    FirstSubtrajectorySelectMover, MultipleSetMinusMover,
    OneWayShootingMover, RandomSubtrajectorySelectMover, SubPathMover,
    EnsembleFilterMover, SelectionMover, FirstAllowedMover,
    LastAllowedMover, OneWayExtendMover, SubtrajectorySelectMover
)

from pathsimulator import (
    PathSimulator, FullBootstrapping, Bootstrapping, PathSampling, MCStep,
    CommittorSimulation, DirectSimulation
)

from sample import Sample, SampleSet

from shooting import ShootingPointSelector, UniformSelector, \
    GaussianBiasSelector, FirstFrameSelector, FinalFrameSelector

from snapshot_modifier import NoModification, RandomVelocities

from storage.storage import Storage, AnalysisStorage

from volume import (
    Volume, VolumeCombination, VolumeFactory, VoronoiVolume,
    EmptyVolume, FullVolume, CVRangeVolume, CVRangeVolumePeriodic,
    IntersectionVolume, UnionVolume, SymmetricDifferenceVolume,
    RelativeComplementVolume, join_volumes
)

from openpathsampling.engines import Trajectory, BaseSnapshot
import openpathsampling.engines.openmm as openmm
import openpathsampling.engines.toy as toy


def git_HEAD():  # pragma: no cover
    from subprocess import check_output
    import os.path
    git_dir = os.path.dirname(os.path.realpath(__file__))
    return check_output(["git", "-C", git_dir, "rev-parse", "HEAD"])[:-1]
    # chops the newline at the end

try:
    import version
except ImportError:  # pragma: no cover
    pass
