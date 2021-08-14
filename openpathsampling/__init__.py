# start by importing version, because for some weird reason, I sometimes get
# networkx as over-riding the version here (no idea why...)
try:
    # should work if installed through normal means: setup.py-based with
    # pip, conda, easy_install, etc.
    from . import version
except ImportError:  # pragma: no cover
    import os
    # should work if someone just set the $PYTHONPATH to include OPS
    directory = os.path.dirname(os.path.realpath(__file__))
    prev_dir = os.path.split(directory)[0]
    setupfile = os.path.join(prev_dir, "setup.py")

    if not os.path.exists(setupfile):
        # now we're screwed
        raise ImportError("Unable to identify OPS version. "
                          + "OPS probably not installed correctly.")

    # continue force-setting version based on `setup.py`
    import imp  # may be Py2 only!
    ops_setup = imp.load_source("ops_setup", setupfile)
    version = imp.new_module("openpathsampling.version")

    version.version = ops_setup.preferences['version']
    version.short_version = ops_setup.preferences['version']
    version.git_version = ops_setup.get_git_version()
    version.full_version = ops_setup.preferences['version']
    if not ops_setup.preferences['released']:
        version.full_version += ".dev-" + version.git_version[:7]
    isrelease = str(ops_setup.preferences['released'])

from .analysis.path_histogram import PathDensityHistogram

from .analysis.replica_network import (
    ReplicaNetwork, trace_ensembles_for_replica,
    trace_replicas_for_ensemble, condense_repeats,
    ReplicaNetworkGraph
)

from .analysis.shooting_point_analysis import (
    ShootingPointAnalysis, SnapshotByCoordinateDict
)

from .analysis.reactive_flux_analysis import (
    ReactiveFluxAnalysis
)

from .analysis.trajectory_transition_analysis import (
    TrajectoryTransitionAnalysis,
    TrajectorySegmentContainer
)

from .analysis.channel_analysis import ChannelAnalysis

from .bias_function import (
    BiasFunction, BiasLookupFunction, BiasEnsembleTable,
    SRTISBiasFromNetwork
)

# from .collectivevariable import (
    # FunctionCV, MDTrajFunctionCV, MSMBFeaturizerCV,
    # InVolumeCV, CollectiveVariable, CoordinateGeneratorCV,
    # CoordinateFunctionCV, CallableCV, PyEMMAFeaturizerCV,
    # GeneratorCV)

from .ensemble import (
    Ensemble, EnsembleCombination,
    EmptyEnsemble, FullEnsemble, PartInXEnsemble,
    AllInXEnsemble, AllOutXEnsemble, WrappedEnsemble,
    SuffixTrajectoryEnsemble, PrefixTrajectoryEnsemble,
    PartOutXEnsemble, LengthEnsemble, NegatedEnsemble,
    ReversedTrajectoryEnsemble, SequentialEnsemble, VolumeEnsemble,
    SequentialEnsemble, IntersectionEnsemble, UnionEnsemble,
    SingleFrameEnsemble, MinusInterfaceEnsemble, TISEnsemble,
    OptionalEnsemble, join_ensembles
)

from .step_visualizer_2D import StepVisualizer2D

from .movechange import (
    EmptyMoveChange, ConditionalSequentialMoveChange,
    NonCanonicalConditionalSequentialMoveChange,
    MoveChange, PartialAcceptanceSequentialMoveChange,
    RandomChoiceMoveChange, SampleMoveChange,
    SequentialMoveChange, KeepLastSampleMoveChange,
    FilterSamplesMoveChange,
    PathSimulatorMoveChange, AcceptedSampleMoveChange,
    RejectedSampleMoveChange, SubMoveChange,
    FilterByEnsembleMoveChange, RejectedNaNSampleMoveChange,
    RejectedMaxLengthSampleMoveChange
)

from .pathmover import Details, MoveDetails, SampleDetails

from .pathmover import (
    RandomChoiceMover, PathMover, ConditionalSequentialMover,
    NonCanonicalConditionalSequentialMover,
    PartialAcceptanceSequentialMover, BackwardShootMover, ForwardShootMover,
    BackwardExtendMover, ForwardExtendMover, MinusMover,
    SingleReplicaMinusMover, PathReversalMover,
    ReplicaExchangeMover, EnsembleHopMover,
    SequentialMover, ConditionalMover,
    PathSimulatorMover, PathReversalSet,
    SampleMover, StateSwapMover, FinalSubtrajectorySelectMover, EngineMover,
    FirstSubtrajectorySelectMover,
    OneWayShootingMover, RandomSubtrajectorySelectMover, SubPathMover,
    EnsembleFilterMover, SelectionMover, FirstAllowedMover,
    LastAllowedMover, OneWayExtendMover, SubtrajectorySelectMover,
    IdentityPathMover, RandomAllowedChoiceMover,
    TwoWayShootingMover, ForwardFirstTwoWayShootingMover,
    BackwardFirstTwoWayShootingMover
)

from .pathsimulators import (
    PathSimulator, FullBootstrapping, Bootstrapping, PathSampling, MCStep,
    CommittorSimulation, ReactiveFluxSimulation, DirectSimulation,
    ShootFromSnapshotsSimulation
)

from .rng import default_rng
from .sample import Sample, SampleSet

from .shooting import (
    ShootingPointSelector, UniformSelector, GaussianBiasSelector,
    FirstFrameSelector, FinalFrameSelector, InterfaceConstrainedSelector,
    BiasedSelector
)

from .snapshot_modifier import (
    NoModification, RandomVelocities, VelocityDirectionModifier,
    SingleAtomVelocityDirectionModifier
)

from .storage.storage import Storage, AnalysisStorage

from .volume import (
    Volume, VolumeCombination,
    EmptyVolume, FullVolume, CVDefinedVolume, PeriodicCVDefinedVolume,
    IntersectionVolume, UnionVolume, SymmetricDifferenceVolume,
    RelativeComplementVolume, join_volumes
)

# from .high_level import move_strategy as strategies
from . import strategies

from .high_level.move_scheme import (
    MoveScheme, DefaultScheme, LockedMoveScheme, SRTISScheme,
    OneWayShootingMoveScheme
)

from .high_level.transition import (
    TISTransition, Transition, TPSTransition, FixedLengthTPSTransition
)

from .high_level.network import (
    MSTISNetwork, TransitionNetwork, MISTISNetwork, TPSNetwork,
    FixedLengthTPSNetwork
)

from .high_level.interface_set import (
    InterfaceSet, VolumeInterfaceSet, PeriodicVolumeInterfaceSet
)

from .high_level.ms_outer_interface import MSOuterTISInterface

from .high_level.part_in_b_tps import (
    PartInBFixedLengthTPSNetwork, PartInBFixedLengthTPSTransition
)

from .ensembles import *
from .pathmovers import *
from .collectivevariables import *
from .pathmovers.move_schemes import *

import openpathsampling.numerics as numerics
import openpathsampling.beta

from openpathsampling.engines import Trajectory, BaseSnapshot

# until engines are proper subpackages, built-ins need to be findable!
import openpathsampling.engines.openmm #as openmm
import openpathsampling.engines.toy #as toy


def git_HEAD():  # pragma: no cover
    from subprocess import check_output
    import os.path
    git_dir = os.path.dirname(os.path.realpath(__file__))
    return check_output(["git", "-C", git_dir, "rev-parse", "HEAD"])[:-1]
    # chops the newline at the end


import os.path

resources_directory = os.path.join(os.path.dirname(__file__), 'resources')
