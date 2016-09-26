# coding: utf-8
import matplotlib.pyplot as plt
import openpathsampling as paths
import numpy as np


# =============================================================================
# ANALYZING THE MSTIS SIMULATION
# =============================================================================
print """ANALYZING THE MSTIS SIMULATION"""

storage = paths.AnalysisStorage("ala_mstis_production.nc")
print "PathMovers:", len(storage.pathmovers)
print "Engines:", len(storage.engines)
print "Samples:", len(storage.samples)
print "Ensembles:", len(storage.ensembles)
print "SampleSets:", len(storage.samplesets)
print "Snapshots:", len(storage.snapshots)
print "Trajectories:", len(storage.trajectories)
print "Networks:", len(storage.networks)

mstis = storage.networks[0]

# Reaction rates
mstis.hist_args['max_lambda'] = { 'bin_width' : 2, 'bin_range' : (0.0, 90) }
mstis.hist_args['pathlength'] = { 'bin_width' : 5, 'bin_range' : (0, 100) }

print 'THE RATE MATRIX'
print mstis.rate_matrix(storage.steps, force=True)

# -----------------------------------------------------------------------------
# Total crossing probability
# -----------------------------------------------------------------------------
print """Total crossing probability"""

stateA = storage.volumes["A0"]
stateB = storage.volumes["B0"]
stateC = storage.volumes["C0"]
tcp_AB = mstis.transitions[(stateA, stateB)].tcp
tcp_AC = mstis.transitions[(stateA, stateC)].tcp
tcp_BC = mstis.transitions[(stateB, stateC)].tcp
tcp_BA = mstis.transitions[(stateB, stateA)].tcp
tcp_CA = mstis.transitions[(stateC, stateA)].tcp
tcp_CB = mstis.transitions[(stateC, stateB)].tcp

plt.plot(tcp_AB.x, tcp_AB)
plt.plot(tcp_CA.x, tcp_CA)
plt.plot(tcp_BC.x, tcp_BC)
plt.plot(tcp_AC.x, tcp_AC) # same as tcp_AB in MSTIS
plt.plot(tcp_AB.x, np.log(tcp_AB))
plt.plot(tcp_CA.x, np.log(tcp_CA))
plt.plot(tcp_BC.x, np.log(tcp_BC))


# -----------------------------------------------------------------------------
# Flux
# -----------------------------------------------------------------------------
print """Flux"""

import pandas as pd
flux_matrix = pd.DataFrame(columns=mstis.states, index=mstis.states)
for state_pair in mstis.transitions:
    transition = mstis.transitions[state_pair]
    flux_matrix.set_value(state_pair[0], state_pair[1], transition._flux)

flux_matrix


# -----------------------------------------------------------------------------
# Conditional transition probability
# -----------------------------------------------------------------------------
print """Conditional transition probability"""

outer_ctp_matrix = pd.DataFrame(columns=mstis.states, index=mstis.states)
for state_pair in mstis.transitions:
    transition = mstis.transitions[state_pair]
    outer_ctp_matrix.set_value(state_pair[0], state_pair[1], transition.ctp[transition.ensembles[-1]])    

outer_ctp_matrix
ctp_by_interface = pd.DataFrame(index=mstis.transitions)
for state_pair in mstis.transitions:
    transition = mstis.transitions[state_pair]
    for ensemble_i in range(len(transition.ensembles)):
        ctp_by_interface.set_value(
            state_pair, ensemble_i,
            transition.conditional_transition_probability(
                storage.steps,
                transition.ensembles[ensemble_i]
        ))
    
    
ctp_by_interface

# Path ensemble properties
hists_A = mstis.transitions[(stateA, stateB)].histograms
hists_B = mstis.transitions[(stateB, stateC)].histograms
hists_C = mstis.transitions[(stateC, stateB)].histograms


# -----------------------------------------------------------------------------
# Interface crossing probabilities
# -----------------------------------------------------------------------------
print """Interface crossing probabilities"""

for hist in [hists_A, hists_B, hists_C]:
    for ens in hist['max_lambda']:
        normalized = hist['max_lambda'][ens].normalized()
        plt.plot(normalized.x, normalized)
# add visualization of the sum
for hist in [hists_A, hists_B, hists_C]:
    for ens in hist['max_lambda']:
        reverse_cumulative = hist['max_lambda'][ens].reverse_cumulative()
        plt.plot(reverse_cumulative.x, reverse_cumulative)
for hist in [hists_A, hists_B, hists_C]:
    for ens in hist['max_lambda']:
        reverse_cumulative = hist['max_lambda'][ens].reverse_cumulative()
        plt.plot(reverse_cumulative.x, np.log(reverse_cumulative))


# -----------------------------------------------------------------------------
# Path length histograms
# -----------------------------------------------------------------------------
print """Path length histograms"""

for hist in [hists_A, hists_B, hists_C]:
    for ens in hist['pathlength']:
        normalized = hist['pathlength'][ens].normalized()
        plt.plot(normalized.x, normalized)
for ens in hists_A['pathlength']:
    normalized = hists_A['pathlength'][ens].normalized()
    plt.plot(normalized.x, normalized)

# Sampling properties



# -----------------------------------------------------------------------------
# Move scheme analysis
# -----------------------------------------------------------------------------
print """Move scheme analysis"""

scheme = storage.schemes[0]
scheme.move_summary(storage.steps)
scheme.move_summary(storage.steps, 'shooting')
scheme.move_summary(storage.steps, 'minus')
scheme.move_summary(storage.steps, 'repex')
scheme.move_summary(storage.steps, 'pathreversal')


# -----------------------------------------------------------------------------
# Replica exchange sampling
# -----------------------------------------------------------------------------
print """Replica exchange sampling"""

repx_net = paths.ReplicaNetwork(scheme, storage.steps)

print """# Replica exchange mixing matrix"""
repx_net.mixing_matrix()

print """# Replica exchange graph"""
repxG = paths.ReplicaNetworkGraph(repx_net)
repxG.draw('spring')

print """# Replica exchange flow"""


# -----------------------------------------------------------------------------
# Replica move history tree
# -----------------------------------------------------------------------------
print """Replica move history tree"""

import openpathsampling.visualize as vis
reload(vis)
from IPython.display import SVG
tree = vis.PathTree(
    [step for step in storage.steps if not isinstance(step.change, paths.EmptyPathMoveChange)],
    vis.ReplicaEvolution(replica=3, accepted=False)
)
tree.options.css['width'] = 'inherit'

SVG(tree.svg())
decorrelated = tree.generator.decorrelated
print "We have " + str(len(decorrelated)) + " decorrelated trajectories."



# -----------------------------------------------------------------------------
# Visualizing trajectories
# -----------------------------------------------------------------------------
print """Visualizing trajectories"""


# Histogramming data (TODO)

