## OpenMM based SR-MS-TIS

### [ Python classes ]
1. Trajectory - a list of snapshots
  member variables
  member functions
  inherited from list
  length (number of snapshots)
  number of atoms
  get list of [snapshot property]
  reverse - reverse direction of trajectory without regard for velocities

2. Snapshot - a trajectory frame in a context, 
  member variables are
  time ?,
  coordinates, 
  velocities ?, 
  simulation box size, 
  potential energy
  kinetic energy 
  member functions to compute 
  total energy, ...
  box volume

3. MSTISSimulation

### Questions ?
- keep number of atoms fixed for a trajectory or allow grand canonical
- keep volume (box size) fixed or allow for NpT
- keep energy fixed or allow for NpT / NVT
  => set a fixed ensemble type for the trajectory
- functions to compute ensemble properties from list of snapshots?
- how much needs to be deep copied (memory efficiency)
- What about Snapshots. Stored with or withour solvant (only necessary for explicit water). Is it enough to always solvate the protein again for next shooting, etc.
- For shooting moves all water is necessary, I guess only to keep chances of the new path to connect states A and B high enough.

### What else do we need?
  State decomposition -Way to define Core sets arbitrarily. Best from MSM state decomposition.
  The trajectory object, do they need to 

### Further considerations
- Every trajectory in TPS or TIS is "valid" there is no bias involved. Trajectories are only rejected if they do not connect specific start and end states, but they are still part of all possible pathways of a specific length. So we could use their data
- TIS is effective because it stops simulations once a certain criterion has been met
- One problem is to identify metastable states as we go along. MSMs can help with this

### State Definition Necessary Properties:
- Easy to measure
- Needs to be extendable / growable by a continuous parameter ?!?
- Both, Voronoi MSM and TIS, RMSD based state definition use cluster centers, which is potentially the most effective way to define regions in configurational space 
- If we combine both ideas then we use the unweighted voronoi idea to decompose close states, while inside voronoi cells a metric-based approach is used. This could also be with a non-uniform metric definition.

### If we use only a MSM then for
- 1 state. Compute Commute Distance for all states and use this to define states of a specific level of stability.
- 1 state. Compute Mean First Passage Time from Center to state x and use this to define states of a specific time distance to identify states.
- 2+ states. Compute Multi-state Committor and use a percentile based way (no timescales are included anymore!)
