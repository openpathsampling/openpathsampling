'''
Created on 08.07.2014

@author: jan-hendrikprinz
'''

import numpy as np

from tis_interfaces import TISInterfaces
from transition_interface_sampling import TransitionInterfaceSampling
from simulator import Simulator

if __name__ == '__main__':
    simulator = Simulator.Alanine_system('auto')

    print "Currently ", simulator.storage.number_of_trajectories(), " simulations in the storage"

    if simulator.storage.number_of_trajectories() == 0:        
        # load initial equilibrate snapshot given by ID #0
        snapshot = simulator.storage.snapshot(0)    
        
        # generate from this snapshot a trajectory with 50 steps
        traj = simulator.generate(snapshot, 50)
        traj.save()
        
        # Save as Multi-Frame pdb  (only alanine, no water !)  
        traj.solute.md().save_pdb('data/mdtraj.pdb', True)    
        
#    print "Clustering"

#    msm = VoronoiTesselation()
#    msm.storage = simulator.storage
#    msm.n_centers = 10
#    msm.atom_indices = atom_indices

    # Use RMSD metric with 5 cluster centers
#    msm.create_rmsd_metric()
#    msm.cluster()
#    msm.assign_storage()
    
#    print msm.assign_all_trajectories()


    
    print "Interfaces"    
    # Create Interaces around first and last frame of the first trajectory and use only solute coordinates
    
    traj = simulator.storage.trajectory(1)[ [0,-1] ].solute
    
    tis = TISInterfaces(simulator.storage, simulator.solute_indices, traj.md(), np.array([[0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080], [0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080]]))
    tis.create_rmsd_metric()

    tis.assign_storage()
    
#    print tis.cores
#    print tis.interfaces
#    print tis.distances
        
#    print msm.assign_index_trajectory(traj.indices())
#    print msm.assign_index_trajectory(traj2.indices())
#    print msm.assign_index_trajectory(traj3.indices())

#    print msm.assign(traj, atom_indices)
#    print msm.assign_first(traj[0], atom_indices)

#    print tis.assign_snapshot(traj[10])
    
#    traj_test = simulator.generate(traj[15], 200, stopping = tis.within_interface(4) )
#    traj_test.save()
#    tis.assign_storage()

#    print traj_test.indices()
#    print len(traj_test)
    
#    tis.assign_all_trajectories()

#    print tis.assign_all_trajectories()

    initial = tis.split_into_connections(simulator.storage.last_trajectory())
    
    sampling = TransitionInterfaceSampling(simulator, tis)
    
    old = initial[0]
    
    for i in range(5):
        new = sampling.sampleTrajectory(old)
        new.save()
             
        old = new
        
    
#    tis.assign_storage()    
#    mat = tis.compute_N_IJ_l_kI()
#    print mat

    simulator.storage.ncfile.close()