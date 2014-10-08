'''
Created on 08.07.2014

@author: Jan-Hendrik Prinz
'''

import numpy as np

from tis_interfaces import TISInterfaces
from transition_interface_sampling import TransitionInterfaceSampling
from Simulator import Simulator
from orderparameter import OP_RMSD_To_Lambda, OP_Multi_RMSD
from volume import LambdaVolume, VoronoiVolume
from ensemble import EnsembleFactory as ef
import time
from pathmover import ForwardShootMover, BackwardShootMover, PathMover, MixedMover
from shooting import UniformSelector
from ensemble import LengthEnsemble, InXEnsemble, OutXEnsemble
from trajectory import Trajectory
from pymbar import MBAR
from snapshot import Snapshot

if __name__ == '__main__':
    simulator = Simulator.Alanine_system('auto')
    PathMover.simulator = simulator

    print "Currently", simulator.storage.number_of_trajectories(), "simulations in the storage"
    print "Currently", simulator.storage.number_of_configurations(), "total frames in the storage"

    if simulator.storage.number_of_trajectories() == 0:        
        # load initial equilibrate snapshot given by ID #0
        snapshot = Snapshot.load(0, 0)
        
        # generate from this snapshot a trajectory with 50 steps
        traj = simulator.generate(snapshot, [LengthEnsemble(slice(0,50))])
        print len(traj)
        traj.save()
        
        # Save as Multi-Frame pdb  (only alanine, no water !)  
        traj.solute.md().save_pdb('data/mdtraj.pdb', True)    
        
    if True:        
        cc = Trajectory.load(1)[ 0 ]
        op = OP_RMSD_To_Lambda('lambda1', cc, 0.00, 1.00, atom_indices=simulator.solute_indices, use_storage=True)

        dd = simulator.storage.trajectory(1)[ 0:50 ]

        lV = LambdaVolume(op, 0.0, 0.06)
        lV2 = LambdaVolume(op, 0.0, 0.08)
        start_time = time.time()
        
#        print time.time() - start_time
    
        # if this uses the same orderparameter it is fast, since the values are cached!
        tis = ef.TISEnsemble(
                       LambdaVolume(op, 0.0, 0.041),
                       LambdaVolume(op, 0.0, 0.041),
                       LambdaVolume(op, 0.0, 0.045),
                       True
                       )

        enAB = ef.A2BEnsemble(
                       LambdaVolume(op, 0.0, 0.041),
                       LambdaVolume(op, 0.0, 0.041),
                       True
                       )
        
        tt = simulator.storage.trajectory(1)[4:18]

        print [ (op(d)) for d in dd ]

        stime = time.time()
        print tis.locate(dd, lazy=True, overlap=50)
        print time.time() - stime
        stime = time.time()
        print tis.locate(dd, lazy=False, overlap=50)
        print time.time() - stime

        stime = time.time()
        print enAB.locate(dd, lazy=True, overlap=50)
        print time.time() - stime
        stime = time.time()
        print enAB.locate(dd, lazy=False, overlap=50)
        print time.time() - stime

        print 'Lazy is about 20 times faster in this example and 10 times with shortcircuit active!!!'
        print 'If Shortcircuit is active general speedup of 30%!!!'

        # This is to cache the values for all snapshots in tt. Makes later access MUCH faster. 
        # Especially because the frames do not have to be read one by one.
        op(tt)
        # print tis
#        print [ s.idx for s in tt]
#        print [ (lV(d)) for d in tt ]
        print [ (op(d)) for d in dd ]
        
        # This tests, if the iteration request works. It basically return True if it makes sense to simulate or if the ensemble cannot
        # be true in the next step. This should be passed to the pathmover to stop simulating for a particular ensemble

        vn = VoronoiVolume(
                OP_Multi_RMSD('Voronoi', tt[[0,10]], atom_indices=simulator.solute_indices),
                state = 0
                )

        print "Iteration test"
        for l in range(0,tt.frames + 0):
            print tis.forward(tt[0:l]), tis(tt[0:l]), lV(tt[l]), lV2(tt[l]), vn(tt[l]), vn.cell(tt[l])

        print "Iteration test"
        for l in range(0,tt.frames + 0):
            print tis.forward(tt[0:l]), tis(tt[0:l]), lV(tt[l]), lV2(tt[l]), vn(tt[l]), vn.cell(tt[l])
            
        print op(tt[0])
        s = Snapshot(coordinates = tt[0].coordinates)
        print op(s)
        s.idx = 0
        print op(s)
        
                           
        
        
                
    
        # Check if the trajectory goes from lambda < 0.06 to lambda >0.08 and back    
        print 'In ensemble'
        print tis(tt)
        
        en = ef.A2BEnsemble(lV, lV, True)
        print en(tt)

        en = InXEnsemble(lV, 0)
        print en(tt)

        en = InXEnsemble(lV, -1)
        print en(tt)

        en = OutXEnsemble(lV, slice(1,-1), lazy = False)
        print en(tt)
        
#        Simulator.op = op
        
        bm = BackwardShootMover(
                selector = UniformSelector(),
                ensemble = tis
                )

        fm = ForwardShootMover(
                selector = UniformSelector(),
                ensemble = tis
                )

        pm = MixedMover([bm, fm])
        
        pm.move(tt)
        print 'ensemble Check :', pm.ensemble(tt)
        
        print pm.ensemble(pm.final)
        print 'Accepted : ', pm.accepted
        print len(pm.final)
        
        print 'Next Check:'
        
        
        en = ef.A2BEnsemble(lV, lV, True)
        print en(pm.final)

        en = InXEnsemble(lV, 0)
        print en(pm.final)

        en = InXEnsemble(lV, -1)
        print en(pm.final)

        en = OutXEnsemble(lV, slice(1,-1), lazy = False)
        print en(pm.final)

        op.save()

        exit()

           
#    print "Clustering"

#    msm = VoronoiTesselation()
#    msm.storage = Simulator.storage
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
    
    # TODO: Change TISInterfaces to use only Trajectory object and not mdtraj.Trajectory
    tis = TISInterfaces(simulator.storage, simulator.solute_indices, traj.md(), np.array([[0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080], [0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080]]))
    tis.create_rmsd_metric()

    tis.assign_storage()
    
    initial = tis.split_into_connections(simulator.storage.last_trajectory())
    
    sampling = TransitionInterfaceSampling(simulator, tis)
    old = initial[0]
    
    for i in range(10):
        new = sampling.sampleTrajectory(old)
        new.save()
             
        old = new    

#    mat = tis.compute_N_IJ_l_kI()
#    print mat

    
    if simulator.storage.number_of_trajectories() > 100:  
        u_kln, N_k = tis.compute_mbar_array(0, 0, infinite_energy = 1e5)
            
        mbar = MBAR(u_kln, N_k, method = 'adaptive', relative_tolerance=1.0e-10, verbose=False)
        
        weights =  mbar.getWeights()
                    
        n = np.zeros(weights.shape, dtype='float')
        
        print np.exp(-mbar.f_k)    
        print N_k
        
        index = 0
        
        for k,n_samples in enumerate(N_k):
            n[index:index+int(n_samples),0:k+1] = 1.0
            index += int(n_samples)
            
        print n
        print weights[:,0]
        print weights[:,2]
    
        for i in range(5):
            print np.sum(n[:,i] * weights[:,0])
        
        print np.sum(n * weights, axis=0)
                    
        print np.sum(weights,axis=0)
                       
#    Simulator.storage.ncfile.close()
#    exit()
 
#    print tis.cores
#    print tis.interfaces
#    print tis.distances
        
#    print msm.assign_index_trajectory(traj.indices())
#    print msm.assign_index_trajectory(traj2.indices())
#    print msm.assign_index_trajectory(traj3.indices())

#    print msm.assign(traj, atom_indices)
#    print msm.assign_first(traj[0], atom_indices)

#    print tis.assign_snapshot(traj[10])
    
#    traj_test = Simulator.generate(traj[15], 200, stopping = tis.within_interface(4) )
#    traj_test.save()
#    tis.assign_storage()

#    print traj_test.indices()
#    print len(traj_test)
    
#    tis.assign_all_trajectories()

#    print tis.assign_all_trajectories()

    for l in range(1,simulator.storage.number_of_trajectories()+1):
        traj = simulator.storage.trajectory(l)

    print "Last"
    ltraj = simulator.storage.last_trajectory()

    initial = tis.split_into_connections(simulator.storage.last_trajectory())
    
    sampling = TransitionInterfaceSampling(simulator, tis)
    old = initial[0]
    
    for i in range(10):
        new = sampling.sampleTrajectory(old)
        new.save()
             
        old = new
    
#    tis.assign_storage()    
#    mat = tis.compute_N_IJ_l_kI()
#    print mat

    simulator.storage.ncfile.close()
