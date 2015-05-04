'''
Created on 08.07.2014

@author: Jan-Hendrik Prinz
'''

import numpy as np

from openpathsampling.Simulator import Simulator
from openpathsampling.collectivevariable import CV_RMSD_To_Lambda, CV_Multi_RMSD
from openpathsampling.volume import LambdaVolume, VoronoiVolume
from openpathsampling.ensemble import EnsembleFactory as ef
import time
from openpathsampling.pathmover import ForwardShootMover, BackwardShootMover, PathMover, MixedMover
from openpathsampling.shooting import UniformSelector
from openpathsampling.ensemble import LengthEnsemble, AllInXEnsemble, AllOutXEnsemble
from openpathsampling.trajectory import Trajectory
from pymbar import MBAR
from openpathsampling.snapshot import Snapshot
from openpathsampling.openmm_simulation import OpenMMSimulation

if __name__ == '__main__':
    simulator = Simulator.Alanine_system('auto')

    PathMover.simulator = simulator
    storage = simulator.storage

    print "Currently", simulator.storage.trajectories.count(), "simulations in the storage"
    print "Currently", simulator.storage.configurations.count(), "total frames in the storage"

    Trajectory.storage = simulator.storage.trajectory

    if simulator.storage.trajectories.count() == 0:
        # load initial equilibrate snapshot given by ID #0
        snapshot = simulator.storage.snapshots.load(0)

        # generate from this snapshot a trajectory with 50 steps
        traj = simulator.generate(snapshot, [LengthEnsemble(slice(0,6))])
        simulator.storage.trajectories.save(traj)

        print len(traj)

        # Save as Multi-Frame pdb  (only alanine, no water ! and overwrite existing)
        traj.solute.md().save_pdb('data/mdtraj.pdb', True)
        
    if True:
        cc = Trajectory.storage.load(1)[ 0 ]
        cc = storage.snapshots.load(0)
        op = CV_RMSD_To_Lambda('lambda1', cc, 0.00, 1.00, atom_indices=simulator.solute_indices)
        storage.collectivevariables.restore(op)
        print op(Trajectory.storage.load(1)[0:2])
        dd = simulator.storage.trajectories.load(1)[ 0:6 ]
        lV = LambdaVolume(op, 0.0, 0.06)
        lV2 = LambdaVolume(op, 0.0, 0.08)

        # if this uses the same collectivevariable it is fast, since the values are cached!
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

        tt = simulator.storage.trajectories.load(1)

#        print [ (op(d)) for d in dd ]
#        op.save()

        stime = time.time()
        print tis.find_valid_slices(dd, lazy=True, overlap=1)
        print time.time() - stime
        stime = time.time()
        print tis.find_valid_slices(dd, lazy=False, overlap=1)
        print time.time() - stime

        stime = time.time()
        print enAB.find_valid_slices(dd, lazy=True, overlap=1)
        print time.time() - stime
        stime = time.time()
        print enAB.find_valid_slices(dd, lazy=False, overlap=1)
        print time.time() - stime

        # This is to cache the values for all snapshots in tt. Makes later access MUCH faster. 
        # Especially because the frames do not have to be read one by one.
        op(tt)
        # print tis
#        print [ s.idx for s in tt]
#        print [ (lV(d)) for d in tt ]
        print "Results :"
        print [ op(d) for d in dd ]
        
        # This tests, if the iteration request works. It basically return True if it makes sense to simulate or if the ensemble cannot
        # be true in the next step. This should be passed to the pathmover to stop simulating for a particular ensemble

        vn = VoronoiVolume(
                CV_Multi_RMSD('Voronoi', tt[[0,2]], atom_indices=simulator.solute_indices),
                state = 0
                )

        print "Iteration test"
        for l in range(0,len(tt) + 0):
            print tis.can_append(tt[0:l]), tis(tt[0:l]), lV(tt[l]), lV2(tt[l]), vn(tt[l]), vn.cell(tt[l])

        print "Iteration test"
        for l in range(0,len(tt) + 0):
            print tis.can_append(tt[0:l]), tis(tt[0:l]), lV(tt[l]), lV2(tt[l]), vn(tt[l]), vn.cell(tt[l])

        print op(tt[0])
        s = Snapshot(coordinates=tt[0].coordinates)
        print op(s)

        # Check if the trajectory goes from lambda < 0.06 to lambda >0.08 and back    
        print 'In ensemble'
        print tis(tt)
        
        en = ef.A2BEnsemble(lV, lV, True)
        print en(tt)

        en = AllInXEnsemble(lV, 0)
        print en(tt)

        en = AllInXEnsemble(lV, -1)
        print en(tt)

        en = AllOutXEnsemble(lV, slice(1,-1), lazy = False)
        print en(tt)

        storage.ensembles.save(en)

        bm = BackwardShootMover(
                selector = UniformSelector(),
                ensemble = tis
                )

        fm = ForwardShootMover(
                selector = UniformSelector(),
                ensemble = tis
                )

        # Not using MixedMover because its implementation is inconsistent with being an actual pathmover at the moment
        rand = np.random.random()
        if rand < 0.5:
            mm = fm
        else:
            mm = bm

        mm = MixedMover([bm, fm])

        tt = storage.trajectories.last()

        pth = mm.move(tt)

        storage.results.save(pth)

        print pth.details.json

        loaded = storage.pathmovers.load(mm.idx[storage])

        print 'ensemble Check :', mm.ensemble(pth.trajectory)

        print mm.ensemble(pth.details.final)
        print 'Accepted : ', pth.details.accepted
        print len(pth.details.final)
        
        print 'Next Check:'

        en = ef.A2BEnsemble(lV, lV, True)
        print en(pth.details.final)

        en = AllInXEnsemble(lV, 0)
        print en(pth.details.final)

        en = AllInXEnsemble(lV, -1)
        print en(pth.details.final)

        en = AllOutXEnsemble(lV, slice(1, -1), lazy = False)
        print en(pth.details.final)

        op.save(storage=storage.collectivevariables)

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

#    tis.assign_storage()    
#    mat = tis.compute_N_IJ_l_kI()
#    print mat

    simulator.storage.ncfile.close()
