"""
@author David W.H. Swenson
"""

import sys
import os

# re only for XYZTranslate.guess_fname_format():
import re
import numpy as np

from simtk.openmm.app.topology import Topology
from simtk.openmm.app.element import Element
from simtk.unit import amu
import mdtraj

# I assume the next directory above us is where the msm-tis classes hide
sys.path.append(os.path.abspath('../'))
import TrajFile
import trajectory
import snapshot
from storage import TrajectoryStorage

class XYZTranslator(object):
    '''
    We consider two data formats: xyz and storage (the latter being native
    for this code). Each format has associated with it trajectory object.
    For xyz, the object is TrajFile.TrajFile. For storage, the object is
    trajectory.Trajectory.

    A few relevant object names:
        self.trajfile is of type TrajFile.TrajFile 
        self.trajectory is of type trajectory.Trajectory

    Once one of these objects has been loaded, we can call self.regularize()
    to generate the other. Then we have an object which can output in either
    format.

    NB: there is no support yet for converting periodic boundary conditions
    between the to formats. In principle, that shouldn't be too hard, but it
    isn't necessary for current purposes (playing with 2D models with no
    periodicity). 
    '''
    
    def __init__(self):
        self.traj = None
        self.trajectory = None
        self.storage = None

    def guess_fname_format(self, fname):
        '''Takes an file name and tries to guess the format from that (under
        the assumption that I'm using my normal convention that the last
        number before the .xyz is the trajectory number).'''
        num = re.search(".*[^0-9]([0-9]+)[^0-9]*.xyz", fname).group(1)
        numfmt = "%0"+str(len(num))+"d"
        fmt = re.sub(num, numfmt, fname)
        return fmt

    def load_trajfile(self, tfile):
        '''Loads xyz file into self.traj, which is a TrajFile object'''
        self.traj = TrajFile.TrajFile().read_xyz(tfile)

    def load_from_storage(self, storage, loadnum=0):
        '''Loads data from storage object into self.trajectory'''
        self.trajectory = trajectory.Trajectory.load(loadnum)

    def trajfile_topology(self, trajfile):
        '''Creates a (fake-ish) OpenMM Topology from an xyz file'''
        topol = Topology()
        # assume that the atoms are constant through the xyz trajectory, so
        # we can just use the first frame:
        myframe = trajfile.frames[0]
        chain = topol.addChain()
        for atom in range(myframe.natoms):
            # make each atom a separate element and a separate residue
            label = "_"+myframe.labels[atom] # underscore to avoid conflicts
            mass = myframe.mass[atom]
            element = Element(   
                                number=atom+1, # abnormal
                                name=label, symbol=label, mass=mass*amu
                             )
            # we need to add the elements to the mdtraj dictionaries, too
            mdtrajelem = mdtraj.element.Element( 
                                number=atom+1,
                                name=label, symbol=label, mass=mass
                            )
            #element = mdtraj.element.Element(atom,label,label,mass)
            residue = topol.addResidue(label, chain)
            topol.addAtom(label, element, residue)
        return topol # note that we can easily make this into mdtraj

    def init_storage(self, fname="test.nc"):
        topol = self.trajfile_topology(self.trajfile)
        self.storage = TrajectoryStorage( topology=topol,
                                          filename=fname, 
                                          mode='auto')
    

    def trajfile2trajectory(self, trajfile):
        '''Converts TrajFile.TrajFile to trajectory.Trajectory'''
        res = trajectory.Trajectory()
        # make sure there's some storage in place
        if (res.storage == None):
            if (self.storage == None):
                
                self.init_storage()
            res.storage = self.storage

        for frame in trajfile.frames:
            mysnap = snapshot.Snapshot( coordinates=np.array(frame.pos),
                                        velocities=np.array(frame.vel) )
            mysnap.save()
            res.append(mysnap)
            res.save()
        return res

    def trajectory2trajfile(self, trajectory):
        '''Converts trajectory.Trajectory to TrajFile.TrajFile'''
        res = TrajFile.TrajFile()
        for snap_i in len(range(trajectory.snapshots())):
            mysnap = trajectory.snapshot(snap_i)
            myframe = TrajFrame()
            myframe.natoms = mysnap.atoms()
            myframe.vel = mysnap.velocities()
            myframe.pos = mysnap.positions()
            #myframe.label = mysnap.labels() # TODO: atomlabels in snap
            myframe.mass = mysnap.masses()
            res.frames.append(myframe)
        return res

    def regularize(self):
        '''Assuming we loaded one of the objects, make the other one'''
        if self.trajectory==None:
            self.trajectory = self.trajfile2trajectory(self.traj)
        if self.trajfile==None:
            self.traj = self.trajectory2trajfile(self.trajectory)


    def output_storage(self, outfname=sys.stdout):
        '''Writes trajectory to `outfname` as NetCDF'''
        pass

    def output_xyz(self, outfname=sys.stdout):
        '''Writes trajectory to `outfname` as .xyz file'''
        self.traj.write_xyz(outfname)


