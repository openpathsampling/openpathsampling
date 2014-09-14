"""
@author David W.H. Swenson
"""

import TrajFile
import trajectory
import storage

class TranslatorXYZ(object):
    '''
    We consider two data formats: xyz and storage (the latter being native
    for this code). Each format has associated with it trajectory object.
    For xyz, the object is TrajFile.TrajFile. For storage, the object is
    trajectory.Trajectory.

    A few relevant object names:
        self.trajfile is of type TrajFile.TrajFile 
        self.trajectory is of type trajectory.Trajectory

    Once one of these objects has been loaded, we can call
    trajfile2trajectory or trajectory2trajfile.
    '''
    def load_trajfile(self, tfile):
        '''Loads xyz file into self.traj, which is a TrajFile object'''
        self.traj = TrajFile.TrajFile.read_xyz(tfile)

    def load_from_storage(self, storage):
        '''Loads data from storage object into self.trajectory'''
        pass

    def trajfile2trajectory(self):
        pass

    def trajectory2trajfile(self):
        pass

    def regularize(self):
        '''Assuming we loaded one of the objects, make the other one'''
        if self.trajectory==None:
            self.trajfile2trajectory()
        if self.trajfile==None:
            self.trajectory2trajfile()


    def output_storage(self, outfname):
        pass

    def output_xyz(self, outfname):
        self.traj.write_xyz(outname)


