"""
@author David W.H. Swenson
"""

import sys
import os

import re
import optparse

import numpy as np

from simtk.openmm.app.topology import Topology
from simtk.openmm.app.element import Element
from simtk.unit import amu, nanometers, picoseconds, Quantity
import mdtraj

# I assume the next directory above us is where the msm-tis classes hide
sys.path.append(os.path.abspath('../'))
import TrajFile
import trajectory
import snapshot
from storage import TrajectoryStorage

class AtomCounter(object):
    '''Let's be honest: that's all we're using the simulation.system object
    for. So I'll duck-punch.'''
    def __init__(self, natoms):
        self.natoms = natoms

    def getNumParticles(self):
        '''QUAAAAACK'''
        return self.natoms

class SimulationDuckPunch(object):
    '''This is what happens when you find a stranger in the Alps.'''
    def __init__(self, topology, system):
        self.system = system
        self.topology = topology
        

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
        self.n_frames_max = 10000
        self.outfile = None
        self.intype = None
        self.outtype = None

    def guess_fname_format(self, fname):
        '''Takes an file name and tries to guess the format from that (under
        the assumption that I'm using my normal convention that the last
        number before the .xyz is the trajectory number).'''
        num = re.search(".*[^0-9]([0-9]+)[^0-9]*.xyz", fname).group(1)
        numfmt = "%0"+str(len(num))+"d"
        fmt = re.sub(num, numfmt, fname)
        return fmt

    # TODO: this should be a class method, right?
    def build_parser(self, parser):
        parser.add_option("-o", "--output", help="output file")
        parser.add_option("-t", "--topology", help="topology (PDB or XYZ)")
        return parser

    def set_infile_outfile(self, infiles, outfile):
        if re.search("\/", outfile):
            directory=re.match("(.*)\/[^\/]*",outfile).group(1)
            if not os.path.isdir(directory):
                os.makedirs(directory)
        self.infiles = []
        for f in infiles:
            if re.match(".*\.xyz$", f):
                self.intype = "xyz"
                self.outtype = "nc"
                self.infiles.append(f)
                self.outfile = (outfile, "nc")
            elif re.match(".*\.nc$", f):
                self.intype = "nc"
                self.outtype = "xyz"
                self.infiles.append(f)
        if self.intype == "xyz" and not re.match(".*\.nc$", outfile): 
            outfile = outfile+".nc"
        elif self.intype == "nc" and not re.match(".*\.xyz$", outfile): 
            outfile = outfile+".xyz"
        self.outfile = outfile
        return

    def topology_from_file(self, fname):
        if re.match(".*\.xyz$", fname):
            topo_f = TrajFile.TrajFile()
            topo_f.read_xyz(fname)
            TrajFile.set_default_mass(topo_f)
            return self.trajfile_topology(topo_f)
        else:
            print "Topology must be .xyz for now"
        return None
        
    def load_trajfile(self, tfile):
        '''Loads xyz file into self.trajfile, which is a TrajFile object'''
        self.trajfile = TrajFile.TrajFile().read_xyz(tfile)

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
            # Check whether atom is already in the topol; add it if not.
            # I hate this approach, since it assume the internal structure
            # of the Element classes
            try:
                element = Element.getBySymbol(label)
            except KeyError:
                element = Element(   
                                    number=atom+1, # abnormal
                                    name=label, symbol=label, mass=mass*amu
                                 )
            # we need to add the elements to the mdtraj dictionaries, too
            try:
                mdtrajelem = mdtraj.element.Element.getBySymbol(label)
            except KeyError:
                mdtrajelem = mdtraj.element.Element( 
                                number=atom+1,
                                name=label, symbol=label, mass=mass
                                )
            residue = topol.addResidue(label, chain)
            topol.addAtom(label, element, residue)
        return topol # note that we can easily make this into mdtraj

    def init_storage(self, fname):
        topol = self.trajfile_topology(self.trajfile)
        system = AtomCounter(self.trajfile.frames[0].natoms)
        self.simulation = SimulationDuckPunch(topol, system)
        self.storage = TrajectoryStorage( topology=topol,
                                          filename=fname, 
                                          mode='auto')
        snapshot.Snapshot.simulator = self
        self.storage.simulator = self
        self.storage.init_classes()
    

    def trajfile2trajectory(self, trajfile):
        '''Converts TrajFile.TrajFile to trajectory.Trajectory'''
        # make sure there's some storage in place
        res = trajectory.Trajectory()
        trajectory.Trajectory.simulator = self
        if (res.storage == None):
            if (self.storage == None):
                self.init_storage(self.outfile)
            res.storage = self.storage

        for frame in trajfile.frames:
            pos = Quantity(np.array(frame.pos),nanometers)
            vel = Quantity(np.array(frame.vel),nanometers/picoseconds)
            mysnap = snapshot.Snapshot( coordinates=pos,
                                        velocities=vel )
            mysnap.save()
            res.append(mysnap)
        res.save()
        return res

    def trajectory2trajfile(self, trajectory):
        '''Converts trajectory.Trajectory to TrajFile.TrajFile'''
        res = TrajFile.TrajFile()
        for snap in trajectory:
            myframe = TrajFile.TrajFrame()
            myframe.natoms = snap.atoms
            vel = snap.velocities / (nanometers/picoseconds)
            pos = snap.coordinates / nanometers
            myframe.vel = vel.tolist()
            myframe.pos = pos.tolist()
            topol = trajectory.simulator.simulation.topology 
            labels = []
            mass = []
            for atom in topol.atoms():
                labels.append(re.sub("^_", "", atom.element.name))
                mass.append(atom.element.mass / amu)
            myframe.mass = mass
            myframe.labels = labels
            res.frames.append(myframe)
        return res

    def translate(self):
        if (self.intype == "xyz"):
            for f in self.infiles:
                print "xyz->nc: "+f
                self.trajfile = TrajFile.TrajFile()
                self.trajfile.read_xyz(f)
                TrajFile.set_default_mass(self.trajfile)
                self.trajectory = self.trajfile2trajectory(self.trajfile)
        elif (self.intype == "nc"):
            self.storage = TrajectoryStorage( 
                                    topology=None,
                                    filename=self.infiles[0],
                                    mode='restore'
                                )
            self.storage.simulator = self
            #self.storage.verbose_root = True # DEBUG
            self.storage._restore_options(self) # not used? (for topol)
            # for now, topology must come from a separate file
            self.storage.topology = self.topology_from_file(self.topol_file)
            ntrajs = self.storage.number_of_trajectories()
            for traj_i in range(ntrajs):
                self.trajectory = self.storage.trajectory(traj_i+1)
                self.trajectory.simulator = self
                natoms = len(self.trajectory[0].coordinates) # ugly
                system = AtomCounter(natoms)
                self.simulation = SimulationDuckPunch(self.storage.topology,
                                                      system)
                self.trajfile = self.trajectory2trajfile(self.trajectory)
                myout = self.outfile
                if re.search("%[0-9]*d", self.outfile):
                    myout = self.outfile % (traj_i+1)
                print "nc->xyz ("+str(traj_i+1)+"/"+str(ntrajs)+"): "+myout
                self.trajfile.write_xyz(myout)

    def regularize(self):
        '''Assuming we loaded one of the objects, make the other one'''
        if self.trajectory==None:
            self.trajectory = self.trajfile2trajectory(self.traj)
        if self.trajfile==None:
            self.traj = self.trajectory2trajfile(self.trajectory)

    def output_xyz(self, outfname=sys.stdout):
        '''Writes trajectory to `outfname` as .xyz file'''
        self.traj.write_xyz(outfname)


if __name__=="__main__":
    parser = optparse.OptionParser()
    translator = XYZTranslator()
    parser = translator.build_parser(parser)
    (opts, args) = parser.parse_args(sys.argv[1:])
    translator.set_infile_outfile(args, opts.output)
    translator.topol_file = opts.topology
    translator.translate()
    

