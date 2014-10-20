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
from storage import Storage

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
    We consider two data formats: xyz and netCDF (the latter being native
    for this code). Each format has associated with it trajectory object.
    For xyz, the object is TrajFile.TrajFile. For netCDF, the object is
    trajectory.Trajectory. This class provides scripts to translate ***

    NB: there is no support yet for converting periodic boundary conditions
    between the to formats. In principle, that shouldn't be too hard, but it
    isn't necessary for current purposes (playing with 2D models with no
    periodicity). 

    Attributes
    ----------
    trajfile : TrajFile.TrajFile
        TrajFile version of our trajectory data
    trajectory : trajectory.Trajectory
        Trajectory version of our trajectory data
    infiles : list of str
        Input files (single .nc, single or multiple .xyz)
    outfile : str
        Output file as format string (i.e., can use notation like %06d to
        label the trajectory number)
    intype : str
        Input format type; either "xyz" or "nc"
    outtype : str
        Outfile format type; either "xyz" or "nc"
    topol_file : str
        Topology file. Required for nc->xyz conversion.
    storage : storage.TrajectoryStorage
        Object that handles writing netCDF format
    n_frames_max : int
        maximum frames per trajectory; required by current version of
        TrajectoryStorage to get XYZTranslator to quack like a Simulator.
    '''
    
    def __init__(self):
        self.trajfile = None
        self.trajectory = None
        self.outfile = None
        self.intype = None
        self.outtype = None

        self.storage = None

        # req'd when we pretend to be a simulator 
        self.n_frames_max = 10000


    def guess_fname_format(self, fname):
        '''Takes an file name and tries to guess the format from that (under
        the assumption that I'm using my normal convention that the last
        number before the .xyz is the trajectory number).
        
        Parameters
        ----------
        fname : str
            A filename that the appropriate format would output.
        
        Returns
        -------
        str
            Format string.'''
        num = re.search(".*[^0-9]([0-9]+)[^0-9]*.xyz", fname).group(1)
        numfmt = "%0"+str(len(num))+"d"
        fmt = re.sub(num, numfmt, fname)
        return fmt

    @staticmethod
    def build_parser(parser):
        '''Adds options to a OptionParser object which we'll need.
        
        Parameters
        ----------
        parser : optparse.OptionParser
        
        Returns
        -------
        optparse.OptionParser
        '''
        parser.add_option("-o", "--output", help="output file")
        parser.add_option("-t", "--topology", help="topology (PDB or XYZ)")
        return parser

    def set_infile_outfile(self, infiles, outfile):
        '''
        Identifies the list of files for input, output.

        Parameters
        ----------
        infiles : list of str
            Files to read as inputs
        outfile : string
            Outfile filename format string (i.e., allows formats like %06d
            for the trajectory number)
        '''
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
                self.outfile = outfile
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
        '''
        Wrapper to generate a topology from an arbitrary (xyz) file.

        Parameters
        ----------
        fname : str
            filename from which we extract the topology

        Returns
        -------
        Topology

        See also
        --------
        trajfile_topology : generates the actuall topology
        '''
        if re.match(".*\.xyz$", fname):
            topo_f = TrajFile.TrajFile()
            topo_f.read_xyz(fname)
            TrajFile.set_default_mass(topo_f)
            return self.trajfile_topology(topo_f)
        else:
            print "Topology must be .xyz for now"
        return None

    def trajfile_topology(self, trajfile):
        '''Creates a (fake-ish) OpenMM Topology from a TrajFile.

        Parameters
        ----------
        trajfile : TrajFile.TrajFile
            first frame from this trajectory generates a topology

        Returns
        -------
        Topology
        '''
        topol = Topology()
        # assume that the atoms are constant through the xyz trajectory, so
        # we can just use the first frame:
        myframe = trajfile.frames[0]
        chain = topol.addChain()
        for atom in range(myframe.natoms):
            # make each atom a separate element and a separate residue
            label = "_"+myframe.labels[atom] # underscore to avoid conflicts
            mass = myframe.mass[atom]
            # Check whether atom is already in the topol; _add_class it if not.
            # I hate this approach, since it assume the internal structure
            # of the Element classes
            try:
                element = Element.getBySymbol(label)
            except KeyError:
                element = Element(   
                                    number=atom+1, # abnormal
                                    name=label, symbol=label, mass=mass*amu
                                 )
            # we need to _add_class the elements to the mdtraj dictionaries, too,
            # because the default conversion functions between mdtraj and
            # openmm topologies don't actually check whether the element
            # dictionaries match up
            try:
                mdtrajelem = mdtraj.element.Element.getBySymbol(label)
            except KeyError:
                mdtrajelem = mdtraj.element.Element( 
                                number=atom+1,
                                name=label, symbol=label, mass=mass
                                )
            residue = topol.addResidue(label, chain)
            topol.addAtom(label, element, residue)
        return topol 

    def init_storage(self, fname):
        '''
        Initializes a storage object with appropriate information.

        Only used if we're converting from xyz to nc: requires that
        self.trajfile be loaded before calling.

        Parameters
        ----------
        fname : str
            filename for the netCDF output
        '''
        topol = self.trajfile_topology(self.trajfile)
        system = AtomCounter(self.trajfile.frames[0].natoms)
        self.simulation = SimulationDuckPunch(topol, system)
        self.storage = Storage( topology_file=topol,
                                          filename=fname, 
                                          mode='auto')
        snapshot.Snapshot.simulator = self
        self.storage.simulator = self
        self.storage.init_classes()
    

    def trajfile2trajectory(self, trajfile):
        '''Converts TrajFile.TrajFile to trajectory.Trajectory
        
        Parameters
        ----------
        trajfile : TrajFile.TrajFile
            The TrajFile object to be translated
        
        Returns
        -------
        Trajectory
            trajectory.Trajectory version of the trajfile; also writes the
            .nc file
        '''
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
        '''Converts trajectory.Trajectory to TrajFile.TrajFile
        
        Parameters
        ----------
        trajectory : trajectory:Trajectory
            The input trajectory to be translated

        Returns
        -------
        TrajFile
            TrajFile.TrajFile version of the trajectory
        '''
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

    def convert_xyz2nc(self):
        '''Entire conversion from a list of xyz files to a nc file. Includes
        loop over trajectories. Can only run when self.infiles and
        self.outfile are set.'''
        for f in self.infiles:
            self.trajfile = TrajFile.TrajFile()
            self.trajfile.read_xyz(f)
            TrajFile.set_default_mass(self.trajfile)
            self.trajectory = self.trajfile2trajectory(self.trajfile)

    def convert_nc2xyz(self):
        '''Entire conversion from a nc file to (possibly multiple) xyz
        files. Includes all trajectoris. Can only run when self.infiles,
        self.outfile, and self.topol_file are set.'''
        self.storage = Storage(
                                topology_file=None,
                                filename=self.infiles[0],
                                mode='restore'
                            )
        self.storage.simulator = self
        #self.storage.verbose_root = True # DEBUG
        self.storage._restore_options(self) # not used? (for topol)
        # for now, topology must come from a separate file
        self.storage.topology = self.topology_from_file(self.topol_file)
        print self.storage
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
            # check for a printf-formatted output file string. If not
            # using such a thing, subsequent trajectories *will*
            # overwrite in the same file
            if re.search("%[0-9]*d", self.outfile):
                myout = self.outfile % (traj_i+1)
            self.trajfile.write_xyz(myout)

    def translate(self):
        '''Convenience method which auto-selects the correct conversion
        method.'''
        if (self.intype == "xyz"):
            self.convert_xyz2nc()
        elif (self.intype == "nc"):
            self.convert_nc2xyz()


if __name__=="__main__":
    parser = optparse.OptionParser()
    translator = XYZTranslator()
    parser = translator.build_parser(parser)
    (opts, args) = parser.parse_args(sys.argv[1:])
    translator.set_infile_outfile(args, opts.output)
    translator.topol_file = opts.topology
    translator.translate()
    

