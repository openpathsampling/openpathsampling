'''
Created on 06.07.2014

@author: JDC Chodera
@author: JH Prinz
'''

import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python
import os.path

import numpy
import simtk.unit as units
import simtk.openmm.app
from simtk.unit import amu
import mdtraj as md

from object_storage import ObjectStorage
from trajectory_store import TrajectoryStorage, SampleStorage
from snapshot_store import SnapshotStorage, ConfigurationStorage, MomentumStorage
from ensemble_store import EnsembleStorage
from opentis.shooting import ShootingPointSelector, ShootingPoint
from opentis.pathmover import PathMover, MoveDetails
from opentis.globalstate import GlobalState
from orderparameter_store import ObjectDictStorage
from opentis.orderparameter import OrderParameter
from opentis.snapshot import Snapshot, Configuration

from opentis.storage.util import ObjectJSON


#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"

#=============================================================================================
# NetCDF Storage for multiple forked trajectories
#=============================================================================================

class Storage(netcdf.Dataset):
    '''
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories '''

    def __init__(self, filename='trajectory.nc', mode=None, n_atoms=None, 
                 topology_file=None, unit_system=None):
        '''
        Create a storage for complex objects in a netCDF file
        
        Parameters
        ----------        
        topology : openmm.app.Topology
            the topology of the system to be stored. Needed for 
        filename : string
            filename of the netcdf file
        mode : string, default: None
            the mode of file creation, one of 'w' (write), 'a' (append) or
            None, which will append any existing files.
        '''

        if mode == None:
            if os.path.isfile(filename):
                mode = 'a'
            else:
                mode = 'w'

        self.filename = filename
        self.links = []

        self.simplifier = ObjectJSON()

        if unit_system is not None:
            self.unit_system = unit_system
        else:
            self.unit_system = units.md_unit_system

        self.unit = dict()

        super(Storage, self).__init__(filename, mode)

        self.trajectory = TrajectoryStorage(self).register()
        self.snapshot = SnapshotStorage(self).register()
        self.configuration = ConfigurationStorage(self).register()
        self.momentum = MomentumStorage(self).register()
        self.ensemble = EnsembleStorage(self).register()
        self.sample = SampleStorage(self).register()
        self.pathmover = ObjectStorage(self, PathMover, named=True, json=True, identifier='json').register()
        self.movedetails = ObjectStorage(self, MoveDetails, named=False, json=True, identifier='json').register()
        self.shootingpoint = ObjectStorage(self, ShootingPoint, named=True, json=True).register()
        self.shootingpointselector = ObjectStorage(self, ShootingPointSelector, named=True, json=True, identifier='json').register()
        self.globalstate = ObjectStorage(self, GlobalState, named=True, json=True, identifier='json').register()
        self.collectivevariable = ObjectDictStorage(self, OrderParameter, Snapshot).register()
        self.cv = self.collectivevariable

        if mode == 'w':
            self._init()

            print topology_file

            if isinstance(topology_file, md.Topology):
                self.topology = topology_file

            elif isinstance(topology_file, simtk.openmm.app.Topology):
                self.topology = md.Topology.from_openmm(topology_file)

            elif type(topology_file) is str:
                self.pdb = md.load(topology_file)
                self.topology = self.pdb.topology
                self.pdb

            print self.topology

            # create a json from the mdtraj.Topology() and store it
            self.write_str('topology', self.simplifier.to_json(self.simplifier.topology_to_dict(self.topology)))

            self.initial_configuration = Configuration(
                coordinates=self.pdb.xyz[0],
                box_vectors=self.pdb.unitcell_vectors,
                potential_energy=0.0
            )

            # Save the initial configuration
            self.configuration.save(self.initial_configuration)

            self.createVariable('initial_configuration_idx', 'i4', 'scalar')
            self.variables['initial_configuration_idx'][:] = self.initial_configuration.idx[self]

            if n_atoms is not None:
                self.atoms = n_atoms
            elif self.topology is not None:
                self.atoms = self.topology.n_atoms
            else:
                raise RuntimeError("Storage given neither n_atoms nor topology")

            self._init_classes()
            self.sync()

        elif mode == 'a':
            self._restore_classes()

            self.topology = self.simplifier.topology_from_dict(self.simplifier.from_json(self.variables['topology'][0]))
            self.atoms = self.topology.n_atoms

            # restore initial configuration
            self.initial_configuration = self.configuration.load(int(self.variables['initial_configuration_idx'][0]))

            # Create a dict of simtk.Unit() instances for all netCDF.Variable()
            for variable_name in self.variables:
                unit = None
                variable = self.variables[variable_name]
                if hasattr(variable, 'unit_simtk'):

                    unit_dict = self.simplifier.from_json(getattr(variable, 'unit_simtk'))
                    if unit_dict is not None:
                        unit = self.simplifier.unit_from_dict(unit_dict)

                self.unit[str(variable_name)] = unit


    def get_unit(self, dimension):
        """
        Return a simtk.Unit instance from the unit_system the is of the specified dimension, e.g. length, time
        """
        return units.Unit({self.unit_system.base_units[units.BaseDimension(dimension)] : 1.0})

    def __getattr__(self, item):
        return self.__dict__[item]

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def _init_classes(self):
        '''
        Run the initialization on all added classes, when the storage is created only!

        Notes
        -----
        Only runs when the storage is created.
        '''

        for storage in self.links:
            # create a member variable which is the associated Class itself
            storage._init()

    def _restore_classes(self):
        '''
        Run restore on all added classes. Usually there is nothing to do.
        '''
#        for storage in self.links:
#            storage._restore()
        pass


    def _init(self):
        """
        Initialize the netCDF file for storage itself.
        """

        # add shared dimension for everyone. scalar and spatial
        if 'scalar' not in self.dimensions:
            self.createDimension('scalar', 1) # scalar dimension
            
        if 'spatial' not in self.dimensions:
            self.createDimension('spatial', 3) # number of spatial dimensions
        
        # Set global attributes.
        setattr(self, 'title', 'Open-Transition-Interface-Sampling')
        setattr(self, 'application', 'Host-Guest-System')
        setattr(self, 'program', 'run.py')
        setattr(self, 'programVersion', __version__)
        setattr(self, 'Conventions', 'Multi-State Transition Interface TPS')
        setattr(self, 'ConventionVersion', '0.1')

        # Create a string to hold the topology
        self.init_str('topology')
        self.write_str('topology', '')

        # Force sync to disk to avoid data loss.
        self.sync()

    def init_object(self, name):
        self.init_str(name)

    def write_as_json(self, name, obj):
        self.write_str(name, self.simplifier.to_json(obj))

    def restore_object(self, name):
        json_string = self.variables[name][0]
        return self.simplifier.from_json(json_string)

    def write_str(self, name, string):
        packed_data = numpy.empty(1, 'O')
        packed_data[0] = string
        self.variables[name][:] = packed_data

    def init_str(self, name):
        self.createVariable(name, 'str', 'scalar')


