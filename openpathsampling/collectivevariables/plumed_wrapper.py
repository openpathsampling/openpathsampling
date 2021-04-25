import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject
import openpathsampling.engines as peng
try:
    import plumed
except ImportError:
    pass
import numpy as np
import warnings
import copy
import os
import sys


class PLUMEDCV(paths.collectivevariable.CoordinateFunctionCV):

    """Make `CollectiveVariable` computed by PLUMED [1]_ according to the
    command `name`: `definition`, where `name` is a PLUMED label and
    `definition` contains all PLUMED keywords.
    Takes an `openpathsampling.engines.trajectory.Trajectory` as input.

    References
    ----------

    .. [1] G.A. Tribello, M. Bonomi, D. Branduardi, C. Camilloni, G. Bussi,
       PLUMED2: New feathers for an old bird, Comp. Phys. Comm. 185, 604
       (2014); https://doi.org/10.1016/j.cpc.2013.09.018

    Examples
    --------
    >>> # To create a `CollectiveVariable` which calculates the dihedral psi
    >>> # formed by the atoms [7,9,15,17] in Ala dipeptide:
    >>> from openpathsampling import PLUMEDCV, PLUMEDInterface
    >>> plmd = PLUMEDInterface(top)
    >>> # top is an `openpathsampling.engines.topology.MDTrajTopology`
    >>> psi_plumed = PLUMEDCV("psi",plmd,"TORSION ATOMS=7,9,15,17")
    >>> print psi_plumed(traj)  # returns psi values for the trajectory
    """

    def __init__(self,
                 name,
                 plmd,
                 definition,
                 components=None,
                 cv_requires_lists=True,
                 cv_wrap_numpy_array=True,
                 cv_scalarize_numpy_singletons=True,
                 **kwargs
                 ):

        """
        Parameters
        ----------
        name : string
            A descriptive name of the PLUMED collective variable,
            equivalent to a `label` in a PLUMED input file.
        plmd : :obj:`openpathsampling.collectivevariable.PLUMEDInterface`
            An interface to the Cython PLUMED wrapper. If the PLUMED collective
            variable is a function of previously defined ones, or if it is defined
            based on group/virtual atoms, the `plmd` interface must be the same
            one that was used for the preceding instantiantions.
        definition : string
            The PLUMED keywords that define the collective variable
            (see http://www.plumed.org/documentation).
        components : list of string
            The components (either default of customized) of the PLUMED
            collective variable (see http://www.plumed.org/documentation).
        cv_requires_lists
        cv_wrap_numpy_array
        cv_scalarize_numpy_singletons
        kwargs
        """

        super(PLUMEDCV, self).__init__(
                name,
                f=PLUMEDCV.compute_cv,
                cv_requires_lists=cv_requires_lists,
                cv_wrap_numpy_array=cv_wrap_numpy_array,
                cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
                **kwargs
                )

        self.plmd = plmd
        self.definition = definition
        if components is None:
            components = []
        self.components = components
        self.topology = plmd.topology
        self.var = self.create_plumed_var(name, definition)

    def create_plumed_var(self, name, definition):
        """Create a PLUMED collective variable.

        Parameters
        ----------
        name : string
            A descriptive name of the PLUMED collective variable,
            equivalent to a `label` in a PLUMED input file.
        definition : string
            The PLUMED keywords that define the collective variable
            (see http://www.plumed.org/documentation).

        Returns
        -------
        data : array
            Array to store the computed PLUMED collective variable.
        """
        self.plmd.cmd("readInputLine", name + ": " + definition)
        if (len(self.components)):
            cdata = []
            for c in self.components:
                shape = np.zeros(1, dtype=np.int_)
                self.plmd.cmd("getDataRank " + name + "." + c, shape)
                data = np.zeros((1))
                self.plmd.cmd("setMemoryForData " + name + "." + c, data)
                cdata.append(data)
            return cdata
        else:
            shape = np.zeros(1, dtype=np.int_)
            self.plmd.cmd("getDataRank " + name, shape)
            data = np.zeros((1))
            self.plmd.cmd("setMemoryForData " + name, data)
            return data

    def compute_cv(self, trajectory):
        """Compute a PLUMED collective variable.

        Parameters
        ----------
        trajectory : :obj:`openpathsampling.engines.trajectory.Trajectory`
            The trajectory along which the collective variable is to be
            computed.

        Returns
        -------
        cv : array
            Computed values of the PLUMED collective variable along the
            `openpathsampling.engines.trajectory.Trajectory`
        """
        cv = []
        try:
            masses = np.array([a.element.mass for a in
                              trajectory.topology.mdtraj.atoms],
                              dtype=np.float64)
        except AttributeError:  # pragma: no cover
            masses = np.ones(self.topology.n_atoms, dtype=np.float64)
            warnings.warn("No masses found in topology. All masses set to one")
        charges = np.zeros(self.topology.n_atoms,
                           dtype=np.float64)  # non-essential
        # warnings.warn("All charges set to zero")
        forces = np.zeros((self.topology.n_atoms, 3))
        virial = np.zeros((3, 3), dtype=np.float64)
        bias = np.zeros((1), dtype=np.float64)  # non-essential
        for step, snapshot in enumerate(trajectory):
            self.plmd.cmd("setStep", step)
            if snapshot.box_vectors != None:
                box = np.array(snapshot.box_vectors, dtype=np.float64)
                self.plmd.cmd("setBox", box)
            positions = snapshot.xyz.astype(np.float64)
            self.plmd.cmd("setBox", box)
            self.plmd.cmd("setPositions", positions)
            self.plmd.cmd("setMasses", masses)
            self.plmd.cmd("setForces", forces)
            self.plmd.cmd("setVirial", virial)
            self.plmd.cmd("setCharges", charges)  # non-essential
            self.plmd.cmd("getBias", bias)   # non-essential
            self.plmd.cmd("calc")
            cv.append(copy.deepcopy(self.var))
        cv = np.array(cv)
        return cv

    def _eval(self, trajectory):
        trajectory = peng.Trajectory(trajectory)
        return self.cv_callable(self, trajectory)

    def to_dict(self):
        return {
            'name': self.name,
            'plmd': self.plmd,
            'definition': self.definition,
            'components': self.components,
            'kwargs': self.kwargs,
            'cv_requires_lists': self.cv_requires_lists,
            'cv_wrap_numpy_array': self.cv_wrap_numpy_array,
            'cv_scalarize_numpy_singletons': self.cv_scalarize_numpy_singletons
        }


class PLUMEDInterface(StorableNamedObject):

    """Interfaces the Cython PLUMED wrapper [1]_ located at
    `/path/to/plumed2/python` and allows to set and get non-`PLUMEDCV`
    commands (i.e., non-outputting).  This includes groups of atoms, centers
    of mass, include files, etc.  Requires PLUMED development version (see
    https://github.com/plumed/plumed2) and sourcing
    `/path/to/plumed2/sourceme.sh`.

    References
    ----------

    .. [1] G.A. Tribello, M. Bonomi, D. Branduardi, C. Camilloni, G. Bussi,
       PLUMED2: New feathers for an old bird, Comp. Phys. Comm. 185, 604
       (2014); https://doi.org/10.1016/j.cpc.2013.09.018

    Examples
    --------
    >>> # To group the atoms [7,9,15,17] corresponding to the dihedral psi
    >>> # in Ala dipeptide:
    >>> from openpathsampling import PLUMEDCV, PLUMEDInterface
    >>> plmd = PLUMEDInterface(top)
    >>> # top is an `openpathsampling.engines.topology.MDTrajTopology`
    >>> plmd.set("group","GROUP ATOMS=7,9,15,17")
    >>> psi_plumed = PLUMEDCV("psi",plmd,"TORSION ATOMS=group")
    >>> print psi_plumed(traj)  # returns psi values for the trajectory
    >>> pld.get()  # returns (('group', 'GROUP ATOMS=7,9,15,17'))

    """

    def __init__(self,
                 topology,
                 pathtoplumed="",
                 timestep=1.,
                 kbt=1.,
                 molinfo="",
                 logfile="plumed.log"):
        """
        Parameters
        ----------
        topology : :obj:`openpathsampling.engines.topology.MDTrajTopology`
        pathtoplumed : string
            path to the PLUMED installation
        timestep : double
            Time step size of the simulation in PLUMED default units (ps).
        kbt : double
            :math:`$k_BT$` in PLUMED default units (kJ/mol).
        molinfo : string
            A PDB file containing information about the molecule. (see
            https://plumed.github.io/doc-v2.4/user-doc/html/_m_o_l_i_n_f_o.html).
        logfile : string
            Name of the PLUMED log file.
        """

        self.interface = plumed.Plumed()  # 8 is default size of real
        self.pathtoplumed = pathtoplumed
        self.topology = topology
        self.timestep = timestep
        self.kbt = kbt
        self.molinfo = molinfo
        self.logfile = logfile
        self._commandlist = []
        self._init_plumed()

    def cmd(self, *args, **kwargs):
        self.interface.cmd(*args, **kwargs)

    def _init_plumed(self):
        if self.pathtoplumed != "":  # pragma: no cover
            #os.system("source " + self.pathtoplumed + "/sourceme.sh")
            sys.path.append(self.pathtoplumed + "/python")
            warnings.warn("Sourced PLUMED from: " + self.pathtoplumed)
        else:
            warnings.warn("Using currently sourced PLUMED from: " +
            plumed.sys.prefix + "/lib/libplumedKernel.dylib")
            #os.environ["PLUMED_KERNEL"][:-30]) not set in conda-forge plumed
        self.cmd("setMDEngine", "python")
        self.cmd("setTimestep", self.timestep)
        self.cmd("setKbT", self.kbt)
        self.cmd("setNatoms", self.topology.n_atoms)
        self.cmd("setLogFile", self.logfile)
        self.cmd("init")
        if (self.molinfo != ""):
            self.cmd("readInputLine", "MOLINFO STRUCTURE=" + self.molinfo)

    def set(self, name, definition):
        """Set a non-outputting command in the `PLUMEDInterface`.

        Parameters
        ----------
        name : string
            Equivalent to a `label` in a PLUMED input file. Not required for
            all commands (can be left empty).
        definition: string
            The PLUMED keywords that define the command
            (see http://www.plumed.org/documentation).
        """

        if (name != ""):
            self.cmd("readInputLine", name + ": " + definition)
        else:
            self.cmd("readInputLine", definition)
        self._commandlist.append((name, definition))

    def get(self):
        """Get commands set in the `PLUMEDInterface`.

        Returns
        -------
        list of tuples
            list of tuples (`label` and `definition`) ran in the
            `PLUMEDInterface` using the `set` function.
        """
        return tuple(self._commandlist)

    def to_dict(self):
        return {
            'pathtoplumed': self.pathtoplumed,
            'topology': self.topology,
            'timestep': self.timestep,
            'kbt': self.kbt,
            'molinfo': self.molinfo,
            'logfile': self.logfile,
            '_commandlist': self._commandlist
        }

    @classmethod
    def from_dict(cls, dct):
        pathtoplumed = dct['pathtoplumed']
        topology = dct['topology']
        timestep = dct['timestep']
        kbt = dct['kbt']
        molinfo = dct['molinfo']
        logfile = dct['logfile']
        obj = cls(pathtoplumed=pathtoplumed,
                  topology=topology,
                  timestep=timestep,
                  kbt=kbt,
                  molinfo=molinfo,
                  logfile=logfile)
        _commandlist = dct['_commandlist']
        for name, definition in _commandlist:
            obj.set(name, definition)
        return obj
