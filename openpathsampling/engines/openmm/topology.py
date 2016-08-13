import mdtraj as md
import numpy as np
import pandas as pd
from simtk.openmm import XmlSerializer

from openpathsampling.engines import Topology


class MDTrajTopology(Topology):
    def __init__(self, mdtraj_topology, subsets=None):
        super(MDTrajTopology, self).__init__(int(mdtraj_topology.n_atoms), 3)
        self.mdtraj = mdtraj_topology
        if subsets is None:
            self.subsets = {}
        else:
            self.subsets = subsets

    def subset(self, list_of_atoms):
        return MDTrajTopology(self.mdtraj.subset(list_of_atoms), self.subsets)

    def to_dict(self):
        out = dict()
        # used_elements = set()

        atom_data = []
        for atom in self.mdtraj.atoms:
            if atom.element is None:
                element_symbol = ""
            else:
                element_symbol = atom.element.symbol

            atom_data.append((atom.serial, atom.name, element_symbol,
                         int(atom.residue.resSeq), atom.residue.name,
                         atom.residue.chain.index))
            # used_elements.add(atom.element)

        out['atom_columns'] = ["serial", "name", "element", "resSeq", "resName", "chainID"]
        out['atoms'] = atom_data
        out['bonds'] = [(a.index, b.index) for (a, b) in self.mdtraj.bonds]

        return {'mdtraj': out, 'subsets': self.subsets}

    @classmethod
    def from_dict(cls, dct):
        # TODO: fix this in a better way. Works for now with mdtraj 1.3.x and 1.4.x
        top_dict = dct['mdtraj']
        # elements = top_dict['elements']

        # for key, el in elements.iteritems():
        #     try:
        #         md.element.Element(
        #                     number=int(el[0]), name=el[1], symbol=el[2], mass=float(el[3])
        #                  )
        #         simtk.openmm.app.Element(
        #                     number=int(el[0]), name=el[1], symbol=el[2], mass=float(el[3])*units.amu
        #                  )
        #     except(AssertionError):
        #         pass

        atoms = pd.DataFrame(top_dict['atoms'], columns=top_dict['atom_columns'])
        bonds = np.array(top_dict['bonds'])

        md_topology = md.Topology.from_dataframe(atoms, bonds)

        return cls(md_topology, dct['subsets'])


class OpenMMSystemTopology(Topology):
    """A Topology that is based on an openmm.system object

    """
    def __init__(self, openmm_system, subsets=None):
        super(OpenMMSystemTopology, self).__init__(
            n_atoms=int(self.system.getNumParticles())
        )
        self.system = openmm_system
        if subsets is None:
            self.subsets = {}
        else:
            self.subsets = subsets

    def subset(self, list_of_atoms):
        return self

    def to_dict(self):
        system_xml = XmlSerializer.serialize(self.system)
        return {'system_xml' : system_xml, 'subsets' : self.subsets}

    @classmethod
    def from_dict(cls, dct):
        system_xml = dct['system_xml']
        subsets = dct['subsets']

        return cls(XmlSerializer.deserialize(system_xml), subsets)
