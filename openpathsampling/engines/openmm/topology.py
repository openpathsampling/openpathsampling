import mdtraj as md
import numpy as np
import pandas as pd
from simtk.openmm import XmlSerializer

from openpathsampling.engines import Topology


class MDTrajTopology(Topology):
    def __init__(self, mdtraj_topology):
        super(MDTrajTopology, self).__init__(int(mdtraj_topology.n_atoms), 3)
        self.mdtraj = mdtraj_topology

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

        return {'mdtraj': out}

    @classmethod
    def from_dict(cls, dct):
        # TODO: fix this in a better way. Works for now with mdtraj 1.3.x and 1.4.x
        top_dict = dct['mdtraj']

        atoms = pd.DataFrame(top_dict['atoms'], columns=top_dict['atom_columns'])
        bonds = np.array(top_dict['bonds'])

        md_topology = md.Topology.from_dataframe(atoms, bonds)

        return cls(md_topology)


class OpenMMSystemTopology(Topology):
    """A Topology that is based on an openmm.system object

    """
    def __init__(self, openmm_system):
        super(OpenMMSystemTopology, self).__init__(
            n_atoms=int(self.system.getNumParticles())
        )
        self.system = openmm_system

    def to_dict(self):
        system_xml = XmlSerializer.serialize(self.system)
        return {'system_xml' : system_xml}

    @classmethod
    def from_dict(cls, dct):
        system_xml = dct['system_xml']

        return cls(XmlSerializer.deserialize(system_xml))
