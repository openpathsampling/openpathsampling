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

        atom_data = []
        for atom in self.mdtraj.atoms:
            if atom.element is None:
                element_symbol = ""
            else:
                element_symbol = atom.element.symbol

            atom_data.append((
                atom.serial, atom.name, element_symbol,
                int(atom.residue.resSeq), atom.residue.name,
                atom.residue.chain.index, atom.segment_id))

        out['atom_columns'] = ["serial", "name", "element", "resSeq",
                               "resName", "chainID", "segmentID"]
        out['atoms'] = atom_data
        out['bonds'] = [(a.index, b.index) for (a, b) in self.mdtraj.bonds]

        return {'mdtraj': out}

    @classmethod
    def from_dict(cls, dct):
        top_dict = dct['mdtraj']

        atoms = pd.DataFrame(top_dict['atoms'], columns=top_dict['atom_columns'])
        bonds = np.array(top_dict['bonds'])

        try:
            md_topology = md.Topology.from_dataframe(atoms, bonds)
            return cls(md_topology)
        except StandardError:
            # we try a fix and add multiples of 10000 to the resSeq

            for ci in np.unique(atoms['chainID']):
                chain_atoms = atoms[atoms['chainID'] == ci]
                indices = chain_atoms.index.tolist()

                old_residue_id = 0
                multiplier = 0
                places = []
                for row, res_id in zip(indices, list(chain_atoms['resSeq'])):
                    if res_id < old_residue_id:
                        if multiplier > 0:
                            atoms.loc[places, 'resSeq'] += 10000 * multiplier

                        places = []
                        multiplier += 1

                    if multiplier > 0:
                        places.append(row)

                    old_residue_id = res_id

                if multiplier > 0:
                    atoms.loc[places, 'resSeq'] += 10000 * multiplier


            # this function is really slow! Reads ~ 1000 atoms per second
            md_topology = md.Topology.from_dataframe(atoms, bonds)

            # that we have successfully created the topology using from_df
            # we remove the wrong multipliers
            # this is weird, but reproduces the current behaviour

            for atom in md_topology.atoms:
                atom.residue.resSeq %= 10000

            return cls(md_topology)


class OpenMMSystemTopology(Topology):
    """A Topology that is based on an openmm.system object

    This uses the XmlSerializer class from OpenMM itself to transform the
    system object into an XML string which is then stored in the JSON.
    This is rather inefficient but works very stable.

    """
    def __init__(self, openmm_system):
        super(OpenMMSystemTopology, self).__init__(
            n_atoms=int(self.system.getNumParticles())
        )
        self.system = openmm_system

    def to_dict(self):
        system_xml = XmlSerializer.serialize(self.system)
        return {'system_xml': system_xml}

    @classmethod
    def from_dict(cls, dct):
        system_xml = dct['system_xml']

        return cls(XmlSerializer.deserialize(system_xml))
