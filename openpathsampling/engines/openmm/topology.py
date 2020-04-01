import mdtraj as md
import numpy as np
import pandas as pd

from openpathsampling.engines import Topology

import logging
logger = logging.getLogger(__name__)


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

            if hasattr(atom, 'segment_id'):
                atom_data.append((
                    atom.serial, atom.name, element_symbol,
                    int(atom.residue.resSeq), atom.residue.name,
                    atom.residue.chain.index, atom.segment_id))
            else:
                atom_data.append((
                    atom.serial, atom.name, element_symbol,
                    int(atom.residue.resSeq), atom.residue.name,
                    atom.residue.chain.index, ''))

        out['atom_columns'] = ["serial", "name", "element", "resSeq",
                               "resName", "chainID", "segmentID"]
        out['atoms'] = atom_data
        out['bonds'] = [(a.index, b.index) for (a, b) in self.mdtraj.bonds]

        return {'mdtraj': out}

    @classmethod
    def from_dict(cls, dct):
        top_dict = dct['mdtraj']

        atoms = pd.DataFrame(
            top_dict['atoms'], columns=top_dict['atom_columns'])
        bonds = np.array(top_dict['bonds'])

        try:
            md_topology = md.Topology.from_dataframe(atoms, bonds)
            return cls(md_topology)
        except Exception:
            # we try a fix and add multiples of 10000 to the resSeq

            logger.info('Normal reconstruction of topology failed. '
                        'Trying a fix to the 10k residue ID problem.')

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

