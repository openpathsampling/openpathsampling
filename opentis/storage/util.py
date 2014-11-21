import json
import yaml
import simtk.unit as units
import mdtraj as md
import numpy as np
import pandas as pd
import simtk.openmm.app

class Simplifier(object):
    """
    A simple implementation of a pickle algorithm to create object that can be converted to json and back
    """
    def __init__(self):
        self.excluded_keys = []

    def simplify(self,obj):
        if type(obj).__module__ != '__builtin__':
            if type(obj) is units.Quantity:
                # This is number with a unit so turn it into a list
                return { '_value' : obj / obj.unit, '_units' : self.unit_to_dict(obj.unit) }
            else:
                return None
        elif type(obj) is list:
            return [self.simplify(o) for o in obj]
        elif type(obj) is tuple:
            return tuple([self.simplify(o) for o in obj])
        elif type(obj) is dict:
            return {key : self.simplify(o) for key, o in obj.iteritems() if type(key) is str and key not in self.excluded_keys}
        else:
            return obj

    def build(self,obj):
        if type(obj) is dict:
            if '_units' in obj and '_value' in obj:
                return obj['_value'] * self.dict_to_unit(obj['_units'])
            else:
                return {key : self.build(o) for key, o in obj.iteritems()}
        elif type(obj) is tuple:
            return tuple([self.build(o) for o in obj])
        elif type(obj) is list:
            return [self.build(o) for o in obj]
        else:
            return obj

    def to_json(self, obj):
        simplified = self.simplify(obj)
        return json.dumps(simplified)

    def from_json(self, json_string):
        simplified = yaml.load(json_string)
        return self.build(simplified)

    def unitsytem_to_list(self, unit_system):
        '''
        Turn a simtk.UnitSystem() into a list of strings representing the unitsystem for serialization
        '''
        return [ u.name  for u in unit_system.units ]

    def unit_system_from_list(self, unit_system_list):
        '''
        Create a simtk.UnitSystem() from a serialialized list of strings representing the unitsystem
        '''
        return units.UnitSystem([ getattr(units, unit_name).iter_all_base_units().next()[0] for unit_name in unit_system_list])

    def unit_to_symbol(self, unit):
        return str(1.0 * unit).split()[1]

    def unit_to_dict(self, unit):
        d = {p.name : int(fac) for p, fac in unit.iter_all_base_units()}
        return d

    def dict_to_unit(self, unit_dict):
        unit = units.Unit({})
        for unit_name, unit_multiplication in unit_dict.iteritems():
            unit *= getattr(units, unit_name)**unit_multiplication

        return unit

    def topology_to_dict(self, topology):
        """Return a copy of the topology

        Returns
        -------
        out : Topology
            A copy of this topology
        """

        out = dict()
        used_elements = set()

        atom_data = []
        for atom in topology.atoms:
            if atom.element is None:
                element_symbol = ""
            else:
                element_symbol = atom.element.symbol

            atom_data.append((atom.serial, atom.name, element_symbol,
                         atom.residue.resSeq, atom.residue.name,
                         atom.residue.chain.index))

            used_elements.add(atom.element)

        out['atom_columns'] = ["serial", "name", "element", "resSeq", "resName", "chainID"]
        out['atoms'] = atom_data
        out['bonds'] = [(a.index, b.index) for (a, b) in topology.bonds]
        out['elements'] = {key: tuple(el) for key, el in md.element.Element._elements_by_symbol.iteritems() if el in used_elements}

        return out

    def topology_from_dict(self, top_dict):
        elements = top_dict['elements']

        for key, el in elements.iteritems():
            try:
                md.element.Element(
                            number=int(el[0]), name=el[1], symbol=el[2], mass=float(el[3])
                         )
                simtk.openmm.app.Element(
                            number=int(el[0]), name=el[1], symbol=el[2], mass=float(el[3])*units.amu
                         )
            except(AssertionError):
                pass


        atoms = pd.DataFrame(top_dict['atoms'], columns=top_dict['atom_columns'])
        bonds = np.array(top_dict['bonds'])

        return md.Topology.from_dataframe(atoms, bonds)

class ObjectSimplifier(Simplifier):
    def __init__(self):
        super(ObjectSimplifier, self).__init__()
        self.excluded_keys = ['idx', 'json', 'identifier']

    def simplify(self,obj):
        if type(obj).__module__ != '__builtin__':
            if hasattr(obj, 'cls'):
                getattr(self.storage, obj.cls).save(obj)
                return { 'idx' : obj.idx[self.storage], 'cls' : obj.cls}

        super(ObjectSimplifier, self).simplifiy(obj)

    def build(self,obj):
        if type(obj) is dict:
            if 'cls' in obj and 'idx' in obj:
                return getattr(self.storage, obj['cls']).load(obj['idx'])

        super(ObjectSimplifier, self).build(obj)