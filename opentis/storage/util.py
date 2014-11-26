import json
import yaml
import simtk.unit as units
import mdtraj as md
import numpy as np
import pandas as pd
import simtk.openmm.app

class ObjectJSON(object):
    """
    A simple implementation of a pickle algorithm to create object that can be converted to json and back
    """
    def __init__(self, unit_system = None):
        self.excluded_keys = []
        self.unit_system = unit_system

    def simplify(self, obj):
        if type(obj).__module__ != '__builtin__':
            if type(obj) is units.Quantity:
                # This is number with a unit so turn it into a list
                if self.unit_system is not None:
                    return { '_value' : obj.value_in_unit_system(self.unit_system), '_units' : self.unit_to_dict(obj.unit.in_unit_system(self.unit_system)) }
                else:
                    return { '_value' : obj / obj.unit, '_units' : self.unit_to_dict(obj.unit) }
            else:
                return None
        elif type(obj) is list:
            return [self.simplify(o) for o in obj]
        elif type(obj) is tuple:
            return tuple([self.simplify(o) for o in obj])
        elif type(obj) is dict:
            result = {key : self.simplify(o) for key, o in obj.iteritems() if type(key) is str and key not in self.excluded_keys }
            return result
        else:
            oo = obj
            return oo


    def build(self,obj):
        if type(obj) is dict:
            if '_units' in obj and '_value' in obj:
                return obj['_value'] * self.unit_from_dict(obj['_units'])
            else:
                return {key : self.build(o) for key, o in obj.iteritems()}
        elif type(obj) is tuple:
            return tuple([self.build(o) for o in obj])
        elif type(obj) is list:
            return [self.build(o) for o in obj]
        else:
            return obj

    def unitsytem_to_list(self, unit_system):
        '''
        Turn a simtk.UnitSystem() into a list of strings representing the unitsystem for serialization
        '''
        return [ u.name  for u in unit_system.units ]

    def unit_system_from_list(self, unit_system_list):
        '''
        Create a simtk.UnitSystem() from a serialialized list of strings representing the unitsystem
        '''
        return units.UnitSystem([ getattr(units, unit_name).iter_base_or_scaled_units().next()[0] for unit_name in unit_system_list])

    def unit_to_symbol(self, unit):
        return str(1.0 * unit).split()[1]

    def unit_to_dict(self, unit):
        unit_dict = {p.name : int(fac) for p, fac in unit.iter_base_or_scaled_units()}
        return unit_dict

    def unit_from_dict(self, unit_dict):
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

    def to_json(self, obj):
        simplified = self.simplify(obj)
        return json.dumps(simplified)

    def from_json(self, json_string):
        simplified = yaml.load(json_string)
        return self.build(simplified)

    def topology_to_json(self, topology):
        return self.to_json(self.topology_to_dict(topology))

    def topology_from_json(self, json_string):
        return self.topology_from_dict(self.from_json(json_string))

    def unit_to_json(self, unit):
        simple = self.unit_to_dict(unit)
        return self.to_json(simple)

    def unit_from_json(self, json_string):
        return self.unit_from_dict(self.from_json(json_string))


class StorableObjectJSON(ObjectJSON):
    def __init__(self, storage, unit_system = None):
        super(StorableObjectJSON, self).__init__(unit_system)
        self.excluded_keys = ['idx', 'json', 'identifier']
        self.storage = storage

    def simplify(self,obj):

        if type(obj).__module__ != '__builtin__':
            if hasattr(obj, 'idx'):
#                getattr(self.storage, obj.cls).save(obj)
                # this also return the base class name used for storage
                base_cls = self.storage.save(obj)
#                return { '_idx' : obj.idx[self.storage], '_cls' : obj.cls}
                return { '_idx' : obj.idx[self.storage], '_base' : base_cls, '_cls' : obj.__class__.__name__ }


        return super(StorableObjectJSON, self).simplify(obj)

    def build(self,obj):
        if type(obj) is dict:
            if '_base' in obj and '_idx' in obj:
                result = self.storage.load(obj['_base'], obj['_idx'])
                # restore also the actual class name

                result.cls = obj['_cls']
                return result
#                return getattr(self.storage, obj['_cls']).load(obj['_idx'])

        return super(StorableObjectJSON, self).build(obj)