import json
import yaml
import base64
import simtk.unit as units
import mdtraj as md
import numpy as np
import pandas as pd
import simtk.openmm.app

import opentis.wrapper as wp

class ObjectJSON(object):
    """
    A simple implementation of a pickle algorithm to create object that can be converted to json and back
    """

    def __init__(self, unit_system = None, class_list = None):
        self.excluded_keys = ['name']
        self.unit_system = unit_system
        if class_list is not None:
            self.class_list = class_list
        else:
            self.class_list = wp.class_list

    def simplify_object(self, obj, base_type = ''):
        return { '_cls' : obj.__class__.__name__, '_dict' : self.simplify(obj.to_dict(), obj.base_cls_name) }

    def simplify(self, obj, base_type = ''):
        if type(obj).__module__ != '__builtin__':
            if type(obj) is units.Quantity:
                # This is number with a unit so turn it into a list
                if self.unit_system is not None:
                    return { '_value' : obj.value_in_unit_system(self.unit_system), '_units' : self.unit_to_dict(obj.unit.in_unit_system(self.unit_system)) }
                else:
                    return { '_value' : obj / obj.unit, '_units' : self.unit_to_dict(obj.unit) }
            elif type(obj) is np.ndarray:
                # this is maybe not the best way to store large numpy arrays!
                return { '_numpy' : self.simplify(obj.shape), '_dtype' : str(obj.dtype), '_data' : base64.b64encode(obj) }
            elif hasattr(obj, 'to_dict'):
                # the object knows how to dismantle itself into a json string so use this
                return { '_cls' : obj.__class__.__name__, '_dict' : self.simplify(obj.to_dict(), obj.base_cls_name) }
            else:
                return None
        elif type(obj) is list:
            return [self.simplify(o, base_type) for o in obj]
        elif type(obj) is tuple:
            return tuple([self.simplify(o, base_type) for o in obj])
        elif type(obj) is dict:
            result = {key : self.simplify(o, base_type) for key, o in obj.iteritems() if type(key) is str and key not in self.excluded_keys }
            return result
        elif type(obj) is slice:
            return { '_slice' : [obj.start, obj.stop, obj.step]}
        else:
            oo = obj
            return oo

    def build(self,obj):
        global class_list
        if type(obj) is dict:
            if '_units' in obj and '_value' in obj:
                return obj['_value'] * self.unit_from_dict(obj['_units'])
            elif '_slice' in obj:
                return slice(*obj['_slice'])
            elif '_numpy' in obj:
                return np.frombuffer(base64.decodestring(obj['_data']), dtype=np.dtype(obj['_dtype'])).reshape(tuple(obj['_numpy']))
            elif '_cls' in obj and '_dict' in obj:
                if obj['_cls'] in self.class_list:
                    attributes = self.build(obj['_dict'])
                    return self.class_list[obj['_cls']].from_dict(attributes)
                else:
                    raise ValueError('Cannot create obj of class "' + obj['_cls']+ '". Class is not registered as creatable!')
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

            atom_data.append((int(atom.serial), atom.name, element_symbol,
                         int(atom.residue.resSeq), atom.residue.name,
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

    def to_json(self, obj, base_type = ''):
        simplified = self.simplify(obj, base_type)
        return json.dumps(simplified)

    def to_json_object(self, obj, base_type = ''):
        simplified = self.simplify_object(obj, base_type)
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
    def __init__(self, storage, unit_system = None, class_list = None):
        super(StorableObjectJSON, self).__init__(unit_system, class_list)
        self.excluded_keys = ['name', 'idx', 'json', 'identifier']
        self.storage = storage

    def simplify(self,obj, base_type = ''):
        if type(obj).__module__ != '__builtin__':
#            print obj.__dict__, hasattr(obj, 'creatable')
            if hasattr(obj, 'idx') and (not hasattr(obj, 'nestable') or (obj.base_cls_name != base_type)):
                # this also returns the base class name used for storage
                # store objects only if they are not creatable. If so they will only be created in their
                # top instance and we use the simplify from the super class ObjectJSON
                base_cls = self.storage.save(obj)
                return { '_idx' : obj.idx[self.storage], '_base' : base_cls, '_cls' : obj.__class__.__name__ }

        return super(StorableObjectJSON, self).simplify(obj, base_type)

    def build(self,obj):
        if type(obj) is dict:
            if '_base' in obj and '_idx' in obj:
                result = self.storage.load(obj['_base'], obj['_idx'])
                # restore also the actual class name

                result.cls = obj['_cls']
                return result

        return super(StorableObjectJSON, self).build(obj)