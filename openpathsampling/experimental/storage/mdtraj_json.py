from ..simstore.custom_json import JSONCodec

try:
    import mdtraj as md
except ImportError:
    md = None
    HAS_MDTRAJ = False
else:
    HAS_MDTRAJ = True

import pandas as pd

def _check_mdtraj():
    if not HAS_MDTRAJ:
        raise RuntimeError("Unable to import MDTraj.")

def traj_to_dict(obj):
    return {'xyz': obj.xyz,
            'topology': obj.topology,
            'time': obj.time,
            'unitcell_lengths': obj.unitcell_lengths,
            'unitcell_angles': obj.unitcell_angles}

def traj_from_dict(dct):
    _check_mdtraj()
    dct = {k: v for k, v in dct.items()
           if k not in ['__class__', '__module__']}
    return md.Trajectory(**dct)

def topology_to_dict(obj):
    dataframe, bonds = obj.to_dataframe()
    return {'atoms': dataframe.to_json(),
            'bonds': bonds}

def topology_from_dict(dct):
    _check_mdtraj()
    return md.Topology.from_dataframe(
        atoms=pd.read_json(dct['atoms']),
        bonds=dct['bonds']
    )

if HAS_MDTRAJ:
    traj_codec = JSONCodec(md.Trajectory, traj_to_dict, traj_from_dict)
    top_codec = JSONCodec(md.Topology, topology_to_dict, topology_from_dict)
    mdtraj_codecs = [traj_codec, top_codec]
else:
    mdtraj_codecs = []

