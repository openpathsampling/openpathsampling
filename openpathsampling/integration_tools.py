"""
Tools for integration with miscellaneous non-required packages.
"""

def error_if_no(name, package_name, has_package):
    if not has_package:
        raise RuntimeError(name + " requires " + package_name
                           + ", which is not installed")

# simtk.unit ########################################################
try:
    from simtk import unit
except ImportError:
    unit = None
    is_simtk_quantity = lambda obj: False
    is_simtk_quantity_type = lambda obj: False
    is_simtk_unit_type = lambda obj: False
    HAS_SIMTK_UNIT = False
else:
    is_simtk_quantity = lambda obj: obj.__class__ is unit.Quantity
    is_simtk_quantity_type = lambda obj: type(obj) is unit.Quantity
    is_simtk_unit_type = lambda obj: type(obj) is unit.Unit
    HAS_SIMTK_UNIT = True

def error_if_no_simtk_unit(name):
    return error_if_no(name, "simtk.unit", HAS_SIMTK_UNIT)

# mdtraj ############################################################
try:
    import mdtraj as md
except ImportError:
    md = None
    HAS_MDTRAJ = False
else:
    HAS_MDTRAJ = True

def error_if_no_mdtraj(name):
    return error_if_no(name, "mdtraj", HAS_MDTRAJ)

# openmm ############################################################
try:
    from simtk import openmm
except ImportError:
    openmm = None
    HAS_OPENMM = False
else:
    HAS_OPENMM = True

def error_if_to_openmm(name):
    return error_if_no(name, "openmm", HAS_OPENMM)
