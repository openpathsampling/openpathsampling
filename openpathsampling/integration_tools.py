"""
Tools for integration with miscellaneous non-required packages.
"""

import importlib
import logging


def error_if_no(name, package_name, has_package):
    if not has_package:
        raise RuntimeError(name + " requires " + package_name
                           + ", which is not installed")

def _chain_import(*packages):
    """
    Import as whichever name is first importable among ``packages``.

    If none exists, raises error based on last package attempted
    """
    error = None
    for package in packages:
        try:
            pkg = importlib.import_module(package)
        except ImportError as e:
            error = e
        else:
            return pkg
    # we raise the last error given
    raise error

# openmm.unit ########################################################
try:
    unit = _chain_import('openmm.unit', 'simtk.unit')
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
    return error_if_no(name, "openmm.unit or simtk.unit", HAS_SIMTK_UNIT)

# mdtraj ############################################################
try:
    # MDTraj currently imports OpenMM from the simtk namespace, leading
    # to warnings being issued that we can't control (and cause our
    # notebook tests for fail). So we need to disable here.
    logging.disable(logging.WARNING)
    import mdtraj as md
    logging.disable(logging.NOTSET)
    # The problem with this is that it will shadow any remaining places
    # we're having this problem -- the simtk import is only done once.
except ImportError:
    md = None
    HAS_MDTRAJ = False
else:
    HAS_MDTRAJ = True

def error_if_no_mdtraj(name):
    return error_if_no(name, "mdtraj", HAS_MDTRAJ)

# openmm ############################################################
try:
    openmm = _chain_import('openmm', 'simtk.openmm')
except ImportError:
    openmm = None
    HAS_OPENMM = False
else:
    HAS_OPENMM = True

def error_if_to_openmm(name):
    return error_if_no(name, "openmm", HAS_OPENMM)
