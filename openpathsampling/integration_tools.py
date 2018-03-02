"""
Tools for integration with miscellaneous non-required packages.
"""

try:
    from simtk import unit as units
except ImportError:
    is_simtk_quantity = lambda obj: False
    HAS_SIMTK_UNIT = False
else:
    is_simtk_quantity = lambda obj: obj.__class__ is units.Quantity
    HAS_SIMTK_UNIT = True

def error_if_no_simtk_unit(name):
    if not HAS_SIMTK_UNIT:
        raise RuntimeError(name + " requires simtk.unit, which is "
                           + "not installed")

