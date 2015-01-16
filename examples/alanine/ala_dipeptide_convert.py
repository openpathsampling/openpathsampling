"""
This is an example file of how to remove water from a netCDF file.

@author: Jan-Hendrik Prinz
"""

from opentis.storage import Storage

if __name__=="__main__":

    storage = Storage(
        filename="test_simple.nc",
        mode='a'
    )

    solute = range(22)
    storage.clone('test_solute.nc', subset=solute)