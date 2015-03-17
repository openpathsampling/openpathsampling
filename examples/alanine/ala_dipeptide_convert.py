"""
This is an example file of how to remove water from a netCDF file.

@author: Jan-Hendrik Prinz
"""

from openpathsampling.storage import Storage

if __name__=="__main__":

    storage = Storage(
        filename="trajectory.nc",
        mode='a'
    )

    solute = range(22)
    storage.clone('test_solute.nc', subset=None)
    storage.close()

    storage2 = Storage('test_solute.nc')
    storage2.clone('test_solute2.nc', subset=solute)