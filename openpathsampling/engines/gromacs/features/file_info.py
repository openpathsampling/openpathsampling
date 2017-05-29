"""
Attributes
----------
file_number : int
file_position : int
"""

variables = ['file_number', 'file_position']

def netcdfplus_init(store):
    store.create_variable(
        'file_number',
        'int',
        description="the file number that this snapshot is stored in"
    )
    store.create_variable(
        'file_position',
        'int',
        description="position within the file to find this trajectory"
    )

