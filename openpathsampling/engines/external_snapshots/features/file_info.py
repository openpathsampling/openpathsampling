"""
Attributes
----------
file_name : str
file_position : int
"""

variables = ['file_name', 'file_position']

schema_entries = [('file_name', 'str'), ('file_position', 'int')]

def netcdfplus_init(store):
    store.create_variable(
        'file_name',
        'str',
        description="the file name that this snapshot is stored in"
    )
    store.create_variable(
        'file_position',
        'int',
        description="position within the file to find this trajectory"
    )

