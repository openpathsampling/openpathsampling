from os.path import join

# file paths

project_path = join(
    '/'
)

project_path = join(
    '/cbio', 'jclab', 'home', 'prinzj', 'projects', 'ops', 'AD_mstis/')

resource_path = join(project_path, 'data/')
data_path = project_path

pdb_file = join(resource_path, 'AD_initial_frame.pdb')

storage_setup = join(data_path, 'ala_mstis_setup.nc')
storage_production = join(data_path, 'ala_mstis_production.nc')

# platform

platform = 'CPU'
