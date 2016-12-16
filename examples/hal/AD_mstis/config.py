import os
from os.path import join
import glob
import re

# file paths

# project_path = join(
#     '/Users', 'jan-hendrikprinz', 'Studium', 'git',
#     'openpathsampling', 'examples', 'hal', 'AD_mstis'
# )
#
# resource_path = project_path

project_path = join(
    '/cbio', 'jclab', 'home', 'prinzj', 'projects', 'ops', 'AD_mstis/')

resource_path = join(project_path, 'data/')
data_path = project_path

# project setting

base_name = os.getenv('OPS_PROJECT_NAME', 'ala_mstis')


project_files = glob.glob(join(data_path, base_name + '*.nc'))

r_pnum = r'\[([0-9]+)\]'
number = 1
for f in project_files:
    try:
        matches = re.findall(r_pnum, f)
        n = int(matches[-1])
        if n >= number:
            number = n + 1
    except:
        pass

number = int(os.getenv('OPS_PROJECT_NUMBER', number))

base_storage_fnc = lambda x: 'ala_mstis[%d]' % x
base_storage_name = base_storage_fnc(number)

bfn = lambda no, part: join(data_path, '%s[%d]_%s.nc' % (base_name, no, part))

pdb_file = join(resource_path, 'AD_initial_frame.pdb')

storage_setup = bfn(number, 'setup')
storage_production = bfn(number, 'production')
storage_resetup = bfn(number + 1, 'setup')
storage_split = bfn(number, 'split')

# platform

platform = 'CUDA'
