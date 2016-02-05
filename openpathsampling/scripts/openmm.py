import simtk.openmm as mm
import os
import os.path
import re
import conda.api
import conda.config


c1 = re.compile('([a-zA-Z0-9 ]*) ([\S]*): dlopen\(([\S]*), ([0-9]*)\): Library not loaded: ([\S]*)([\s]*)Referenced from: ([\S]*)([\s]*)Reason: ([a-zA-Z0-9 ]*)')
subdirs = ['DIR', 'PLUGIN_DIR', 'LIB_PATH']
openmm_envs = [env for env in os.environ if 'openmm' in env.lower()]
err = []

if not os.path.isdir(mm.Platform.getDefaultPluginsDirectory()):
    err += ['The default Plugin path "{0}" does not exist.'.format(
        mm.Platform.getDefaultPluginsDirectory()
        )]

# for env in openmm_envs:
#     if 'dir' in env.lower() or 'path' in env.lower():
#         path = os.environ[env]
#         if not os.path.isdir(path):
#             err += ['Error in environment variable "{0}": Path "{1}" does not exist!'.format(env, path)]

for sub in subdirs:
    env = 'OPENMM_' + sub
    if env not in os.environ:
        err += ['Needed environment variable "{0}" not set!'.format(env)]
    else:
        path = os.environ[env]

        if not os.path.isdir(path):
            err += ['Error in environment variable "{0}": Path "{1}" does not exist!'.format(env, path)]
        else:
            files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
            okay = False
            for f in files:
                if 'openmm' in os.path.join(path, f).lower():
                    okay = True

            if not okay:
                err = ['Folder "{0}" in environment variable "{1}" does not contain files that are openmm specific. Check if the path is correct!'.format(path, env)]

for s in mm.Platform.getPluginLoadFailures():
    res = c1.search(s)
    groups = res.groups()
    if not os.path.isfile(groups[1]):
        err += ['While loading a dynamic library the requested file "{0}" could not be found!'.format(groups[1])]
    if not os.path.isfile(groups[2]) and groups[1] != groups[2]:
        err += ['While loading a dynamic library the requested file "{0}" could not be found!'.format(groups[1])]

if not err:
    print 'No apparent errors detected'
else:
    print 'The following errors have been detected'
    for e in err:
        print e


print 'Found the following installed platforms'

possible_platforms = ['Reference', 'CPU', 'OpenCL', 'CUDA']

platforms = []

for pf_idx in range(mm.Platform.getNumPlatforms()):
    pf = mm.Platform.getPlatform(pf_idx)
    platforms += [(pf.getName(), pf.getSpeed(), '{0:20s} : {1:3.0f} x'.format(pf.getName(), pf.getSpeed()))]

for pf in sorted(platforms, key=lambda x : x[1]):
    print pf[2]

if 'CUDA' not in [pf[0] for pf in platforms]:
    print 'CUDA not working. If your machine has an NVidia graphics card make sure you have CUDA SDK Toolkit 7.0 and the CUDA Drivers installed.'

if 'CPU' not in [pf[0] for pf in platforms]:
    print 'CPU not working. This means probably that your Plugin directory is not set correctly. Check the above error messages.'

if err:
    conda_dir = conda.config.default_prefix
    pkgs_dir = conda.config.pkgs_dirs[0]

    for d in conda.config.pkgs_dirs:
        for f in os.listdir(d):
            if f[0] !='_' and f[-4:] != '.txt' and f[-4:] != '.bz2':
                vals = f.split('-')
                if len(vals) >= 3:
                    version = vals[-2]
                    last = vals[-1].split('_')
                    name = '-'.join(vals[0:-2])

                    sub_ver = last[-1]
                    if len(last) > 1:
                        specs = re.findall('([a-zA-z]+[0-9]+)', last[0])

                    path = os.path.join(d,f)

                    if name == 'openmm':
                        txt = """
# PATHS for OPENMM
export OPENMM_DIR="%s"
export OPENMM_PLUGIN_DIR="${OPENMM_DIR}/lib/plugins"
export OPENMM_INCLUDE_PATH="${OPENMM_DIR}/include"
export OPENMM_LIB_PATH="${OPENMM_DIR}/lib"
                        """

                        print "To fix the environment variables you could try adding the following passage to the end of your .bashrc"
                        print txt % path
