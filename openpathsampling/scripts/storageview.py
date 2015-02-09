import sys
import argparse
import os
import json


import openpathsampling as paths

from openpathsampling.storage import Storage

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Analyze a file.')
    parser.add_argument('file', metavar='file.nc', help='an integer for the accumulator')

    args = parser.parse_args()
    file = args.file

    if not os.path.isfile(file):
        print file, 'does not exist ! ENDING!'
        exit()

    storage = Storage(
        filename = file,
        mode = 'a'
    )

    class ReadableObjectJSON(paths.todict.ObjectJSON):
        def __init__(self, unit_system = None, class_list = None):
            super(ReadableObjectJSON, self).__init__(unit_system, class_list)
            self.excluded_keys = ['name', 'idx', 'json', 'identifier']
            self.storage = storage

        def build(self,obj):
            if type(obj) is dict:
                if '_cls' in obj and '_idx' in obj:
                    return obj['_cls'] + '[#' + str(obj['_idx']) + ']'

                if '_cls' in obj:
                    return { obj['_cls'] : self.build(obj['_dict'])}

                return { key: self.build(value) for key, value in obj.iteritems() }

            return super(ReadableObjectJSON, self).build(obj)

    def indent(s, width=4):
        spl = s.split('\n')
        spl = [' '*width + p for p in spl]
        return '\n'.join(spl)

    def headline(s):
        print
        print "###############################################################################"
        print "##", s.upper()
        print "###############################################################################"
        print

    def line(a, b):
        print '    {:<32} : {:<30}'.format(a,b)

    def nline(n, a, b):
        print '     {:>4}] {:<25} : {:<30}'.format(n,a,b)

    def format_json(json_str):
        obj = json.loads(json_str)
        return json.dumps(obj, sort_keys=True,
                  indent=2, separators=(',', ': '))

    def format_by_json(obj):
        return json.dumps(obj, sort_keys=True,
                  indent=2, separators=(',', ': '))

    simplifier = ReadableObjectJSON()

    headline("General")

    line("Filename", file)
    line("Size", str(os.path.getsize(file) / 1024 / 1024) + " MB")

    headline("Content")

    line("Number of trajectories", storage.trajectory.count())
    line("Number of configurations", storage.configuration.count())
    line("Number of momenta", storage.momentum.count())

    headline("Topology")

    topology = storage.topology

    line("Number of Atoms", topology.n_atoms)
    line("Number of Dimensions", topology.n_spatial)

    if type(topology) is paths.MDTrajTopology:
        line('MDTraj Topology','')

#    md_topology = topology.md

#    counterion_indices = [ a.index for a in md_topology.atoms if a.residue.name[-1] == '+']
#    solvent_indices = [ a.index for a in md_topology.atoms if a.residue.name == 'HOH']
#    protein_indices = [ a.index for a in md_topology.atoms if a.residue.name[-1] != '+' and a.residue.name != 'HOH']

#    line("Number of waters", len(solvent_indices) / 3)
#    line("Number of protein atoms", len(protein_indices))

    headline("Snapshot Zero")
    # load initial equilibrate snapshot given by ID #0
    snapshot = storage.snapshot.load(0)

    line("Potential Energy",str(snapshot.potential_energy))
    line("Kinetic Energy",str(snapshot.kinetic_energy))

    headline("Ensembles")

    for e_idx in range(0, storage.ensemble.count()):
        ensemble = storage.ensemble.load(e_idx)
        nline(e_idx,ensemble.cls, '')
#        print indent(str(ensemble),16)
        print indent(format_by_json(simplifier.from_json(ensemble.json)), 16)

    headline("PathMovers")

    for p_idx in range(0, storage.pathmover.count()):
        pathmover = storage.pathmover.load(p_idx)
        nline(p_idx,pathmover.name, '')
        print indent(format_by_json(simplifier.from_json(pathmover.json)), 16)

    headline("ShootingPointSelector")

    for p_idx in range(0, storage.shootingpointselector.count()):
        obj = storage.shootingpointselector.load(p_idx)
        nline(p_idx, obj.cls, '')
#        print indent(format_by_json(simplifier.from_json(obj.json)), 16)

    headline("ShootingPoints (" + str(storage.shootingpoint.count()) + ")")

#    for p_idx in range(0, storage.shootingpoint.count()):
#        obj = storage.shootingpoint.load(p_idx)
#        nline(p_idx,obj.json, obj.cls)

    headline("Orderparameters (" + str(storage.collectivevariable.count()) + ")")

    for p_idx in range(0, storage.collectivevariable.count()):
        obj = storage.collectivevariable.load(p_idx)
        add = ''
        if len(obj.storage_caches[storage])>0:
            add = '{ %d : %f, ... } ' % obj.storage_caches[storage].iteritems().next()
        nline(p_idx,obj.name, str(len(obj.storage_caches[storage])) + ' entries ' + add)

    headline("MovePaths")

    for p_idx in range(0, storage.movepaths.count()):
        obj = storage.movepaths.load(p_idx)
        nline(p_idx, '', '')
        print indent(str(obj),16)


    headline("SampleSets")

    for p_idx in range(0, storage.sampleset.count()):
        obj = storage.sampleset.load(p_idx)
        nline(p_idx, str(len(obj.samples)) + ' sample(s)', [storage.idx(sample) for sample in obj.samples ])
        print indent(str(obj.movepath),16)


    headline("Samples")

    def shortened_dict(d):
        keys = sorted(d.keys())
        old_idx = -2
        count = 0
        for idx in keys:
            if idx == old_idx + 1 or idx == old_idx - 1:
                count += 1
            else:
                if count > 1:
                    sys.stdout.write(" <" + str(count - 1) + ">")
                if old_idx >= 0 and count > 0:
                    sys.stdout.write(" " + str(old_idx))
                sys.stdout.write(" " + str(idx))
                count = 0
            old_idx = idx

        if count > 1:
            sys.stdout.write(" <" + str(count - 1) + "> ")
        if count > 0:
            sys.stdout.write(" " + str(old_idx))

    def format_traj(traj_obj):
        s = ''
        traj = storage.trajectory.snapshot_indices(traj_obj.idx[storage])
        old_idx = -2
        count = 0
        for idx in traj:
            if idx == old_idx + 1 or idx == old_idx - 1:
                count += 1
            else:
                if count > 1:
                    s += " <" + str(count - 1) + ">"
                if old_idx >= 0 and count > 0:
                    s += " " + str(old_idx)
                s += " " + str(idx)
                count = 0
            old_idx = idx

        if count > 1:
            s += " <" + str(count - 1) + "> "
        if count > 0:
            s += " " + str(old_idx)

        return s

    def print_traj(name, traj_obj):
        traj = storage.trajectory.snapshot_indices(traj_obj.idx[storage])
        sys.stdout.write("      {:>10}:  {:>5} frames [".format(name, str(len(traj))) + ' ]')
        print format_traj(traj_obj)


    for o_idx in range(0, storage.sample.count()):

        sample = storage.sample.load(o_idx)
        nline(o_idx, 'trajectory #' + str(storage.idx(sample.trajectory)),'[ ' + str(len(sample.trajectory)) + ' frames ] ' + format_traj(sample.trajectory))

    headline("Trajectories")

    for t_idx in range(0, storage.trajectory.count()):
        traj = storage.trajectory.snapshot_indices(t_idx)
        sys.stdout.write("  {:>4} [{:>5} frames] : ".format(str(t_idx),str(len(traj))))
        old_idx = -2
        count = 0
        for idx in traj:
            if idx == old_idx + 1 or idx == old_idx - 1:
                count += 1
            else:
                if count > 1:
                    sys.stdout.write(" <" + str(count - 1) + ">")
                if old_idx >= 0 and count > 0:
                    sys.stdout.write(" " + str(old_idx))
                sys.stdout.write(" " + str(idx))
                count = 0
            old_idx = idx

        if count > 1:
            sys.stdout.write(" ... <" + str(count - 1) + "> ...")
        if count > 0:
            sys.stdout.write(" " + str(old_idx))


        sys.stdout.write("\n")
