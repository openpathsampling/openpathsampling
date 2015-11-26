import sys
import argparse
import os
import json

import openpathsampling as paths
from openpathsampling.storage import AnalysisStorage

from openpathsampling.netcdfplus import ObjectJSON

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Analyze a file.')
    parser.add_argument('file', metavar='file.nc', help='an integer for the accumulator')

    args = parser.parse_args()
    file_name = args.file

    if not os.path.isfile(file_name):
        print file_name, 'does not exist ! ENDING!'
        exit()

    storage = AnalysisStorage(
        filename=file_name
    )

    class ReadableObjectJSON(ObjectJSON):
        def __init__(self, unit_system=None):
            super(ReadableObjectJSON, self).__init__(unit_system)
            self.excluded_keys = ['name', 'idx', 'json', 'identifier']
            self.storage = storage

        def build(self, obj):
            if type(obj) is dict:
                if '_cls' in obj and '_idx' in obj:
                    return obj['_cls'] + '[#' + str(obj['_idx']) + ']'

                if '_cls' in obj:
                    return {obj['_cls']: self.build(obj['_dict'])}

                return {key: self.build(value) for key, value in obj.iteritems()}

            return super(ReadableObjectJSON, self).build(obj)


    def indent(s, width=4):
        spl = s.split('\n')
        spl = [' ' * width + p for p in spl]
        return '\n'.join(spl)


    def headline(s):
        print
        print "###############################################################################"
        print "##", s.upper()
        print "###############################################################################"
        print


    def line(a, b):
        print '    {:<32} : {:<30}'.format(a, b)


    def nline(n, a, b):
        print '     {:>4}] {:<25} : {:<30}'.format(n, a, b)


    def format_json(json_str):
        obj = json.loads(json_str)
        return json.dumps(obj, sort_keys=True,
                          indent=2, separators=(',', ': '))


    def format_by_json(obj):
        return json.dumps(obj, sort_keys=True,
                          indent=2, separators=(',', ': '))


    simplifier = ReadableObjectJSON()

    headline("General")

    line("Filename", file_name)
    line("Size", str(os.path.getsize(file_name) / 1024 / 1024) + " MB")

    headline("Content")

    line("Number of trajectories", len(storage.trajectories))
    line("Number of snapshots", len(storage.snapshots))
    line("Number of configurations", len(storage.configurations))
    line("Number of momenta", len(storage.momenta))

    headline("Topology")

    topology = storage.topology

    line("Number of Atoms", topology.n_atoms)
    line("Number of Dimensions", topology.n_spatial)

    if type(topology) is paths.MDTrajTopology:
        line('MDTraj Topology', '')

    headline("Snapshot Zero")
    # load initial equilibrate snapshot given by ID #0
    snapshot = storage.snapshots.load(0)

    line("Potential Energy", str(snapshot.potential_energy))
    line("Kinetic Energy", str(snapshot.kinetic_energy))

    headline("Ensembles")

    for e_idx in range(0, len(storage.ensembles)):
        ensemble = storage.ensembles.load(e_idx)
        nline(e_idx, ensemble.cls, '')
        print indent(str(ensemble), 16)
    #        print indent(format_by_json(ensemble.to_dict()), 16)

    headline("PathMovers")

    for p_idx in range(0, len(storage.pathmovers)):
        pathmover = storage.pathmovers.load(p_idx)
        nline(p_idx, pathmover.name, '')
    #        print indent(format_by_json(pathmover.to_dict()), 16)

    headline("ShootingPointSelector")

    for p_idx in range(0, len(storage.shootingpointselectors)):
        obj = storage.shootingpointselectors.load(p_idx)
        nline(p_idx, obj.cls, '')
    #        print indent(format_by_json(simplifier.from_json(obj.json)), 16)

    # headline("ShootingPoints (" + str(len(storage.shootingpoints)) + ")")

    #    for p_idx in range(0, storage.shootingpoints.count()):
    #        obj = storage.shootingpoints.load(p_idx)
    #        nline(p_idx,obj.json, obj.cls)

    headline("CollectiveVariables (" + str(len(storage.collectivevariables)) + ")")

    #    all_snapshot_traj = storage.snapshots.all()

    for p_idx in range(0, len(storage.collectivevariables)):
        obj = storage.collectivevariables.load(p_idx)
#        nline(p_idx, obj.name, '')
        add = ''

        found_values = [ (idx, value) for idx, value in enumerate(obj._store_dict.value_store[:]) if value is not None ]
        if len(found_values) > 0:
           add = '{ %s, ... } ' % ( ', '.join(map(lambda val : '%d : %f' % (val[0], val[1]), found_values[0:5] )))

        nline(p_idx,obj.name, str(len(found_values)) + ' entries ' + add)

    headline("MCSteps")

    for p_idx in range(0, len(storage.steps)):
        obj = storage.steps.load(p_idx)
        nline(p_idx, '', '')
        print indent(str(obj.change), 16)

    headline("SampleSets")

    for p_idx in range(0, len(storage.samplesets)):
        obj = storage.samplesets.load(p_idx)
        nline(p_idx, str(len(obj)) + ' sample(s)', [storage.idx(sample) for sample in obj])
        # print indent(str(obj.movepath),16)

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
        traj = storage.trajectories.snapshot_indices(traj_obj.idx(storage.trajectories))
        old_idx = -2
        count = 0
        for idx in traj:
            if idx / 2 == old_idx / 2 + 1 or idx / 2 == old_idx / 2 - 1:
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
            s += " " + str(old_idx / 2) + ('-' if old_idx % 2 == 0 else '+')

        return s


    def print_traj(name, traj_obj):
        traj = storage.trajectories.snapshot_indices(traj_obj.idx[storage.trajectories])
        sys.stdout.write("      {:>10}:  {:>5} frames [".format(name, str(len(traj))) + ' ]')
        print format_traj(traj_obj)


    for o_idx in range(0, len(storage.samples)):
        sample = storage.samples.load(o_idx)
        nline(o_idx, 'trajectory #' + str(storage.idx(sample.trajectory)),
              '[ ' + str(len(sample.trajectory)) + ' frames ] ' + format_traj(sample.trajectory))

    headline("Trajectories")

    for t_idx in range(0, len(storage.trajectories)):
        traj = storage.trajectories.snapshot_indices(t_idx)
        sys.stdout.write("  {:>4} [{:>5} frames] : ".format(str(t_idx), str(len(traj))))
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
