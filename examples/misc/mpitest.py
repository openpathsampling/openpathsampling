from mpi4py import MPI

import openpathsampling as paths
import openpathsampling.engines.lammps as eng

script = open('simple-lj.lammps','r').read()
engine = eng.Engine(
    inputs=script,
    options={
        'n_steps_per_frame': 10
    }
)

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
print "Proc %d out of %d procs has" % (me,nprocs)

traj = engine.generate(engine.current_snapshot, [paths.LengthEnsemble(1000).can_append])


MPI.Finalize()
