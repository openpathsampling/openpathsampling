import openpathsampling as paths

from openpathsampling.experimental.storage.dask_integration import \
        SerialScheduler

from openpathsampling.experimental.simstore.tools import none_to_default


class NewCommittor(paths.PathSimulator):
    def __init__(self, states, engine, storage, scheduler=None):
        super().__init__(storage)
        # NOTE: velocity randomization is performed with the Gromacs
        # gen_vel parameter
        self.scheduler = none_to_default(scheduler, SerialScheduler())
        self.states = states
        all_states = paths.join_volumes(states)
        self.init_ens = (paths.LengthEnsemble(1)
                         & paths.AllOutXEnsemble(all_states))
        self.final_ens = paths.SequentialEnsemble([
            paths.AllOutXEnsemble(all_states),
            paths.AllInXEnsemble(all_states) & paths.LengthEnsemble(1)
        ])
        self.engine = engine
        self.mover = paths.ForwardExtendMover(
            ensemble=self.init_ens,
            target_ensemble=self.final_ens,
            engine=engine
        )

    def to_dict(self):
        return {'states': self.states,
                'init_ens': self.init_ens,
                'final_ens': self.final_ens,
                'engine': self.engine,
                'mover': self.mover,
                'storage': self.storage}

    @classmethod
    def from_dict(cls, dct):
        obj = cls(dct['states'], dct['engine'], dct['storage'])
        obj.init_ens = dct['init_ens']
        obj.final_ens = dct['final_ens']
        obj.mover = dct['mover']
        return obj

    def committor_task(self, snap, step_num):
        samples = paths.Sample(trajectory=paths.Trajectory([snap]),
                               ensemble=self.init_ens,
                               replica=0)

        change = self.mover.move_core([samples])
        step = paths.MCStep(
            simulation=self,
            mccycle=step_num,
            active=paths.SampleSet(change.trials[0]),
            change=change
        )
        return step

    def run(self, snaps, n_steps):
        task = self.scheduler.wrap_task(self.committor_task)
        for snap_num, snap in enumerate(snaps):
            for shot_num in range(n_steps):
                step_num = snap_num * len(snaps) + shot_num
                result = task(snap, step_num)
                if self.storage is not None:
                    self.scheduler.store_results(self.storage, result)

        self.scheduler.finalize()

