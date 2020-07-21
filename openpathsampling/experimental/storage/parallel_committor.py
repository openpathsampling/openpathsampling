import openpathsampling as paths


class NewCommittor(paths.PathSimulator):
    def __init__(self, vol, engine, storage):
        super().__init__(storage)
        self.vol = vol
        self.init_ens = paths.LengthEnsemble(1) & paths.AllInXEnsemble(vol)
        self.final_ens = paths.LengthEnsemble(5) & paths.AllInXEnsemble(vol)
        self.engine = engine
        self.mover = paths.ForwardExtendMover(
            ensemble=self.init_ens,
            target_ensemble=self.final_ens,
            engine=engine
        )


    def to_dict(self):
        return {'vol': self.vol,
                'init_ens': self.init_ens,
                'final_ens': self.final_ens,
                'engine': self.engine,
                'mover': self.mover,
                'storage': self.storage}

    @classmethod
    def from_dict(cls, dct):
        obj = cls(dct['vol'], dct['engine'], dct['storage'])
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
        for snap in snaps:
            for step_num in range(n_steps):
                step = self.committor_task(snap, step_num)
