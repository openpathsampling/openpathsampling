from openpathsampling.netcdfplus import VariableStore
from openpathsampling.pathsimulator import MCStep


class MCStepStore(VariableStore):
    def __init__(self):
        super(MCStepStore, self).__init__(
            MCStep,
            ['change', 'active', 'previous', 'simulation', 'mccycle']
        )

    def _init(self, units=None):
        super(MCStepStore, self)._init()

        # New short-hand definition
        self.init_variable('change', 'obj.pathmovechanges', chunksizes=(1,))
        self.init_variable('active', 'obj.samplesets', chunksizes=(1,))
        self.init_variable('previous', 'obj.samplesets', chunksizes=(1,))
        self.init_variable('simulation', 'obj.pathsimulators', chunksizes=(1,))
        self.init_variable('mccycle', 'int', chunksizes=(1,))