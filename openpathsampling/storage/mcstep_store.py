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
        self.init_variable('change', 'obj.pathmovechanges')
        self.init_variable('active', 'obj.samplesets')
        self.init_variable('previous', 'obj.samplesets')
        self.init_variable('simulation', 'obj.pathsimulators')
        self.init_variable('mccycle', 'int')