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
        self.create_variable('change', 'obj.pathmovechanges')
        self.create_variable('active', 'obj.samplesets')
        self.create_variable('previous', 'obj.samplesets')
        self.create_variable('simulation', 'obj.pathsimulators')
        self.create_variable('mccycle', 'int')