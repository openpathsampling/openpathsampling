from openpathsampling.netcdfplus import VariableStore
from openpathsampling.pathsimulators import MCStep


class MCStepStore(VariableStore):
    def __init__(self):
        super(MCStepStore, self).__init__(
            MCStep,
            ['simulation', 'mccycle', 'previous', 'active', 'change']
        )

    def initialize(self, units=None):
        super(MCStepStore, self).initialize()

        # New short-hand definition
        self.create_variable('change', 'obj.movechanges')
        self.create_variable('active', 'obj.samplesets')
        self.create_variable('previous', 'obj.samplesets')
        self.create_variable('simulation', 'obj.pathsimulators')
        self.create_variable('mccycle', 'int')
