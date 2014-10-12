"""
Stubs and other tricks used across many tests to get make things quack like
a duck.

@author David W.H. Swenson
"""

class CallIdentity(object):
    '''Stub for a callable that returns itself'''
    def __init__(self):
        self.name = "Id"
    def __call__(self, value):
        return value


class AtomCounter(object):
    '''Let's be honest: that's all we're using the simulation.system object
    for. So I'll duck-punch.'''
    def __init__(self, natoms):
        self.natoms = natoms

    def getNumParticles(self):
        '''QUAAAAACK'''
        return self.natoms

class SimulationDuckPunch(object):
    '''This is what happens when you find a stranger in the Alps.'''
    def __init__(self, topology, system):
        self.system = system
        self.topology = topology
        


