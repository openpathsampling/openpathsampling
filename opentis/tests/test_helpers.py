"""
Stubs and other tricks used across many tests to get make things quack like
a duck.

@author David W.H. Swenson
"""

import os
from pkg_resources import resource_filename

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
        
def prepend_exception_message(e, failmsg):
    """Modifies an exception by prepending failmsg"""
    if not e.args:
        e.args = [failmsg]
    else:
        arg0 = failmsg+e.args[0]
        e.args = tuple([arg0] + list(e.args[1:]))

def data_filename(fname, subdir='test_data'):
    return resource_filename('opentis', 
                             os.path.join('tests', subdir, fname))

