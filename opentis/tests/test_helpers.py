"""
Stubs and other tricks used across many tests to get make things quack like
a duck.

@author David W.H. Swenson
"""

import os
from pkg_resources import resource_filename

from nose.tools import assert_items_equal

def assert_equal_array_array(truth, beauty):
    for (t_atom, b_atom) in zip(truth, beauty):
        assert_items_equal(t_atom, b_atom)

def assert_not_equal_array_array(list_a, list_b):
    exist_diff = False
    for (alpha, beta) in zip(list_a, list_b):
        for (elem_a, elem_b) in zip(alpha, beta):
            if elem_a != elem_b:
                exist_diff = True
    return exist_diff



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

def true_func(value):
    return True
