"""
Created on 03.09.2014

@author: Jan-Hendrik Prinz, David W.H. Swenson
"""

from openpathsampling.volume import Volume


class Block(Volume):
    """
    Define a volume in discrete state space
    """

    def __init__(self, states, inverted=False):
        super(Volume, self).__init__()
        if type(states) is int:
            self.states = [states]
        elif type(states) is slice:
            self.states = range(states.indices(states.stop))
        elif type(states) is bool:
            self.states = states

        self.states = frozenset(states)

        self.inverted = inverted

    def __call__(self, snapshot):
        if type(self.states) is bool:
            return self.states

        return snapshot.state in self.states
                
    def __str__(self):
        """
        Returns a string representation of the volume
        """
        if len(self.states) > 10:
            s = 'State[%d, ..., %d]' % (min(self.states), max(self.states))
        elif len(self.states) > 0:
            s = 'State[%s]' % (','.join(map(str, self.states)))
        else:
            s = 'State[]'

        if self.inverted:
            s = 'not ' + s

        return s

    def __eq__(self, other):
        if type(other) is Block:
            return self.states == other.states

        return False

    def __and__(self, other):
        if isinstance(other, Block):
            if self.inverted and other.inverted:
                return Block(self.states | other.states, True)
            elif self.inverted:
                return Block(other.states ^ self.states)
            elif other.inverted:
                return Block(other.states ^ self.states)
            else:
                return Block(self.states & other.states)

    def __or__(self, other):
        if isinstance(other, Block):
            if self.inverted and other.inverted:
                return Block(self.states & other.states, True)
            elif self.inverted:
                return Block(self.states - other.states, True)
            elif other.inverted:
                return Block(other.states - self.states, True)
            else:
                return Block(self.states | other.states)

    def __sub__(self, other):
            if self.inverted and other.inverted:
                return Block(self.states ^ other.states)
            elif self.inverted:
                return Block(other.states | self.states, True)
            elif other.inverted:
                return Block(other.states & self.states)
            else:
                return Block(self.states - other.states)

    def __xor__(self, other):
            if self.inverted and other.inverted:
                return Block(self.states | other.states)
            elif self.inverted:
                return Block(self.states ^ other.states, True)
            elif other.inverted:
                return Block(other.states ^ self.states, True)
            else:
                return Block(self.states ^ other.states)

    def __neg__(self):
        return Block(self.states, not self.inverted)
