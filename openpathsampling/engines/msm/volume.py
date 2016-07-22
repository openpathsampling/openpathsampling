'''
Created on 03.09.2014

@author: Jan-Hendrik Prinz, David W.H. Swenson
'''

from openpathsampling.volume import Volume


class Block(Volume):
    """
    Define a volume in discrete state space
    """

    def __init__(self, states):
        super(Volume, self).__init__()
        if type(states) is int:
            self.states = [states]
        elif type(states) is slice:
            self.states = range(*q.indices(q.stop))
        elif type(states) is bool:
            self.states = states

        self.states = frozenset(states)

    def __call__(self, snapshot):
        if type(self.states) is bool:
            return self.states

        return snapshot.state in self.states
                
    def __str__(self):
        """
        Returns a string representation of the volume
        """
        if len(self.states) == 1:
            return 'State[%s]' % self.name
        elif len(self.states) > 1:
            return 'State[%d, ..., %d]' % min(self.states), max(self.states)
        else:
            return 'State[]'

    def __eq__(self, other):
        if type(other) is Block:
            return self.states == other.states

        return False

    def __and__(self, other):
        if isinstance(self, Block):
            return Block()

    def __or__(self, other):
        if isinstance(self, Block):

    def __sub__(self, other):
        if isinstance(self, Block):

