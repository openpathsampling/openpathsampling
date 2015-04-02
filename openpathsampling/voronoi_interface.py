'''
Created on 10.07.2014

@author: jan-hendrikprinz
'''


class Object(object):
            pass

class VoronoiInterface(object):
    '''
    Allows to define an interface using Voronoi cells with given metric and a cutoff
    '''


    def __init__(self, tesselation, cutoff = 0.005, generators = None):
        '''
        Constructor
        '''

        self.tesselation = tesselation
        self.cutoff = cutoff
        if generators is None:
            # Use all generators as __call__ cells
            self.generators = range(tesselation.size)
        else:
            self.generators = generators

    def __call__(self, snapshot):
        self.tesselation.assign(snapshot)


@restores_as_full_object
class VoronoiVolume(Volume):
    '''
    Volume given by a Voronoi cell specified by a set of centers

    Parameters
    ----------
    orderparameter : OP_Multi_RMSD
        must be an OP_Multi_RMSD orderparameter that returns several RMSDs
    state : int
        the index of the center for the chosen voronoi cell

    Attributes
    ----------
    orderparameter : orderparameter
        the orderparameter object
    state : int
        the index of the center for the chosen voronoi cell

    '''

    def __init__(self, orderparameter, state):
        super(VoronoiVolume, self).__init__()
        self.orderparameter = orderparameter
        self.state = state

    def cell(self, snapshot):
        '''
        Returns the index of the voronoicell snapshot is in

        Parameters
        ----------
        snapshot : Snapshot
            the snapshot to be tested

        Returns
        -------
        int
            index of the voronoi cell
        '''
        distances = self.orderparameter(snapshot)
        min_val = 1000000000.0
        min_idx = -1
        for idx, d in enumerate(distances):
            if d < min_val:
                min_val = d
                min_idx = idx

        return min_idx

    def __call__(self, snapshot, state = None):
        '''
        Returns `True` if snapshot belongs to voronoi cell in state

        Parameters
        ----------
        snapshot : Snapshot
            snapshot to be tested
        state : int or None
            index of the cell to be tested. If `None` (Default) then the internal self.state is used

        Returns
        -------
        bool
            returns `True` is snapshot is on the specified voronoi cell

        '''

        # short but slower would be

        if state is None:
            state = self.state

        return self.cell(snapshot) == state