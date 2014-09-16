"""
@author David W.H. Swenson
"""

from hull import Hull
from volume import Volume

class HullVolume(Volume):
    '''Defines a Volume based on a Hull object.

    Attributes
    ----------
    hull : Hull
        a hull object
    orderparameter : orderparameter
        object that takes a system and return a point in the space within
        which the hull is defined
    '''

    # TODO: if we find ourselves looking for both Voronoi cells and Delaunay
    # simplices, one can be sped up by first doing the other

    def __init__(self, hull, orderparameter):
        super(HullVolume, self).__init__()
        self.hull = hull
        self.orderparameter = orderparameter

    def cell(self, snapshot):
        '''Returns the simplex number for the snapshot
        
        Parameters
        ----------
        snapshot : Snapshot
            the snapshot to test

        Returns
        -------
        int
            the number of the simplex (in the order used in the hull)
        '''
        return self.hull.find_simplex(self.orderparameter(snapshot))

    def __call__(self, snapshot):
        '''
        Returns `True` is the snapshot is within the hull

        Parameters
        ----------
        snapshot : Snapshot
            snapshot to be tested

        Returns
        -------
        bool
            `True` if snapshot is within the hull
        '''
        # with enough simplices, a lookup table is_in_volume[simplex] would
        # be faster here, but this'll do the job
        return (self.cell(snapshot) in self.hull.hull_volume_simplices())
