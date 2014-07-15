'''
Created on 10.07.2014

@author: jan-hendrikprinz
'''

#from __future__ import print_function, absolute_import, division
import numpy as np
from msmbuilder import clustering

from msmbuilder.metrics import RMSD

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
            # Use all generators as inside cells
            self.generators = range(tesselation.size)
        else:
            self.generators = generators

    def inside(self, snapshot):
        self.tesselation.assign(snapshot)