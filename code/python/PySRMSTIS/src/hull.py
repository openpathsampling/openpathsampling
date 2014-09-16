'''
@author David W.H. Swenson
'''
import scipy.spatial as qhull
import numpy as np


class Edge(object):
    '''
    Object to mark edges in a simplicial decomposition. Really just a
    convenience struct: Edge.edge gives the pair of vertex id number that
    make the edge, and Edge.simplices gives a list of simplex numbers which
    share the edge.
    '''
    def __init__(self, edge, simplices):
        ''' Obvious initialization.  '''
        self.edge = edge
        self.simplices = simplices


class Hull(object):
    '''
    Hulls start with a Delaunay triangulation, and define a volume as a set
    of contiguous simplices. 
    '''

    def __init__(self, vertices=[], dist_test=None, initpt=[], autogen=True):
        '''
        Initialization. If dist_test and autopt are defined, and if autogen
        is true, then the initizaliation will also generate the hull
        according to the algorithm in generate_hull().
        '''
        self.vertices = np.array(vertices)
        self.initpt = initpt
        self.dist_test = dist_test
        self.edges  = []  # one for every edge
        self.labels = []  # one for every simplex
        self.regions = { -2: [], 0: []} # default regions
        if (self.dist_test and self.initpt != [] and autogen):
            self.generate_hull()

    def convex_hull(self, points=[]):
        '''
        Calculates the triangulation, and makes sure all the relevant data
        structures are set up: vertices, simplices, edges, neighbors for
        each simplex, etc.
        '''
        if (points==[]): points=self.vertices
        self.hull = qhull.Delaunay(points)
        # abstract out the qhull object as much as possible
        self.find_simplex = self.hull.find_simplex 
        self.simplices = self.hull.simplices
        self.simplex_neighbors = self.hull.neighbors

        # set up the edges list
        edges_ptr = self.hull.vertex_neighbor_vertices[0]
        edges_val = self.hull.vertex_neighbor_vertices[1]
        edges = []
        for i in range(len(self.vertices)):
            for j in edges_val[edges_ptr[i]:edges_ptr[i+1]]:
                if (j>i):
                    edges.append([i, j])
        # set up vertex_simplices and edge_simplices
        vertex_simplices = []
        for i in range(len(self.vertices)):
            vertex_simplices.append([])
            for s in range(len(self.simplices)):
                if i in self.simplices[s]:
                    vertex_simplices[i].append(s)
        edge_simplices = []
        for edge in edges:
            myedgesimplices = []
            for s in vertex_simplices[edge[0]]:
                if s in vertex_simplices[edge[1]]:
                    myedgesimplices.append(s)
            edge_simplices.append(myedgesimplices)
        for i in range(len(edges)):
            self.edges.append(Edge(edges[i], edge_simplices[i]))

        # set up labels
        for i in range(len(self.simplices)):
            self.labels.append(0)
            self.regions[0].append(i)

    def move_simplex_to_region(self, simplex, region):
        '''Reassigns a simplex from its old region to a new region.'''
        old_region = self.labels[simplex]
        self.labels[simplex] = region
        self.regions[old_region].remove(simplex)
        self.regions[region].append(simplex)

    def check_edge(self, edge, verbose=False):
        '''
        For each edge, in checks whether the edge satisfies
        dist_test.test(). If not, all simplices associated with the edge are
        marked as out-of-state.
        '''
        dq = self.vertices[edge.edge[0]] - self.vertices[edge.edge[1]]
        testv = self.dist_test.test(dq)
        if (not testv):
            for simplex in edge.simplices:
                self.move_simplex_to_region(simplex, -2)
        if (verbose):
            passfail = "passes" if (testv) else "fails"
            print "Edge from", self.vertices[edge.edge[0]], "to", \
                    self.vertices[edge.edge[1]], passfail

    def check_edges(self):
        ''' Simple loop over all edges. See Hull.check_edge().  '''
        for edge in self.edges:
            self.check_edge(edge)

    def label_connected(self, initial, newregion=1):
        ''' Labels all regions connected to `initial` as belonging to region
        `newregion`.'''
        oldlabel = self.labels[initial]
        to_scan = [initial]
        if (not (newregion in self.regions.keys())):
            self.regions[newregion] = []
        while (len(to_scan) > 0):
            test = to_scan[0]
            self.move_simplex_to_region(test,newregion)
            for neighbor in self.simplex_neighbors[test]:
                if (neighbor>=0 and self.labels[neighbor] == oldlabel):
                    to_scan.append(neighbor)
            to_scan.pop(0)

    def region_neighbors(self, region):
        ''' Returns a list of regions numbers which are neighbors to
        `region`'''
        initial = self.regions[region][0]
        label = self.labels[initial]
        to_scan = [initial]
        rneighbors = []
        visited = []
        while (len(to_scan) > 0):
            test = to_scan[0]
            visited.append(test)
            for neighbor in self.simplex_neighbors[test]:
                neighbor_label=self.labels[neighbor] if neighbor>=0 else -1
                if (neighbor_label == label):
                    if not (neighbor in visited):
                        to_scan.append(neighbor)
                elif (not (neighbor_label in rneighbors)):
                    rneighbors.append(neighbor_label)
            to_scan.pop(0)
        return rneighbors

    def label_all_regions(self, initial, undefined_regions=[0]):
        '''Relabels regions according to connectivity, for any region in the
        set `undefined_regions`'''
        regionnum=sorted(self.regions.keys())[-1] + 1
        self.label_connected(initial,regionnum)
        for r in undefined_regions:
            while len(self.regions[r]) > 0:
                regionnum += 1
                self.label_connected(self.regions[r][0], regionnum)

    def remove_holes(self):
        '''Assigns holes in the interface (region 1) to the interface.'''
        rset = [r for r in self.regions.keys() if r>1]
        rneighbors = {}
        isHole = {}
        for region in rset:
            rneighbors[region] = self.region_neighbors(region)
            isHole[region] = -1

        hole_groups = []
        for region in rset:
            scanned = []
            to_scan = [region]
            while (isHole[region]==-1):
                test = to_scan[0]
                scanned.append(test)
                r_to_scan = [r for r in rneighbors[test] \
                    if r!=1 and not r in scanned]
                if -1 in r_to_scan:
                    for r in scanned:
                        isHole[r] = 0
                to_scan.extend(r_to_scan)
                to_scan.pop(0)
                if r_to_scan==[]:
                    for r in scanned:
                        isHole[r] = 1
                    hole_groups.append(scanned)
        print isHole
        print hole_groups

        for hole in [r for r in isHole.keys() if isHole[r]==1]:
            self.label_connected(self.regions[hole][0], 1)

    def hull_volume_simplices(self):
        ''' Returns simplices that are within the interface.'''
        return self.regions[1]

    def hull_facets(self):
        '''Returns vertex numbers that form the facets of the outside of
        the hull.'''
        # The trick to this algorithm is in the ordering of the
        # self.simplex_neighbors and self.simplices lists. Assume we have a
        # simplex s defined by vertices [a b c]. It has neighboring simplices
        # [na nb nc]. These lists are constructed such that simplex na is
        # the neighbor "across" from vertex a; i.e., the facet shared by
        # simplex s and simplex na includes all the vertices of simplex s
        # *except* vertex a. The Qhull-based Delaunay triangulation builds
        # structures in this format.
        facets = []
        for s in self.regions[1]:
            for i in range(len(self.simplex_neighbors[s])): 
                if self.simplex_neighbors[s][i] not in self.regions[1]:
                    facets.append([v for v in self.simplices[s] 
                                        if v!=self.simplices[s][i] ])
        # TODO: test for correctness (results look credible, tho)
        return facets

    def vertex_list_to_points(self, vlist):
        ''' Takes a list of vertices (arbitrary shape) and returns the
        spatial points associated with those vertices (same shape with added
        dimension of length d, dimensionality of space)
        '''
        nplist = np.array(vlist)
        pts = []
        ndim = len(self.vertices[0])
        for p in nplist.flat:
            pts.append(self.vertices[p])
        pts = np.array(pts).reshape( nplist.shape + tuple([ndim]) )
        return pts
                

    def generate_hull(self,points=[]):
        ''' Default algorithm to build the interface hull. '''
        if (points==[]): points=self.vertices
        self.convex_hull(points)
        self.check_edges()
        initial = self.hull.find_simplex(self.initpt)
        self.label_all_regions(initial)
        self.label_all_regions(self.regions[-2][0],[-2])
        self.remove_holes()


class MaxDist(object):
    '''
    Example of an appropriate distance test object for use with Hull. This
    one sets a maximum distance for each collective variable.
    '''
    def __init__(self, maxima):
        ''' Initialization. Requires maximum value for each CV as a list. '''
        self.maxima = maxima

    def test(self, val):
        ''' Distance tests must have a method called test().  '''
        res=True
        i=0
        while (res and (i<len(self.maxima))):
            res = (self.maxima[i] >= abs(val[i]))
            i += 1
        return res
