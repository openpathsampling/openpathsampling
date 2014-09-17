#!/usr/bin/env python

import sys
import math

# Classes for handling basic input, output, and manipulation of frames.
# Currently supports input and output of xyz-vel data (VMD-compatible XYZ
# files with velocities).

# TODO: fix so that we don't assume that all masses are the same in the
# randomization of velocities

class NoPBC(object):
    
    def __init__(self):
        pass

    def distance_vector(self, pos1, pos2):
        diff = []
        mag2 = 0.0
        for d in range(len(pos1)):
            dx = pos2[d]-pos1[d]
            mag2 += dx*dx
            diff.append(dx)
        return (diff, math.sqrt(mag2))


class Cubic_PBC(object):
    
    def __init__(self):
        self.dx = []

    def distance_vector(self, pos1, pos2):
        from math import sqrt
        diff = []
        mag2 = 0.0
        for d in range(len(self.dx)):
            dx = pos2[d]-pos1[d]
            if (dx > 0.5*self.dx[d]):
                dx -= self.dx[d]
                dx *= -1.0
            elif (dx < -0.5*self.dx[d]):
                dx += self.dx[d]
                dx *= -1.0
            diff.append(dx)
            mag2 += dx*dx
        return (diff, sqrt(mag2))

def matrix_vect_mul(m, v):
    if (len(v)==2 and len(m)==4):
        v0 = m[0]*v[0] + m[1]*v[1]
        v1 = m[2]*v[0] + m[3]*v[1]
        res = [v0, v1]
    elif (len(v)==3 and len(m)==9):
        v0 = m[0]*v[0] + m[1]*v[1] + m[2]*v[2]
        v1 = m[3]*v[0] + m[4]*v[1] + m[5]*v[2]
        v2 = m[6]*v[0] + m[7]*v[1] + m[8]*v[2]
        res = [v0, v1, v2]
    else:
        print "problem with fast mv multiplier: m,v=",m,v
        quit()
    return res

def matrix_inv(m):
    if (len(m)==4):
        rdet = 1.0/(m[0]*m[3]-m[1]*m[2])
        minv = [ m[3]*rdet, -m[1]*rdet, -m[2]*rdet, m[0]*rdet ]
    else:
        print "only have fast inversion for 2D systems"
        quit()
    return minv


class Parallelipiped_PBC(object):
    
    def __init__(self):
        self.sides = []
        self.angles = []

        self.h = []
        self.hinv = []
        self.len_abc = []
        return

    def set_h(self, h):
        self.h = h
        self.hx = []
        if (len(self.h) == 4):
            len0=math.sqrt(self.h[0]*self.h[0]+self.h[2]*self.h[2])
            len1=math.sqrt(self.h[1]*self.h[1]+self.h[3]*self.h[3])
            self.len_abc = [len0, len1]
            dd = 2

        for i in range(len(self.h)):
            self.hx.append(0.0)

        for d in range(dd):
            rr = 0.0
            for i in range(dd):
                rr += self.h[i*dd+d]*self.h[i*dd+d]
            rr = 1.0/math.sqrt(rr)
            for i in range(dd):
                self.hx[i*dd+d] = self.h[i*dd+d]*rr

        self.hxinv = matrix_inv(self.hx)
        return

    def distance_vector(self, pos1, pos2):
        dd = len(self.len_abc)
        dx = []
        for d in range(dd):
            dx.append(pos2[d]-pos1[d])
        dx = matrix_vect_mul(self.hxinv, dx)
        for d in range(dd):
            if (dx[d] > 0.5*self.len_abc[d]):
                dx[d] -= self.len_abc[d]
            elif (dx[d] < -0.5*self.len_abc[d]):
                dx[d] += self.len_abc[d]
        dx = matrix_vect_mul(self.hx, dx)
        dr2 = 0.0
        for d in range(dd):
            dr2 += dx[d]*dx[d]
        dist = math.sqrt(dr2)
        rdr = 1.0/dist
        for d in range(dd):
            dx[d] *= rdr
        return (dx, dist)

    def nearest_replica(pos1, pos2):
        dd = len(self.len_abc)
        replica = []
        pass
            
            

def load_pbc_from_file(fname):
    from FileReader import FileReader
    freader = FileReader(fname, "s")
    arr = freader.next_line_as_array()
    if (len(arr)>0):
        if (arr[0] == "PARALLELIPIPED"):
            h = arr[1:]
            pbc = Parallelipiped_PBC()
            pbc.set_h(h)
        elif (arr[0] == "RECTANGULAR"):
            dx = arr[1:]
            pbc = Cubic_PBC()
            pbc.dx = dx
    else:
        pbc = NoPBC()
    
    return pbc

    pass
    


class TrajFrame(object):

    def __init__(self):
        self.vel = []
        self.pos = []
        self.labels = []
        self.mass = []
        self.natoms = 0
        self.molecule = ""
        self.dd=3
        self.cubic_pbc=1000000

    def total_linear_momentum(self):
        p_lin = []
        for d in range(self.dd):
            p_lin.append(0.0)

        for atom in range(self.natoms):
            for d in range(self.dd):
                p_lin[d] += self.vel[atom][d]*self.mass[atom]

        return p_lin

    def tare_linear_momentum(self):
        p_lin = self.total_linear_momentum()
        inv_natoms = 1.0/self.natoms
        for d in range(self.dd):
            p_lin[d] *= inv_natoms
        
        for atom in range(self.natoms):
            for d in range(self.dd):
                self.vel[atom][d] -= p_lin[d]/self.mass[atom]
        return

    # NOTE: this assumes a cubic box
    # we keep this one around for compatibility with old scripts
    def distance_vector(self, atom1, atom2):
        pos1 = self.pos[atom1]
        pos2 = self.pos[atom2]
        return self.pbc.distance_vector(pos1, pos2)

    def ke(self):
        ke = 0.0;
        for atom in range(self.natoms):
            atom_ke = 0.0
            for d in range(self.dd):
                atom_ke += self.vel[atom][d]**2
            atom_ke *= self.mass[atom]
            ke += atom_ke
        ke *= 0.5
        return ke

    # TODO: move this into an external script
    def random_vel_fixed_E(self, kinetic_energy):
        from random import random
        from math import sqrt
        new_ke = 0.0
        for atom in range(self.natoms):
            for d in range(self.dd):
                self.vel[atom][d] = 2.0*random()-1.0

        self.tare_linear_momentum()
        new_ke = self.ke()
        scaling = sqrt(kinetic_energy/new_ke) # FIXME this assumes ident mass
        for atom in range(self.natoms):
            for d in range(self.dd):
                self.vel[atom][d] *= scaling
        return

    def print_frame_header(self, f=sys.stdout):
        if (len(self.labels)==0): # if we don't have labels, use atom number
            for i in range(self.natoms):
                self.labels.append(str(i+1))

        f.write(str(self.natoms)+"\n")
        f.write(self.molecule)
        return

    def reverse_momenta(self):
        for atom in range(self.natoms):
            for d in range(self.dd):
                self.vel[atom][d] *= -1.0
        return

    def print_frame_xyzvel(self, f=sys.stdout):
        self.print_frame_header(f)
        for i in range(self.natoms):
            f.write(str(self.labels[i]))
            for d in range(self.dd): 
                f.write(" "+str(self.pos[i][d]))
            if (len(self.vel[i]) == self.dd):
                for d in range(self.dd):
                    f.write(" "+str(self.vel[i][d]))
            f.write("\n")
                
        return       


def parse_line_xyz_vel(line):
    from re import sub, split
    splitter = split('\s+', line) # whitespace separated
    while (splitter[0] == ''):
        del splitter[0]
    split_float = [splitter[0]] # label sits here
    # TODO: handle more circumstances
    #   1. only positions
    #   2. positions, velocities, and masses
    while (len(splitter) > 7):  # sometimes end w/ empty items
        del splitter[-1]
    if (len(splitter) < 7):     # same thing if there are no velocities
        while (len(splitter) > 4):
            del splitter[-1]

    if ((len(splitter) == 7) or (len(splitter) == 4)):
        for val in splitter[1:]: #make everything else into floats
            split_float.append(float(val))

    return split_float


class TrajFile(object):
    def __init__(self):
        self.frames = []
        self.f=0
        return

    def read_next_frame(self, fname):
        natoms=0 # stays this way if we hit EOF, right?
        if (self.f == 0):
            self.f = sys.stdin if fname=="" else open(fname, 'r')
        natoms_line = self.f.readline()
        if (natoms_line != ''):
            natoms=int(natoms_line)
        myframe = TrajFrame()
        myframe.natoms = natoms
        nlines = natoms+1
        for i in range(nlines):
            line = self.f.readline()
            if (i==0):
                myframe.molecule = line
            else:
                vals = parse_line_xyz_vel(line)
                myframe.labels.append(vals[0])
                myframe.pos.append(vals[1:4])
                myframe.vel.append(vals[4:7])
        if (natoms==0):
            self.f.close()
        return myframe

    def read_xyz(self, fname):
        self.f = sys.stdin if fname=="" else open(fname, 'r')
        newframe = self.read_next_frame(fname)
        while (newframe.natoms != 0):
            self.frames.append(newframe)
            newframe = self.read_next_frame(fname)
        return


    def write_xyz(self, fname):
        f = sys.stdout if fname=="" else open(fname, 'w')
        for frame in self.frames:
            frame.print_frame_xyzvel(f)
        if (f!=sys.stdout):
            f.close()
        return

    # NOTE: This hasn't actually been tested. But it looks like it should work
    def reverse(self):
        reverse_frames = []
        nframes = len(self.frames)
        for i in range(nframes):
            reverse_frames.append(self.frames[nframes-1-i])
            reverse_frames[i].reverse_momenta()
        self.frames = reverse_frames
        return

def mass_parse(str_mass, natoms):
    mass = []
    if (str_mass == "unit"):
        for i in range(natoms):
            mass.append(1.0)
    else:
        mass = eval(str_mass)
        for i in range(len(mass)):
            mass[i] = float(mass[i])
    return mass

def frames_equal(f1, f2):
    isEqual = True
    if f1.natoms!=f2.natoms: isEqual = False
    i=0
    while (i<range(len(f1.natoms)) and isEqual==True):
        if (f1.mass[i] != f2.mass[i]): isEqual = False
        if (f1.labels[i] != f2.labels[i]): isEqual = False
        if (f1.dd != f2.dd): isEqual = False
        d=0
        while (d<f1.dd and isEqual==True):
            if (f1.pos[i][d] != f2.pos[i][d]): isEqual=False
            if (f1.vel[i][d] != f2.vel[i][d]): isEqual=False
            d+=1
        i+=1

    return isEqual

def trajs_equal(t1, t2):
    isEqual = True
    i=0
    while (i<range(len(t1.frames)) and isEqual==True):
        if (not frames_equal(t1.frames[i], t2.frames[i])):
            isEqual = False
        i += 1
    return isEqual


if __name__ == "__main__":
    myfile = TrajFile()
    if len(sys.argv)==2:
        sys.argv.append("")
    myfile.read_xyz(str(sys.argv[1]))
    myfile.write_xyz(str(sys.argv[2]))

