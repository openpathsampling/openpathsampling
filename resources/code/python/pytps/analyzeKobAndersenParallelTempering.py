#!/usr/local/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Analyze parallel tempering simulations conducted with 'ParallelTempering.py' driver.

DESCRIPTION


COPYRIGHT

@author John D. Chodera <jchodera@gmail.com>

This source file is released under the GNU General Public License.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import os.path
import sys
import math

import numpy

#import Scientific.IO.NetCDF as netcdf # pynetcdf
#import scipy.io.netcdf as netcdf # scipy netcdf
import netCDF4 as netcdf # enthought netCDF4 package

import simtk.unit as units
#import simtk.chem.openmm as openmm

import pymbar
import timeseries

#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: $"

#=============================================================================================
# GLOBAL CONSTANTS
#=============================================================================================


#=============================================================================================
# SUBROUTINES
#=============================================================================================

def write_file(filename, contents):
    """Write the specified contents to a file.

    ARGUMENTS
      filename (string) - the file to be written
      contents (string) - the contents of the file to be written

    """

    outfile = open(filename, 'w')

    if type(contents) == list:
        for line in contents:
            outfile.write(line)
    elif type(contents) == str:
        outfile.write(contents)
    else:
        raise "Type for 'contents' not supported: " + repr(type(contents))

    outfile.close()

    return

def read_file(filename):
    """Read contents of the specified file.

    ARGUMENTS
      filename (string) - the name of the file to be read

    RETURNS
      lines (list of strings) - the contents of the file, split by line
    """

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    return lines

def read_pdb(filename):
    """
    Read the contents of a PDB file.

    ARGUMENTS

    filename (string) - name of the file to be read

    RETURNS

    atoms (list of dict) - atoms[index] is a dict of fields for the ATOM residue

    """
    
    # Read the PDB file into memory.
    pdbfile = open(filename, 'r')

    # Extract the ATOM entries.
    # Format described here: http://bmerc-www.bu.edu/needle-doc/latest/atom-format.html
    atoms = list()
    for line in pdbfile:
        if line[0:6] == "ATOM  ":
            # Parse line into fields.
            atom = dict()
            atom["serial"] = line[6:11]
            atom["atom"] = line[12:16]
            atom["altLoc"] = line[16:17]
            atom["resName"] = line[17:20]
            atom["chainID"] = line[21:22]
            atom["Seqno"] = line[22:26]
            atom["iCode"] = line[26:27]
            atom["x"] = line[30:38]
            atom["y"] = line[38:46]
            atom["z"] = line[46:54]
            atom["occupancy"] = line[54:60]
            atom["tempFactor"] = line[60:66]
            atoms.append(atom)
            
    # Close PDB file.
    pdbfile.close()

    # Return dictionary of present residues.
    return atoms

def write_pdb(filename, trajectory, atoms):
    """Write out replica trajectories as multi-model PDB files.

    ARGUMENTS
       filename (string) - name of PDB file to be written
       trajectory
       atoms (list of dict) - parsed PDB file ATOM entries from read_pdb() - WILL BE CHANGED
    """

    # Create file.
    outfile = open(filename, 'w')

    nframes = trajectory.shape[0]

    # Write trajectory as models
    for frame_index in range(nframes):
        outfile.write("MODEL     %4d\n" % (frame_index+1))

        # Write ATOM records.
        for (index, atom) in enumerate(atoms):
            atom["x"] = "%8.3f" % trajectory[frame_index,index,0]
            atom["y"] = "%8.3f" % trajectory[frame_index,index,1]
            atom["z"] = "%8.3f" % trajectory[frame_index,index,2]
            outfile.write('ATOM  %(serial)5s %(atom)4s%(altLoc)c%(resName)3s %(chainID)c%(Seqno)5s   %(x)8s%(y)8s%(z)8s\n' % atom)

        outfile.write("ENDMDL\n")
        
    # Close file.
    outfile.close()

    return

def show_mixing_statistics(ncfile, show_transition_matrix=False):
    """
    Print summary of mixing statistics.

    """

    print "Computing mixing statistics..."

    states = ncfile.variables['states'][:,:].copy()

    # Determine number of iterations and states.
    [niterations, nstates] = ncfile.variables['states'][:,:].shape
    
    # Compute statistics of transitions.
    Nij = numpy.zeros([nstates,nstates], numpy.float64)
    for iteration in range(niterations-1):
        for ireplica in range(nstates):
            istate = states[iteration,ireplica]
            jstate = states[iteration+1,ireplica]
            Nij[istate,jstate] += 0.5
            Nij[jstate,istate] += 0.5
    Tij = numpy.zeros([nstates,nstates], numpy.float64)
    for istate in range(nstates):
        Tij[istate,:] = Nij[istate,:] / Nij[istate,:].sum()

    if show_transition_matrix:
        # Print observed transition probabilities.
        PRINT_CUTOFF = 0.001 # Cutoff for displaying fraction of accepted swaps.
        print "Cumulative symmetrized state mixing transition matrix:"
        print "%6s" % "",
        for jstate in range(nstates):
            print "%6d" % jstate,
        print ""
        for istate in range(nstates):
            print "%-6d" % istate,
            for jstate in range(nstates):
                P = Tij[istate,jstate]
                if (P >= PRINT_CUTOFF):
                    print "%6.3f" % P,
                else:
                    print "%6s" % "",
            print ""

    # Estimate second eigenvalue and equilibration time.
    mu = numpy.linalg.eigvals(Tij)
    mu = -numpy.sort(-mu) # sort in descending order
    if (mu[1] >= 1):
        print "Perron eigenvalue is unity; Markov chain is decomposable."
    else:
        print "Perron eigenvalue is %9.5f; state equilibration timescale is ~ %.1f iterations" % (mu[1], 1.0 / (1.0 - mu[1]))

    return

def check_positions(ncfile):
    """Make sure no positions have gone 'nan'.

    ARGUMENTS
       ncfile (NetCDF) - NetCDF file object for input file
    """

    # Get current dimensions.
    niterations = ncfile.variables['positions'].shape[0]
    nstates = ncfile.variables['positions'].shape[1]
    natoms = ncfile.variables['positions'].shape[2]

    # Compute torsion angles for each replica
    for iteration in range(niterations):
        for replica in range(nstates):
            # Extract positions
            positions = array(ncfile.variables['positions'][iteration,replica,:,:])
            # Check for nan
            if any(isnan(positions)):
                # Nan found -- raise error
                print "Iteration %d, state %d - nan found in positions." % (iteration, replica)
                # Report coordinates
                for atom_index in range(natoms):
                    print "%16.3f %16.3f %16.3f" % (positions[atom_index,0], positions[atom_index,1], positions[atom_index,2])
                    if any(isnan(positions[atom_index,:])):
                        raise "nan detected in positions"

    return

def estimate_free_energies(ncfile, ndiscard = 0, nuse = None):
    """Estimate free energies of all sampled states.

    ARGUMENTS
       ncfile (NetCDF) - input YANK netcdf file

    OPTIONAL ARGUMENTS
       ndiscard (int) - number of iterations to discard to equilibration
       nuse (int) - maximum number of iterations to use (after discarding)

    TODO: Automatically determine 'ndiscard'.
    """

    # Get current dimensions.
    niterations = ncfile.variables['energies'].shape[0]
    nstates = ncfile.variables['energies'].shape[1]
    natoms = ncfile.variables['energies'].shape[2]

    # Extract energies.
    print "Reading energies..."
    energies = ncfile.variables['energies']
    u_kln_replica = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    for n in range(niterations):
        u_kln_replica[:,:,n] = energies[n,:,:]
    print "Done."

    # Deconvolute replicas
    print "Deconvoluting replicas..."
    u_kln = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = energies[iteration,:,:]
    print "Done."

    # Compute total negative log probability over all iterations.
    u_n = numpy.zeros([niterations], numpy.float64)
    for iteration in range(niterations):
        u_n[iteration] = sum(numpy.diagonal(u_kln[:,:,iteration]))
    #print u_n

    # DEBUG
    outfile = open('u_n.out', 'w')
    for iteration in range(niterations):
        outfile.write("%8d %24.3f\n" % (iteration, u_n[iteration]))
    outfile.close()

    # Discard initial data to equilibration.
    u_kln_replica = u_kln_replica[:,:,ndiscard:]
    u_kln = u_kln[:,:,ndiscard:]
    u_n = u_n[ndiscard:]

    # Truncate to number of specified conforamtions to use
    if (nuse):
        u_kln_replica = u_kln_replica[:,:,0:nuse]
        u_kln = u_kln[:,:,0:nuse]
        u_n = u_n[0:nuse]
    
    # Subsample data to obtain uncorrelated samples
    N_k = numpy.zeros(nstates, numpy.int32)
    indices = timeseries.subsampleCorrelatedData(u_n) # indices of uncorrelated samples
    indices = range(0,u_n.size) # DEBUG - assume samples are uncorrelated
    N = len(indices) # number of uncorrelated samples
    N_k[:] = N      
    u_kln[:,:,0:N] = u_kln[:,:,indices]
    print "number of uncorrelated samples:"
    print N_k
    print ""

    #===================================================================================================
    # Estimate free energy difference with MBAR.
    #===================================================================================================   
   
    # Initialize MBAR (computing free energy estimates, which may take a while)
    print "Computing free energy differences..."
    #mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'self-consistent-iteration') # use slow self-consistent-iteration (the default)
    mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'Newton-Raphson', initialize='BAR', relative_tolerance=1.0e-13) # use faster Newton-Raphson solver

    # Get matrix of dimensionless free energy differences and uncertainty estimate.
    print "Computing covariance matrix..."
    (Deltaf_ij, dDeltaf_ij) = mbar.getFreeEnergyDifferences()
   
#    # Matrix of free energy differences
    print "Deltaf_ij:"
    for i in range(nstates):
        for j in range(nstates):
            print "%8.3f" % Deltaf_ij[i,j],
        print ""        
    
#    print Deltaf_ij
#    # Matrix of uncertainties in free energy difference (expectations standard deviations of the estimator about the true free energy)
    print "dDeltaf_ij:"
    for i in range(nstates):
        for j in range(nstates):
            print "%8.3f" % dDeltaf_ij[i,j],
        print ""        

    # Return free energy differences and an estimate of the covariance.
    return (Deltaf_ij, dDeltaf_ij)

def estimate_enthalpies(ncfile, ndiscard = 0, nuse = None):
    """Estimate enthalpies of all sampled states.

    ARGUMENTS
       ncfile (NetCDF) - input YANK netcdf file

    OPTIONAL ARGUMENTS
       ndiscard (int) - number of iterations to discard to equilibration
       nuse (int) - number of iterations to use (after discarding) 

    TODO: Automatically determine 'ndiscard'.
    TODO: Combine some functions with estimate_free_energies.
    """

    # Get current dimensions.
    niterations = ncfile.variables['energies'].shape[0]
    nstates = ncfile.variables['energies'].shape[1]
    natoms = ncfile.variables['energies'].shape[2]

    # Extract energies.
    print "Reading energies..."
    energies = ncfile.variables['energies']
    u_kln_replica = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    for n in range(niterations):
        u_kln_replica[:,:,n] = energies[n,:,:]
    print "Done."

    # Deconvolute replicas
    print "Deconvoluting replicas..."
    u_kln = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = energies[iteration,:,:]
    print "Done."

    # Compute total negative log probability over all iterations.
    u_n = numpy.zeros([niterations], numpy.float64)
    for iteration in range(niterations):
        u_n[iteration] = sum(numpy.diagonal(u_kln[:,:,iteration]))
    #print u_n

    # DEBUG
    outfile = open('u_n.out', 'w')
    for iteration in range(niterations):
        outfile.write("%8d %24.3f\n" % (iteration, u_n[iteration]))
    outfile.close()

    # Discard initial data to equilibration.
    u_kln_replica = u_kln_replica[:,:,ndiscard:]
    u_kln = u_kln[:,:,ndiscard:]
    u_n = u_n[ndiscard:]
    
    # Truncate to number of specified conformations to use
    if (nuse):
        u_kln_replica = u_kln_replica[:,:,0:nuse]
        u_kln = u_kln[:,:,0:nuse]
        u_n = u_n[0:nuse]

    # Subsample data to obtain uncorrelated samples
    N_k = numpy.zeros(nstates, numpy.int32)
    indices = timeseries.subsampleCorrelatedData(u_n) # indices of uncorrelated samples
    indices = range(0,u_n.size) # DEBUG - assume samples are uncorrelated
    N = len(indices) # number of uncorrelated samples
    N_k[:] = N      
    u_kln[:,:,0:N] = u_kln[:,:,indices]
    print "number of uncorrelated samples:"
    print N_k
    print ""

    # Compute average enthalpies.
    H_k = numpy.zeros([nstates], numpy.float64) # H_i[i] is estimated enthalpy of state i
    dH_k = numpy.zeros([nstates], numpy.float64)
    for k in range(nstates):
        H_k[k] = u_kln[k,k,:].mean()
        dH_k[k] = u_kln[k,k,:].std() / sqrt(N)

    return (H_k, dH_k)

def construct_atom_list(N, NA):
    """Write out replica trajectories as multi-model PDB files.

    ARGUMENTS
       atoms (list of dict) - parsed PDB file ATOM entries from read_pdb() - WILL BE CHANGED
       filename (string) - name of PDB file to be written
       title (string) - the title to give each PDB file
       ncfile (NetCDF) - NetCDF file object for input file       

    """

    NB = N - NA
    atoms = list()

    index = 1
    for n in range(NA):
        atom = dict()
        atom['serial'] = index
        atom['atom'] = ' Ar '
        atom['altLoc'] = ' '
        atom['resName'] = 'Ar '
        atom['chainID'] = ' '
        atom['Seqno'] = '%5d' % index        
        index += 1
        atoms.append(atom)        
    for n in range(NB):
        atom = dict()
        atom['serial'] = index
        atom['atom'] = ' He '
        atom['altLoc'] = ' '
        atom['resName'] = 'He '
        atom['chainID'] = ' '
        atom['Seqno'] = '%5d' % index        
        index += 1
        atoms.append(atom)
        
    return atoms

#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

if __name__ == "__main__":
    netcdf_filename = 'tempering.nc' # name of NetCDF file to analyze

    # Open NetCDF file for reading.
    #ncfile = netcdf.NetCDFFile(netcdf_filename, 'r')
    ncfile = netcdf.Dataset(netcdf_filename, 'r')    

    # Number of particles
    N = 150
    NA = 120

    # Analyze acceptance probabilities.
    show_mixing_statistics(ncfile, show_transition_matrix=True)

    # Estimate free energies
    ndiscard = 0 # number of iterations to discard
    estimate_free_energies(ncfile, ndiscard=ndiscard)
    
    # Get current dimensions.
    niterations = ncfile.variables['states'].shape[0]
    nstates = ncfile.variables['states'].shape[1]

    print "%d iterations, %d states" % (niterations, nstates)
    
    #===================================================================================================    
    # Write states sampled by each replica
    #===================================================================================================
    
    state_nk = ncfile.variables['states'][:,:].copy()
    outfile = open('states.out', 'w')
    for n in range(niterations):
        for k in range(nstates):
            outfile.write("%6d" % state_nk[n,k])
        outfile.write("\n")        
    outfile.close()

    #===================================================================================================    
    # Write trajectories as PDB files.
    #===================================================================================================

    write_trajectories = False # WARNING: THESE FILES ARE LARGE! 
    map_into_box = False    
    if write_trajectories:
        print "Writing trajectories..."    
        atom_list = construct_atom_list(N, NA)
        for replica_index in range(nstates):
            print "replica %d / %d" % (replica_index, nstates)
            filename = "replica-%d.pdb" % replica_index
            trajectory = ncfile.variables['positions'][:,replica_index,:,:].copy() * 10.0
            print "trajectory.shape = %s" % str(trajectory.shape)
            if (map_into_box):
                # map coordinates back into box
                for frame_index in range(nframes):
                    for atom_index in range(natoms):
                        for k in range(3):
                            while (trajectory[frame_index,atom_index,k] < 0.0):
                                trajectory[frame_index,atom_index,k] += (length / units.angstroms)
                            while (trajectory[frame_index,atom_index,k] >= length / units.angstroms):
                                trajectory[frame_index,atom_index,k] -= (length / units.angstroms)
            # Write PDB
            write_pdb(filename, trajectory, atom_list)

        print "Done."

    #===================================================================================================
    # Clean up.
    #===================================================================================================   

    # Close NetCDF file.
    ncfile.close()
        
    print "Done."
    
    
