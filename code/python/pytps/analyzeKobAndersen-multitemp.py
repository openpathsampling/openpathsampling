
#!/usr/local/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Analyze replica-exchange transition path sampling in the s-field of a Kob-Andersen system.

DESCRIPTION


REFERENCES

[1] Hedges LO, Jack RL, Garrahan JP, and Chandler D. Dynamic order-disorder in atomic models
of structural glass-formers. Science 323:1309, 2009.

[2] Minh DDL and Chodera JD. Optimal estimators and asymptotic variances for nonequilibrium
path-ensemble averages. JCP 131:134110, 2009.

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
import time

import numpy

# import Scientific.IO.NetCDF as netcdf # from Scientific
import netCDF4 as netcdf # from netCDF4 available in enthought python

import pymbar # requires pymbar: http://simtk.org/home/pymbar
import timeseries
import simtk.unit as units

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

def analyze_acceptance_probabilities(ncfile, cutoff = 0.0):
    """Analyze acceptance probabilities.

    ARGUMENTS
       ncfile (NetCDF) - NetCDF file to be analyzed.

    OPTIONAL ARGUMENTS
       cutoff (float) - cutoff for showing acceptance probabilities as blank (default: 0.1)
    """

    # Get current dimensions.
    niterations = ncfile.variables['mixing'].shape[0]
    nstates = ncfile.variables['mixing'].shape[1]

    # Compute mean.
    Pij = numpy.mean(ncfile.variables['mixing'][:,:,:], 0)

    # Write title.
    print "Average state-to-state acceptance probabilities"
    print "(Probabilities less than %(cutoff)f shown as blank.)" % vars()
    print ""

    # Write header.
    print "%4s" % "",
    for j in range(nstates):
        print "%6d" % j,
    print ""

    # Write rows.
    for i in range(nstates):
        print "%4d" % i, 
        for j in range(nstates):
            if Pij[i,j] > cutoff:
                print "%6.3f" % Pij[i,j],
            else:
                print "%6s" % "",
            
        print ""

    return

def extract_activities(ncfile):
    """
    Extract activities from replica-exchange netcdf file and write to text file for plotting.mat

    ARGUMENTS
       ncfile (NetCDF) - input replica-exchange netcdf file
    """

    # Get current dimensions.
    (niterations, nstates) = ncfile.variables['energies'].shape

    # Extract energies.
    print "Reading energies..."
    energies = ncfile.variables['energies']
    u_kln_replica = zeros([nstates, nstates, niterations], float64)
    for n in range(niterations):
        u_kln_replica[:,:,n] = energies[n,:,:]
    print "Done."

    # Deconvolute replicas
    print "Deconvoluting replicas..."
    u_kln = zeros([nstates, nstates, niterations], float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = energies[iteration,:,:]
    print "Done."

    # Show all self-energies
    print 'all self-energies for all replicas'
    for iteration in range(niterations):
        for replica in range(nstates):
            state = int(ncfile.variables['states'][iteration,replica])
            print '%12.1f' % energies[iteration, replica, state],
        print ''

    # If no energies are 'nan', we're clean.
    if not any(isnan(energies[:,:,:])):
        return

    # There are some energies that are 'nan', so check if the first iteration has nans in their *own* energies:
    u_k = diag(energies[0,:,:])
    if any(isnan(u_k)):
        print "First iteration has exploded replicas.  Check to make sure structures are minimized before dynamics"
        print "Energies for all replicas after equilibration:"
        print u_k
        sys.exit(1)

    # There are some energies that are 'nan' past the first iteration.  Find the first instances for each replica and write PDB files.
    first_nan_k = zeros([nstates], int32)
    for iteration in range(niterations):
        for k in range(nstates):
            if isnan(energies[iteration,k,k]) and first_nan_k[k]==0:
                first_nan_k[k] = iteration
    if not all(first_nan_k == 0):
        print "Some replicas exploded during the simulation."
        print "Iterations where explosions were detected for each replica:"
        print first_nan_k
        print "Writing PDB files immediately before explosions were detected..."
        for replica in range(nstates):            
            if (first_nan_k[replica] > 0):
                state = ncfile.variables['states'][iteration,replica]
                iteration = first_nan_k[replica] - 1
                filename = 'replica-%d-before-explosion.pdb' % replica
                title = 'replica %d state %d iteration %d' % (replica, state, iteration)
                write_pdb(atoms, filename, iteration, replica, title, ncfile)
                filename = 'replica-%d-before-explosion.crd' % replica                
                write_crd(filename, iteration, replica, title, ncfile)
        sys.exit(1)

    # There are some energies that are 'nan', but these are energies at foreign lambdas.  We'll just have to be careful with MBAR.
    # Raise a warning.
    print "WARNING: Some energies at foreign lambdas are 'nan'.  This is recoverable."
        
    return

def check_positions(ncfile):
    """Make sure no positions have gone 'nan'.

    ARGUMENTS
       ncfile (NetCDF) - NetCDF file object for input file
    """

    # Get current dimensions.
    (niterations, nstates, natoms) = ncfile.variables['positions'].shape

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
    """Estimate free energies of all alchemical states.

    ARGUMENTS
       ncfile (NetCDF) - input YANK netcdf file

    OPTIONAL ARGUMENTS
       ndiscard (int) - number of iterations to discard to equilibration
       nuse (int) - maximum number of iterations to use (after discarding)

    TODO: Automatically determine 'ndiscard'.
    """

    # Get current dimensions.
    (niterations, nstates, natoms) = ncfile.variables['energies'].shape

    # Extract energies.
    print "Reading energies..."
    energies = ncfile.variables['energies']
    u_kln_replica = zeros([nstates, nstates, niterations], float64)
    for n in range(niterations):
        u_kln_replica[:,:,n] = energies[n,:,:]
    print "Done."

    # Deconvolute replicas
    print "Deconvoluting replicas..."
    u_kln = zeros([nstates, nstates, niterations], float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = energies[iteration,:,:]
    print "Done."

    # Compute total negative log probability over all iterations.
    u_n = zeros([niterations], float64)
    for iteration in range(niterations):
        u_n[iteration] = sum(diagonal(u_kln[:,:,iteration]))
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
    N_k = zeros(nstates, int32)
    indices = timeseries.subsampleCorrelatedData(u_n) # indices of uncorrelated samples
    #indices = range(0,u_n.size) # DEBUG - assume samples are uncorrelated
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
    mbar = pymbar.MBAR(u_kln, N_k, verbose = False, method = 'self-consistent-iteration') # use slow self-consistent-iteration (the default)
    #mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'Newton-Raphson') # use faster Newton-Raphson solver

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
    """Estimate enthalpies of all alchemical states.

    ARGUMENTS
       ncfile (NetCDF) - input YANK netcdf file

    OPTIONAL ARGUMENTS
       ndiscard (int) - number of iterations to discard to equilibration
       nuse (int) - number of iterations to use (after discarding) 

    TODO: Automatically determine 'ndiscard'.
    TODO: Combine some functions with estimate_free_energies.
    """

    # Get current dimensions.
    (niterations, nstates, natoms) = ncfile.variables['positions'].shape

    # Extract energies.
    print "Reading energies..."
    energies = ncfile.variables['energies']
    u_kln_replica = zeros([nstates, nstates, niterations], float64)
    for n in range(niterations):
        u_kln_replica[:,:,n] = energies[n,:,:]
    print "Done."

    # Deconvolute replicas
    print "Deconvoluting replicas..."
    u_kln = zeros([nstates, nstates, niterations], float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = energies[iteration,:,:]
    print "Done."

    # Compute total negative log probability over all iterations.
    u_n = zeros([niterations], float64)
    for iteration in range(niterations):
        u_n[iteration] = sum(diagonal(u_kln[:,:,iteration]))
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
    N_k = zeros(nstates, int32)
    indices = timeseries.subsampleCorrelatedData(u_n) # indices of uncorrelated samples
    #indices = range(0,u_n.size) # DEBUG - assume samples are uncorrelated
    N = len(indices) # number of uncorrelated samples
    N_k[:] = N      
    u_kln[:,:,0:N] = u_kln[:,:,indices]
    print "number of uncorrelated samples:"
    print N_k
    print ""

    # Compute average enthalpies.
    H_k = zeros([nstates], float64) # H_i[i] is estimated enthalpy of state i
    dH_k = zeros([nstates], float64)
    for k in range(nstates):
        H_k[k] = u_kln[k,k,:].mean()
        dH_k[k] = u_kln[k,k,:].std() / sqrt(N)

    return (H_k, dH_k)

def computeKasFunctionOfS(ncfile, mbar, slice_start=None, slice_end=None):
    """
    ARGUMENTS

    ncfile (NetCDF file handle) - the netcdf file handle
    mbar (pymbar) - initialized MBAR

    """

    #===================================================================================================
    # Extract data from netcdf file.
    #===================================================================================================
    print "extracting data"
    # Get current dimensions.
    (niterations, nstates) = ncfile.variables['states'].shape
    (nstates, nframes, natoms, ndim) = ncfile.variables['trajectory_coordinates'].shape

    # Get factor weighting sK contribution to log probability.
    # sKfactor = N * t_obs / delta_t
    sKfactor = getattr(ncfile, 'sKfactor')

    # Get s-values for states.
    s_k = ncfile.variables['fields'][:].copy()

    # Get activities (stored by replica).
    activities = ncfile.variables['activities']
    K_kn_replica = numpy.zeros([nstates,niterations], numpy.float64)
    for n in range(niterations):
        K_kn_replica[:,n] = activities[n,:]

    # Compute reduced potentials.
    u_kln_replica = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    log_probabilities = ncfile.variables['log_probabilities']
    for n in range(niterations):
        u_kln_replica[:,:,n] = - log_probabilities[n,:,:]

    print "Deconvolute Replica...."
    # Deconvolute replicas
    u_kln = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    K_kn = numpy.zeros([nstates, niterations], numpy.float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = u_kln_replica[:,:,iteration]
        K_kn[state_indices,iteration] = K_kn_replica[:,iteration]

    # Compute total negative log probability over all iterations.
    u_n = numpy.zeros([niterations], numpy.float64)
    for iteration in range(niterations):
        u_n[iteration] = sum(numpy.diagonal(u_kln[:,:,iteration]))

    #===================================================================================================
    # Discard data not in our desired slice.
    #===================================================================================================

    if slice_start is None:
        slice_start = 0
    if slice_end is None:
        slice_end = niterations

    u_kln_replica = u_kln_replica[:,:,slice_start:slice_end]
    u_kln = u_kln[:,:,slice_start:slice_end]
    u_n = u_n[slice_start:slice_end]
    K_kn = K_kn[:,slice_start:slice_end]

    #===================================================================================================
    # Subsample data to obtain uncorrelated samples
    #===================================================================================================

    print "Subsample..."
    N_k = numpy.zeros(nstates, numpy.int32)
    indices = timeseries.subsampleCorrelatedData(u_n) # indices of uncorrelated samples
    N = len(indices) # number of uncorrelated samples
    N_k[:] = N
    print N_k
    u_kln = u_kln[:,:,indices]
    K_kn = K_kn[:,indices]

    #===================================================================================================
    # Compute K as a function of s.
    #===================================================================================================

    # inital setup
    svalues = numpy.linspace(-0.01, 0.08, 100) # values of reduced s at which to evaluate reduced <K>_s
    Kvalues = 0.0 * svalues # corresponding estimates of reduced <K>_s
    dKvalues = 0.0 * svalues # uncertainty in estimates of reduced <K>_s
    print "Initialize.."
    # Initialize MBAR for slice data.
    mbar_slice = pymbar.MBAR(u_kln, N_k, verbose=False, method='self-consistent-iteration', relative_tolerance=1.0e-15, initial_f_k=mbar.f_k)
    print "Compute K"
    # Compute s values
    for i in range(svalues.size):
        s = svalues[i]
        # Compute effective energy (notice weight is exp[-sK], so contribution to reduced potential is +sK)
        u_kn = u_kln[:,0,:] + s * K_kn * sKfactor
        [K, dK] = mbar_slice.computePerturbedExpectation(u_kn, K_kn)
        Kvalues[i] = K
        dKvalues[i] = dK

    return [svalues, Kvalues, dKvalues]


def computeObservables(mbar, u_kln, K_kn):
    """
    Arguments are:

    mbar -> instance of mbar
    u_kln -> the reduced potential
    K_kn -> the list of reduced activities

    Observables valculated:

    s* star is calculated by estimatin d<K>/ds and then K* star is found to serve as a cut off


    """

    print "Now in function computeObservables "
     #===================================================================================================
    # Compute dK/ds.
    #===================================================================================================

    print "Computing d<K>/ds..."

    svalues = numpy.linspace(-0.01, 0.08, 100) # values of reduced s at which to evaluate reduced d<K>/ds
    chi_values = 0.0 * svalues # corresponding values of reduced - d<K>/ds
    # Arrange samples of K in n-indexing.
    K_nk = numpy.array(numpy.matrix(K_kn).T)
    Kvalues = 0.0 * svalues # corresponding estimates of reduced <K>_s
    dKvalues = 0.0 * svalues # uncertainty in estimates of reduced <K>_s
    K_n = K_kn[mbar.indices]
    for i in range(svalues.size):
        s = svalues[i]
        # Compute effective energy (notice weight is exp[-sK], so contribution to reduced potential is +sK)
        u_kn = u_kln[:,0,:] + s * K_kn[:,:] * sKfactor
        #compute average K at a certain value of s
        [K, dK] = mbar.computePerturbedExpectation(u_kn, K_kn)
        Kvalues[i] = K
        dKvalues[i] = dK
        # Compute weights w_n for this value of s.
        log_w_kn = mbar._computeUnnormalizedLogWeights(u_kn)
        log_w_n = log_w_kn[mbar.indices] # unnormalized log weights, single n indexing
        fs = - pymbar.logsum(log_w_n) # reduced free energy for thermodynamic state at field s
        w_n = numpy.exp(log_w_n + fs)
        # Compute derivative dw/ds
        dw_n = w_n[:] * (-K_n[:]  + (w_n[:] * K_n[:]).sum())
        # Compute d<K>/ds
        chi_values[i] = - (dw_n[:] * K_n[:]).sum()
    maxChiValue = max(chi_values)
    
    #Writing some information to a file
    outfile = open('chi.out', 'w')
    for i in range(svalues.size):
        outfile.write("%12.6f %12.6f\n" % (svalues[i], chi_values[i]))
    outfile.close()

    k_star=0.0
    # now find the corresponding average K value for it which, would then be K* and serves as a cut off
    for i in range (chi_values.size):
        print svalues[i]
        if chi_values[i]==maxChiValue:
            s_star= svalues[i]
            k_star = Kvalues[i]
            print Kvalues[i]
            break
    print "This is the cutoff K*" + str(k_star)
    print "S* is equal to " +str(s_star)

    #===================================================================================================
    # Compute trajectory entropy difference between inactive and active phases at s = 0.
    #===================================================================================================

    compute_trajectory_entropy = True
    if compute_trajectory_entropy:
        print "Computing trajectory entropy difference between inactive and active phases at s = 0..."

        print "minimum K = %f, maximum K = %f" % (K_kn.min(), K_kn.max())

        Kcutoff = k_star # below cutoff is inactive phase; active phase is above cutoff
        LOG_INFINITY = 1000.0 # substitute for log of infinity

        # Compute 'augmented' u_kln from U.
        u_kln_augmented = numpy.zeros([nstates+2, nstates+2, N], numpy.float64)
        u_kln_augmented[0:nstates,0:nstates,0:N] = u_kln[0:nstates,0:nstates,0:N]
        for k in range(nstates):
            inactive_mask = K_kn[k,:] < Kcutoff
            active_mask = K_kn[k,:] >= Kcutoff

            print "k = %5d : %5d inactive, %5d active" % (k, sum(inactive_mask), sum(active_mask))

            u_kln_augmented[k,nstates+0,:] = u_kln[k,0,:]
            u_kln_augmented[k,nstates+1,:] = u_kln[k,0,:]

            u_kln_augmented[k,nstates+0,inactive_mask] += LOG_INFINITY
            u_kln_augmented[k,nstates+1,active_mask] += LOG_INFINITY

        N_k_augmented = numpy.zeros([nstates+2], numpy.int32)
        N_k_augmented[0:nstates] = N

        # Construct estimate of state free energies.
        initial_f_k = numpy.zeros([nstates+2], numpy.float64)
        initial_f_k[0:nstates] = mbar.f_k[0:nstates]
        mbar_augmented = pymbar.MBAR(u_kln_augmented, N_k_augmented, verbose=False, method='self-consistent-iteration', relative_tolerance=1.0e-15, initial_f_k=initial_f_k)
        #mbar = pymbar.MBAR(u_kln_augmented, N_k_augmented, verbose = True, method = 'Newton-Raphson', initialize='BAR', relative_tolerance=1.0e-13)
        [Delta_f_ij, dDelta_f_ij, Delta_u_ij, dDelta_u_ij, Delta_s_ij, dDelta_s_ij] = mbar_augmented.computeEntropyAndEnthalpy(uncertainty_method='svd-ew')

        print "Delta f = %.3f +- %.3f" % (Delta_f_ij[nstates+0,nstates+1], dDelta_f_ij[nstates+0,nstates+1])
        print "Delta u = %.3f +- %.3f" % (Delta_u_ij[nstates+0,nstates+1], dDelta_u_ij[nstates+0,nstates+1])
        print "Delta s = %.3f +- %.3f" % (Delta_s_ij[nstates+0,nstates+1], dDelta_s_ij[nstates+0,nstates+1])
        print "Now printing in function Compute observables "
        print "printing %5d %5d " %(N, nstates)

    #===================================================================================================
    # Compute trajectory entropy difference between inactive and active phases at s = s*.
    #===================================================================================================

    print "Computing trajectory entropy differences between inactive and active phases at s=s*"
    Kcutoff = k_star
    LOG_INFINITY = 1000.0

    computeTrajEntropy = True
    if computeTrajEntropy:

        print "Computing trajectory entropy difference between inactive and active phases at s = 0..."
       
        # Compute 'augmented' u_kln from U.
        u_kln_augmented_star = numpy.zeros([nstates+2, nstates+2, N], numpy.float64)
        u_kln_augmented_star[0:nstates,0:nstates,0:N] = u_kln[0:nstates, 0:nstates, 0:N]
        for k in range(nstates):
            inactive_mask = K_kn[k,:] < Kcutoff
            active_mask = K_kn[k,:] >= Kcutoff

            u_kln_augmented_star[k,nstates+0,:] = u_kln[k,0,:] + s_star * K_kn[k,:] * sKfactor
            u_kln_augmented_star[k,nstates+1,:] = u_kln[k,0,:] + s_star * K_kn[k,:] * sKfactor

            u_kln_augmented_star[k,nstates+0,inactive_mask] += LOG_INFINITY
            u_kln_augmented_star[k,nstates+1,active_mask] += LOG_INFINITY

        N_k_augmented_star = numpy.zeros([nstates+2], numpy.int32)
        N_k_augmented_star[0:nstates] = N

        # Construct estimate of state free energies.
        initial_f_k = numpy.zeros([nstates+2], numpy.float64)
        initial_f_k[0:nstates] = mbar.f_k[0:nstates]
        mbar_augmented_star = pymbar.MBAR(u_kln_augmented_star, N_k_augmented_star, verbose = True, method = 'self-consistent-iteration', relative_tolerance=1.0e-15, initial_f_k=initial_f_k)
        #mbar_augmented = pymbar.MBAR(u_kln_augmented, N_k_augmented, verbose=False, method='self-consistent-iteration', relative_tolerance=1.0e-15)
        [Delta_f_ij_star, dDelta_f_ij_star, Delta_u_ij_star, dDelta_u_ij_star, Delta_s_ij_star, dDelta_s_ij_star] = mbar_augmented_star.computeEntropyAndEnthalpy(uncertainty_method='svd-ew')

        print "Delta f = %.3f +- %.3f" % (Delta_f_ij_star[nstates+0,nstates+1], dDelta_f_ij_star[nstates+0,nstates+1])
        print "Delta u = %.3f +- %.3f" % (Delta_u_ij_star[nstates+0,nstates+1], dDelta_u_ij_star[nstates+0,nstates+1])
        print "Delta s = %.3f +- %.3f" % (Delta_s_ij_star[nstates+0,nstates+1], dDelta_s_ij_star[nstates+0,nstates+1])


 #===================================================================================================
    # Test entropy calculation
    #===================================================================================================
    testingEntropyCalculations = False
    if testingEntropyCalculations:

        print "Testing entropy calculation..."
        [Delta_f_ij, dDelta_f_ij, Delta_u_ij, dDelta_u_ij, Delta_s_ij, dDelta_s_ij] = mbar.computeEntropyAndEnthalpy(uncertainty_method='svd-ew')

        print "Delta_s_ij"
        print Delta_s_ij
        print dDelta_s_ij

        print "Delta f_16 = %8.1f +- %8.1f kT" % (Delta_f_ij[0,nstates-1], dDelta_f_ij[0,nstates-1])
        print "Delta u_16 = %8.1f +- %8.1f kT" % (Delta_u_ij[0,nstates-1], dDelta_u_ij[0,nstates-1])
        print "Delta s_16 = %8.1f +- %8.1f kT" % (Delta_s_ij[0,nstates-1], dDelta_s_ij[0,nstates-1])

    #===================================================================================================
    # Compute PMF in K at s = 0
    #===================================================================================================

    print "Computing PMF in K at s = 0..."

    nbins = 100
    sigma = 0.002
    Kvalues = numpy.linspace(K_n.min() - 0.01, K_n.max() + 0.01, nbins) # values of reduced K at which to estimate PMF
   
    # Compute 'perturbed' u_kln from U.
    u_kln_perturbed = numpy.zeros([nstates, nbins, N], numpy.float64)

    for k in range(nstates):
        for l in range(nbins):
            u_kln_perturbed[k,l,:] = u_kln[k,0,:] + (K_kn[k,:] - Kvalues[l])**2 / (2.0 * sigma**2)

    [delta_fij, ddelta_fij] = mbar.computePerturbedFreeEnergies(u_kln_perturbed)

    # DEBUG
    outfile = open('pmf-K.out', 'w')
    index = delta_fij[0,:].argmin()
    for i in range(Kvalues.size):
        outfile.write("%12.6f %12.6f %12.6f\n" % (Kvalues[i], delta_fij[index,i], ddelta_fij[index,i]))
    outfile.close()






#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

if __name__ == "__main__":
    netcdf_filename = 'repex-Ts.nc' # name of NetCDF file to analyze

    #ncfile = netcdf.NetCDFFile(netcdf_filename, 'r') # for Scientific.IO.NetCDF
    ncfile = netcdf.Dataset(netcdf_filename, 'r') # for netCDF4

    analyze_acceptance_probabilities(ncfile)
        
    # Get current dimensions.
    (niterations, nstates) = ncfile.variables['states'].shape
    (nstates, nframes, natoms, ndim) = ncfile.variables['trajectory_coordinates'].shape

    # Get factor weighting sK contribution to log probability.
    # sKfactor = N * t_obs / delta_t
    sKfactor = getattr(ncfile, 'sKfactor')

    # Get s-values for states.
    s_k = ncfile.variables['fields'][:].copy()
    
    # Get activities (stored by replica).
    activities = ncfile.variables['activities']
    K_kn_replica = numpy.zeros([nstates,niterations], numpy.float64)
    for n in range(niterations):
        K_kn_replica[:,n] = activities[n,:]
    
    # Compute reduced potentials.
    print "Computing reduced potentials..."
    u_kln_replica = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    log_probabilities = ncfile.variables['log_probabilities']    
    for n in range(niterations):
        u_kln_replica[:,:,n] = - log_probabilities[n,:,:]
    print "Done."

    # Deconvolute replicas
    print "Deconvoluting replicas..."
    u_kln = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    K_kn = numpy.zeros([nstates, niterations], numpy.float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = u_kln_replica[:,:,iteration]
        K_kn[state_indices,iteration] = K_kn_replica[:,iteration]
    print "Done."

    # DEBUG
    print u_kln[0,:,0]
    print u_kln[0,0,0] + s_k[:] * K_kn[0,0] * sKfactor
    #raise Exception("Stop.")

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
    ndiscard = 1000
    u_kln_replica = u_kln_replica[:,:,ndiscard:]
    u_kln = u_kln[:,:,ndiscard:]
    u_n = u_n[ndiscard:]
    K_kn = K_kn[:,ndiscard:]
    
    # Truncate to number of specified conformations to use
    nuse = None
    if (nuse is not None):
        u_kln_replica = u_kln_replica[:,:,0:nuse]
        u_kln = u_kln[:,:,0:nuse]
        u_n = u_n[0:nuse]
        K_kn = K_kn[:,0:nuse]

    #===================================================================================================
    # Write sampled K for each replica and state.
    #===================================================================================================    

    outfile = open('K-replica.out', 'w')
    for n in range(K_kn.shape[1]):
        for k in range(nstates):
            outfile.write("%12.6f" % K_kn_replica[k,n])
        outfile.write("\n")        
    outfile.close()            

    outfile = open('K-state.out', 'w')
    for n in range(K_kn.shape[1]):
        for k in range(nstates):
            outfile.write("%12.6f" % K_kn[k,n])
        outfile.write("\n")        
    outfile.close()            


    #===================================================================================================
    # Compute K at simulated values of s (using full dataset).
    #===================================================================================================    
    
    svalues = s_k
    Kvalues = 0.0 * svalues
    dKvalues = 0.0 * svalues
    for i in range(svalues.size):
        s = svalues[i]
        # Compute effective energy (notice weight is exp[-sK], so contribution to reduced potential is +sK)
        Kvalues[i] = K_kn[i,:].mean()
        g = timeseries.statisticalInefficiency(K_kn[i,:])
        dKvalues[i] = K_kn[i,:].std() / math.sqrt(float(niterations) / g)
        print "g at s = %.3f is %.1f" % (s, g)

    # DEBUG
    outfile = open('simKs.out', 'w')
    for i in range(svalues.size):
        outfile.write("%12.6f %12.6f %12.6f\n" % (svalues[i], Kvalues[i], dKvalues[i]))
    outfile.close()            

    #===================================================================================================
    # Subsample data to obtain uncorrelated samples
    #===================================================================================================
    
    N_k = numpy.zeros(nstates, numpy.int32)
    indices = timeseries.subsampleCorrelatedData(u_n) # indices of uncorrelated samples
    #indices = range(u_n.size) # DEBUG: Assume data are uncorrelated
    N = len(indices) # number of uncorrelated samples
    print "N = %d" % N
    N_k[:] = N      
    u_kln = u_kln[:,:,indices]
    K_kn = K_kn[:,indices]    
    print "number of uncorrelated samples:"
    print N_k
    print ""
    
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

    outfile = open('replica-fields.out', 'w')
    for n in range(niterations):
        for k in range(nstates):
            outfile.write("%12.6f" % s_k[state_nk[n,k]])
        outfile.write("\n")        
    outfile.close()
    
    #===================================================================================================    
    # Write un-wrapped trajectories.
    #===================================================================================================

    write_trajectories = False
    map_into_box = False
    if write_trajectories:
        print "Writing trajectories..."    
        nparticles = 150
        NA = int(0.8 * nparticles)
        principal_component_density=0.96
        sigma       = 0.3405 * units.nanometers # arbitrary reference lengthscale        
        volume = NA * sigma**3 / principal_component_density
        length = volume**(1./3.)    
        atoms = construct_atom_list(nparticles, NA)
        for replica_index in range(nstates):
            filename = "replica-%d.pdb" % replica_index
            trajectory = ncfile.variables['trajectory_coordinates'][replica_index,:,:,:].copy() * 10.0
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
            write_pdb(filename, trajectory, atoms)

        print "Done."

    #===================================================================================================    
    # Write wrapped trajectories.
    #===================================================================================================

    write_trajectories = False
    map_into_box = True
    if write_trajectories:
        print "Writing trajectories..."    
        nparticles = 150
        NA = int(0.8 * nparticles)
        principal_component_density=0.96
        sigma       = 0.3405 * units.nanometers # arbitrary reference lengthscale        
        volume = NA * sigma**3 / principal_component_density
        length = volume**(1./3.)    
        atoms = construct_atom_list(nparticles, NA)
        for replica_index in range(nstates):
            filename = "replica-wrapped-%d.pdb" % replica_index
            trajectory = ncfile.variables['trajectory_coordinates'][replica_index,:,:,:].copy() * 10.0
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
            write_pdb(filename, trajectory, atoms)

        print "Done."

    #===================================================================================================
    # Initialize MBAR.
    #===================================================================================================   
   
    # Initialize MBAR (computing free energy estimates, which may take a while)
    print "Computing free energy differences..."
    #mbar = pymbar.MBAR(u_kln, N_k, verbose=False, method='self-consistent-iteration', relative_tolerance=1.0e-15) # use slow self-consistent-iteration (the default)
    #mbar = pymbar.MBAR(u_kln, N_k, verbose=True, method='Newton-Raphson', initialize='BAR', relative_tolerance=1.0e-9)

    mbar = pymbar.MBAR(u_kln, N_k, verbose=False, method='self-consistent-iteration', initialize='BAR', relative_tolerance=1.0e-2)
    initial_f_k = mbar.f_k
    mbar = pymbar.MBAR(u_kln, N_k, verbose=True, method='Newton-Raphson', initial_f_k=initial_f_k, relative_tolerance=1.0e-9)

    # Get matrix of dimensionless free energy differences and uncertainty estimate.
    print "Computing covariance matrix..."
    (Deltaf_ij, dDeltaf_ij) = mbar.getFreeEnergyDifferences()

    print "Deltaf_ij = "
    print Deltaf_ij
    print "dDeltaf_ij = "
    print dDeltaf_ij
    
    #===================================================================================================
    # Compute <K>_s using MBAR for various slices of iterations, to check convergence
    #===================================================================================================
    computeKs = False
    if computeKs == True:
        slice_size = 1000 # size of slice
        maxiterations = ncfile.variables['states'].shape[0]
        for slice_start in range(0, maxiterations, slice_size):
            slice_end = slice_start + slice_size

            [svalues, Kvalues, dKvalues] = computeKasFunctionOfS(ncfile, mbar, slice_start, slice_end)

            # Write out <K>_s for this slice
            filename = "Ks."+str(slice_start)+"-"+str(slice_end)+".out"
            outfile = open(filename, 'w')
            for i in range(svalues.size):
                outfile.write("%12.6f %12.6f %12.6f\n" % (svalues[i], Kvalues[i], dKvalues[i]))
            outfile.close()


    #===================================================================================================
    # Calculate s* and all observables for the value of s*
    #===================================================================================================
    computeObservables(mbar,u_kln, K_kn)

   
    K_n = K_kn[mbar.indices]
    #===================================================================================================
    # Compute trajectory entropy as a function of s.
    #===================================================================================================    
    computeTrajEntropy = False
    if computeTrajEntropy:

        print "Computing trajectory entropy as a function of s..."

        nbins = 100
        svalues = numpy.linspace(-0.01, 0.08, nbins) # values of reduced s at which to evaluate trajectory entropy

        # Compute 'augmented' u_kln from U.
        u_kln_augmented = numpy.zeros([nstates+nbins, nstates+nbins, N], numpy.float64)
        u_kln_augmented[0:nstates,0:nstates,:] = u_kln
        for k in range(nstates):
            for l in range(nbins):
                u_kln_augmented[k,nstates+l,:] = u_kln[k,0,:] + svalues[l] * K_kn[k,:] * sKfactor

        N_k_augmented = numpy.zeros([nstates+nbins], numpy.int32)
        N_k_augmented[0:nstates] = N_k[0:nstates]
        print "initialize mbar_augmented....."
        mbar_augmented = pymbar.MBAR(u_kln_augmented, N_k_augmented, verbose = True, method = 'Newton-Raphson', initialize='BAR', relative_tolerance=1.0e-13)
        #mbar_augmented = pymbar.MBAR(u_kln_augmented, N_k_augmented, verbose=False, method='self-consistent-iteration', relative_tolerance=1.0e-15)
        [Delta_f_ij, dDelta_f_ij, Delta_u_ij, dDelta_u_ij, Delta_s_ij, dDelta_s_ij] = mbar_augmented.computeEntropyAndEnthalpy(uncertainty_method='svd-ew')

        print "f_ij"
        print Delta_f_ij
        print dDelta_f_ij

        print "u_ij"
        print Delta_u_ij
        print dDelta_u_ij

        print "s_ij"
        print Delta_s_ij
        print dDelta_s_ij

        # DEBUG
        outfile = open('entropy-s.out', 'w')
        index = Delta_s_ij[0,:].argmin()
        for i in range(svalues.size):
            outfile.write("%12.6f %12.6f %12.6f\n" % (svalues[i], Delta_s_ij[index,nstates+i], dDelta_s_ij[index,nstates+i]))
        outfile.close()

        # DEBUG
        outfile = open('f-s.out', 'w')
        index = Delta_f_ij[0,:].argmin()
        for i in range(svalues.size):
            outfile.write("%12.6f %12.6f %12.6f\n" % (svalues[i], Delta_f_ij[index,nstates+i], dDelta_f_ij[index,nstates+i]))
        outfile.close()

    #===================================================================================================
    # Compute trajectory entropy as a function of K.
    #===================================================================================================    

    print "Computing trajectory entropy as a function of K at s = 0..."

    nbins = 100
    sigma = 0.002
    Kvalues = numpy.linspace(K_n.min(), K_n.max(), nbins) # values of reduced K at which to estimate PMF

    # Compute 'augmented' u_kln from U.
    u_kln_augmented = numpy.zeros([nstates+nbins, nstates+nbins, N], numpy.float64)
    u_kln_augmented[0:nstates,0:nstates,:] = u_kln
    for k in range(nstates):
        for l in range(nbins):
            u_kln_augmented[k,nstates+l,:] = u_kln[k,0,:] + (K_kn[k,:] - Kvalues[l])**2 / (2.0 * sigma**2)            
    
    N_k_augmented = numpy.zeros([nstates+nbins], numpy.int32)
    N_k_augmented[0:nstates] = N

    mbar_augmented = pymbar.MBAR(u_kln_augmented, N_k_augmented, verbose=False, method='self-consistent-iteration', relative_tolerance=1.0e-15) 
    [Delta_f_ij, dDelta_f_ij, Delta_u_ij, dDelta_u_ij, Delta_s_ij, dDelta_s_ij] = mbar_augmented.computeEntropyAndEnthalpy(uncertainty_method='svd-ew')

    print "f_ij"
    print Delta_f_ij
    print dDelta_f_ij

    print "u_ij"
    print Delta_u_ij
    print dDelta_u_ij

    print "s_ij"
    print Delta_s_ij
    print dDelta_s_ij

    # DEBUG
    outfile = open('entropy-K.out', 'w')
    index = Delta_s_ij[nstates,nstates:].argmin()
    for i in range(Kvalues.size):
        outfile.write("%12.6f %12.6f %12.6f\n" % (Kvalues[i], Delta_s_ij[nstates+index,nstates+i], dDelta_s_ij[nstates+index,nstates+i]))
    outfile.close()

    #===================================================================================================
    # Clean up.
    #===================================================================================================   

    # Close NetCDF file.
    ncfile.close()
        
    print "Work done  :) go home and have a beer!"
    
    
