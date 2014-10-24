#!/bin/tcsh

# This is needed for ICC variables, for some reason.
resoft

source ~/.cshrc


#  Sample Batch Script for a MVAPICH-Intel job
#
# $HOME/.soft contains:
#
#  @teragrid-basic
#  @globus-4.0
#  @teragrid-dev
#
# $HOME/.mpd.conf contains:
#
#  MPD_SECRETWORD=XXXXXXX     # random alphanumeric chars
#                             # (MUST contain at least one alphabetic char)
#
# (make sure the file .mpd.conf has permissions 700)
#
#  Submit this script using the command: qsub <script_name>
#
#  Use the "qstat" command to check the status of a job.
#
# The following are embedded QSUB options. The syntax is #PBS (the # does
# _not_  denote that the lines are commented out so do not remove).
#
# account : which account to charge (e.g. dpm)
#PBS -A dxi
#
# specify queue
#PBS -q lincoln
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=72:00:00
#
# nodes: number of 8-core nodes
#   ppn: how many cores per node to use (1 through 8)
#       (you are always charged for the entire node)
#PBS -l nodes=5:ppn=2
#
# shell to use
#PBS -S /bin/tcsh
#
# export all my environment variables to the job
#PBS -V
#
# job name (default = name of script file)
#PBS -N kobandersen-AB
#
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
#PBS -o stdout
#
# filename for standard error (default = <job_name>.e<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
#PBS -e stderr
#
# End of embedded QSUB options
#
# set echo               # echo commands before execution; use for debugging
#

# change to job scratch directory, 
# use cdjob <jobid> to go to this directory once the job has started
cd /u/ac/jchodera/code/pytps/kob-andersen/python

# get executable and input files from mass storage
#msscmd "cd dir1, get a.out, mget *.input"

# mss doesn't keep executable bit set, so need to set it on program
#chmod +x a.out

mvapich2-start-mpd
setenv NP `wc -l ${PBS_NODEFILE} | cut -d'/' -f1`
echo "number of processors:"
echo $NP
echo "nodes:"
cat ${PBS_NODEFILE}

setenv MV2_SRQ_SIZE 4000
echo "starting via mpirun..."
mpirun -machinefile ${PBS_NODEFILE} -np ${NP} python KobAndersen-mpi-AB.py

# The mpirun syntax above will assign ranks to nodes in round-robin fashion.
# To get ranks *packed* into nodes, use this syntax:
#
#  mpirun  -machinefile ${PBS_NODEFILE} -np ${NP} a.out

mpdallexit

# save output files back to mass storage
#msscmd "cd dir1, mput *.output"     

