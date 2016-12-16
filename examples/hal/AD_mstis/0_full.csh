#!/bin/tcsh
#  Batch script for single thread GPU job.
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=36:00:00
#
# memory : the amount of available memory
#PBS -l mem=4G
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify queue
#PBS -q batch
#
# nodes: number of nodes
#   ppn: number of processes per node
#  gpus: number of gpus per node
#  GPUs are in 'exclusive' mode by default, but 'shared' keyword sets them to shared mode.
#PBS -l nodes=1:ppn=1:gpus=1:shared
#
# export all my environment variables to the job
#PBS -V
#
# job name (default = name of script file)
#PBS -N AD_mstis_0_full
#
# mail settings (one or more characters)
# email is sent to local user, unless another email address is specified with PBS -M option
# n: do not send mail
# a: send mail if job is aborted
# b: send mail when job begins execution
# e: send mail when job terminates
#PBS -m n
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
##PBS -o myoutput

# Change to working directory used for job submission
cd $PBS_O_WORKDIR

# Launch my program.
pwd
setenv OPS_PROJECT_NUMBER 1
python /cbio/jclab/home/prinzj/projects/ops/AD_mstis/AD_mstis_A1_resetup.py
setenv OPS_PROJECT_NUMBER 2
python /cbio/jclab/home/prinzj/projects/ops/AD_mstis/AD_mstis_2_run.py
python /cbio/jclab/home/prinzj/projects/ops/AD_mstis/AD_mstis_5_split.py
