#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=10:00:00
#SBATCH --job-name=GW

#################################################################
####################    TURBOMOLE       #########################
#################################################################
# parallel calculation
export PARNODES=$SLURM_NPROCS
export PARA_ARCH=SMP
export OMP_NUM_THREADS=$SLURM_NPROCS
export MKL_NUM_THREADS=$SLURM_NPROCS
export TM_PAR_OMP=on
module load turbomole/7.4.1

#Single point with RI (for G0W0)
ridft > ridft.out

#Single point (for GW)
#dscf > dscf.out

# TDDFT, RPA, GW...
escf > escf.out
