#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=00:10:00
#SBATCH --job-name=TDDFT

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

#Single point with RI
ridft > ridft.out

# TDDFT , RPA  Absorptionsspektrum
escf > escf.out
