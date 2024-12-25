#!/bin/bash

 

#SBATCH --partition=highmem
##SBATCH --exclusive
#SBATCH -n 128
#SBATCH --ntasks-per-node=64

#SBATCH --mem=500000

 

source ~/.bashrc
module load icc/20.4
module load openmpi/4.1.2

 

export OMPI_NUM_THREADS=1
mpirun ./a.out
