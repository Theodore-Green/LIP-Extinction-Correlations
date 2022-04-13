#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=03:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=theodore.green@princeton.edu ## Add own email address

cd $SLURM_SUBMIT_DIR
module load openmpi/4.0.1-gnu4.8.5

## Choose to run LIP correlations or impact correlations
mpiexec -np 16 julia ./MPIReduceCorrelationHistogramLIPs.jl
## mpiexec -np 16 julia ./MPIReduceCorrelationHistogramImpacts.jl
