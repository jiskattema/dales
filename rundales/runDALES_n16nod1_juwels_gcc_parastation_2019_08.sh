#!/bin/bash -x
#SBATCH --account=hku28
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --output=./slurm_dales-%j.out
#SBATCH --error=./slurm_error-%j.out
#SBATCH --time=24:00:00
#SBATCH --partition=batch

# Note: The current working directory at this point is
# the directory where sbatch was executed.

echo "Starting run"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# load modules
module use /gpfs/software/juwels/otherstages
module load Stages/2018b
module load GCC/8.2.0 # module load Intel/2019.0.117-GCC-7.3.0
module load ParaStationMPI/5.2.1-1  # module load ParaStationMPI
module load netCDF-Fortran/4.4.4-serial  # module load HDF
module load HDF    # module load netCDF
module load netCDF # module load netCDF-Fortran/4.4.4
module load CMake # module load CMake
# module load Intel/2019.0.117-GCC-7.3.0
# module load ParaStationMPI
# module load HDF
# module load netCDF
# module load netCDF-Fortran/4.4.4
# module load CMake


# set a few variables
    export DALES_EXECUTABLE='dales4'
    export busyfile='runDALES_busy.txt'
    export RUNLOG='./runlog.txt'

echo "# of nodes, tasks:" $SLURM_NNODES $SLURM_NTASKS

# adding a busy file 
echo "runDALES job ${SLURM_JOB_ID} is busy ... don't interfere!" > $busyfile


srun ${DALES_EXECUTABLE}>$RUNLOG
#or:  srun -n $SLURM_NTASKS ${DALES_EXECUTABLE}>$RUNLOG

rm $busyfile

