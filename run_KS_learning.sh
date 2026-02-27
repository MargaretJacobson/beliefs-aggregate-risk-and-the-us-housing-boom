#!/bin/bash
#SBATCH --ntasks=1                  #number of times to execute commands
#SBATCH  --cpus-per-task=9         #number of cpus
#SBATCH --job-name=Jacobson_GE      #Job name (not that important unless running multiple jobs)
#SBATCH --output=out_learn.txt             #Can specify an output and error file, otherwise it just puts them in one file that has the job
#SBATCH --error=error_learn.txt            #number as a name
#SBATCH --mem-per-cpu=50G         #Memory per CPU 18G
#SBATCH --time=3:20:00                #Time needed to execute program
# SBATCH --partition=vhm
#SBATCH --nice==-10

cd `pathf $SLURM_SUBMIT_DIR`

#Remove old files
rm -f *.o *.mod


#Load software we need, called modules

###################################################################
#        Intel Fortran
###################################################################
# Need to run Intel Fortran instead of Gfortran because arrays have rank great than 7 which Gfortran can't handle
. /etc/profile.d/modules.sh
module load comp/intel/2019.1
module load lib/mkl/2019.1


## Parallelization parameters
# Without beliefs, set out to 1
# With beliefs, use nested parallelization
let outer=9 # Threads in outer loop
let threads=20 #Threads in inner loop
export OMP_NUM_THREADS=$outer,$threads #$PBS_NUM_NODES*$PBS_NUM_PPN    #sets number of threads, otherwise would just be ppn
export OMP_PlACES=CORES
export OMP_NESTED=TRUE                                          #turns on nested parallelism
export OMP_DYNAMIC=TRUE
ulimit -s unlimited                                             #sets unlimited stack size

make

#export KMP_HANDLE_SIGNALS=0
KMP_STACKSIZE=100G OMP_NUM_THREADS=$outer,$threads ./main.o    


