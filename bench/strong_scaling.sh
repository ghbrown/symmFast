#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --nodes=40
#SBATCH --ntasks-per-node=16
#SBATCH --partition=secondary

module load openmpi/4.1.0-gcc-7.2.0-pmi2 gcc

use_real=true

if [ "$use_real" == true ]; then
    # REAL
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/qmchamm/shared/jacobflib/petsc/arch-linux-c-opt/lib
    executable=ex2
    ofile='output_real.txt'
else
    # COMPLEX
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/qmchamm/shared/jacobflib/petsc/arch-linux-c-complex-opt/lib
    executable=ex2_complex
    ofile='output_complex.txt'
fi

#matrix dimension
#matrix memory N=15000, complex entries: 1.8 gigabytes
N=15000

#number of processes along side of 2-D mesh
rs=(1 2 3 4  5  6  8  10 15  20  25  30  35)

#total number of processes, r*(r+1)/2
ps=(1 3 6 10 15 21 36 55 120 210 325 465 630)

#total number of trials/iterations for each case
nits=(10 10 20 100 100 100 100 100 100 100 100 100 100)

for i in ${!ps[@]}; do
  printf "%d %d %d %d\n" $N ${ps[$i]} ${rs[$i]} ${nits[$i]} >> ${ofile}
  srun -n ${ps[$i]} ./${executable} -symmfastrow_size $N -symmfastcol_size $N -denserow_size $N -densecol_size $N -symmfastcomm_per_dim ${rs[$i]} -num_iterations ${nits[$i]} 1>> ${ofile}
done



