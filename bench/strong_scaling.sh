#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=secondary

#module load openmpi/4.1.0-gcc-7.2.0-pmi2 gcc

#execname=

#file to which processor number and 
pfile='n_procs_strong.txt'
timefile='times_strong.txt'
ofile='output.txt'

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/qmchamm/shared/jacobflib/petsc/arch-linux-c-opt/lib

#matrix dimension
#matrix memory N=15000, complex entries: 1.8 gigabytes
#N=15000
N=150

#number of processes along side of 2-D mesh
rs=(1 2 4  8  10 25  50)

#total number of processes, r*(r+1)/2
ps=(1 3 10 36 55 325 1275)

for i in ${!ps[@]}; do
  printf "%d %d %d\n" $N ${ps[$i]} ${rs[$i]} >> ${ofile}
  #srun execname 1>> timefile
  #command line input for matrix size?
  ./ex2 -symmfastrow_size $N -symmfastcol_size $N -denserow_size $N -densecol_size $N -symmfastcomm_per_dim ${rs[$i]} 1>> ${ofile}
done

#add delimiters between runs in case of multiple runs
#without clearing/saving
echo "---" >> pfile
#echo "---" >> timefile
