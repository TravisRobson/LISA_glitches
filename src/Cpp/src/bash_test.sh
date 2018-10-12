#!/bin/sh
#SBATCH -J test
#SBATCH -o test_output-%j.dat
#SBATCH -e test_error-%j.dat
#SBATCH -p priority
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --ntasks=1
#SBATCH --mem 64000
#SBATCH --mail-user travis.robson@montana.edu
#SBATCH --mail-type ALL

make

srun ./g_b --xonly 1 --nmcmc 10000 --nburn 3000 --ptmcmc --inputfile Input_580.dat &

wait