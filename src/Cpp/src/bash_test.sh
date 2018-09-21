#!/bin/sh
#SBATCH -J test
#SBATCH -o test_output-%j.dat
#SBATCH -e test_error-%j.dat
#SBATCH -p xlarge
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --ntasks=2
#SBATCH --mem 64000
#SBATCH --mail-user travis.robson@montana.edu
#SBATCH --mail-type ALL

make

srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_003.dat &
srun ./g_b --nmcmc 100000 --nburn 1000 --ptmcmc --inputfile Input_004.dat &

wait