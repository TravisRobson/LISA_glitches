#!/bin/sh
#SBATCH -J GliBur1
#SBATCH -o GliBur1_output-%j.dat
#SBATCH -e GliBur1_error-%j.dat
#SBATCH -p xlarge
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --ntasks=13
#SBATCH --mem 64000
#SBATCH --mail-user travis.robson@montana.edu
#SBATCH --mail-type ALL
make

srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_010.dat &

srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_011.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_012.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_013.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_014.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_015.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_016.dat &

srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_017.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_018.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_019.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_020.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_021.dat &
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_022.dat &

wait
