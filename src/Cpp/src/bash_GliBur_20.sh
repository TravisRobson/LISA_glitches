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

srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_20.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_21.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_22.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_23.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_24.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_25.dat&

srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_26.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_27.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_28.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_29.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_30.dat&
srun ./g_b --nmcmc 10000 --nburn 1000 --ptmcmc --inputfile Input_31.dat&

wait
