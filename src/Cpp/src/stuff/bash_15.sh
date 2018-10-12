#!/bin/sh
#SBATCH --workdir="/mnt/lustrefs/work/travis.robson/CppGlitchBurst"
#SBATCH -J GB_15
#SBATCH -o GB_15_output-%j.dat
#SBATCH -e GB_15_error-%j.dat
#SBATCH -p priority
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --ntasks=1
#SBATCH --mem 64000
#SBATCH --mail-user travis.robson@montana.edu
#SBATCH --mail-type ALL

make

python MakeInput.py 15 burst glitch_OP12 8 15e-3 0.05*HOUR
wait

srun ./g_b --xonly 0 --nmcmc 10000 --nburn 3000 --ptmcmc --inputfile Input_15.dat

wait

file="Output/Burst_logL_15.dat"
read -r firstline<"$file"

for word in $firstline
do
	srun ./thermo_error Output/Burst_logL_15.dat $word 10000 BF_15.dat
	break
done
