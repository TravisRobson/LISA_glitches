#!/bin/sh
#SBATCH --workdir="/mnt/lustrefs/work/travis.robson/CppGlitchBurst"
#SBATCH -J GB_30
#SBATCH -o GB_30_output-%j.dat
#SBATCH -e GB_30_error-%j.dat
#SBATCH -p priority
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --ntasks=1
#SBATCH --mem 64000
#SBATCH --mail-user travis.robson@montana.edu
#SBATCH --mail-type ALL

make

python MakeInput.py 30 burst glitch_OP12 5 15e-3 0.05*HOUR
wait

srun ./g_b --xonly 1 --nmcmc 10000 --nburn 3000 --ptmcmc --inputfile Input_30.dat

wait

file="Output/Burst_logL_30.dat"
read -r firstline<"$file"

for word in $firstline
do
	srun ./thermo_error Output/Burst_logL_30.dat $word 10000 BF_30.dat
	break
done
