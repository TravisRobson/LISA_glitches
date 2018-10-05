#!/bin/sh
#SBATCH -J bf
#SBATCH -o bf_output-%j.dat
#SBATCH -e bf_error-%j.dat
#SBATCH -p express
#SBATCH -N 1
#SBATCH -t 20:00
#SBATCH --ntasks=10
#SBATCH --mem 64000
#SBATCH --mail-user travis.robson@montana.edu
#SBATCH --mail-type ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

gcc -O2 thermo_error.c -lm -llapack -lblas -o thermo_error



file="Output/Burst_logL_220.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_220.dat $word 10000 BF_220.dat &
    break
done

file="Output/Burst_logL_240.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_240.dat $word 10000 BF_240.dat &
    break
done

file="Output/Burst_logL_260.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_260.dat $word 10000 BF_260.dat &
    break
done

file="Output/Burst_logL_280.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_280.dat $word 10000 BF_280.dat &
    break
done

file="Output/Burst_logL_300.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_300.dat $word 10000 BF_300.dat &
    break
done

## file="Output/Burst_logL_320.dat"
## read -r firstline<"$file"
## 
## for word in $firstline
## do  
##     srun ./thermo_error Output/Burst_logL_320.dat $word 10000 BF_320.dat &
##     break
## done

######################################
file="Output/Burst_logL_221.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_221.dat $word 10000 BF_221.dat &
    break
done

file="Output/Burst_logL_241.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_241.dat $word 10000 BF_241.dat &
    break
done

file="Output/Burst_logL_261.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_261.dat $word 10000 BF_261.dat &
    break
done

file="Output/Burst_logL_281.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_281.dat $word 10000 BF_281.dat &
    break
done

file="Output/Burst_logL_301.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_301.dat $word 10000 BF_301.dat &
    break
done

## file="Output/Burst_logL_321.dat"
## read -r firstline<"$file"
## 
## for word in $firstline
## do  
##     srun ./thermo_error Output/Burst_logL_321.dat $word 10000 BF_321.dat &
##     break
## done

wait

rm fit*.dat

