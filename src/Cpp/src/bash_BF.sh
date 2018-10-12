#!/bin/sh
#SBATCH -J bf
#SBATCH -o bf_output-%j.dat
#SBATCH -e bf_error-%j.dat
#SBATCH -p priority
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --ntasks=10
#SBATCH --mem 64000
#SBATCH --mail-user travis.robson@montana.edu
#SBATCH --mail-type ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

gcc -O2 thermo_error.c -lm -llapack -lblas -o thermo_error

num=5
mum=6


file="Output/Burst_logL_"$num"20.dat"
echo $file
echo BF_"$num"20.dat
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$num"20.dat $word 50000 BF_"$num"20.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat


file="Output/Burst_logL_"$num"40.dat"
read -r firstline<"$file"

echo $file
echo BF_"$num"40.dat

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$num"40.dat $word 50000 BF_"$num"40.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat



file="Output/Burst_logL_"$num"60.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$num"60.dat $word 50000 BF_"$num"60.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat




file="Output/Burst_logL_"$num"80.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$num"80.dat $word 50000 BF_"$num"80.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat


file="Output/Burst_logL_"$mum"00.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$mum"00.dat $word 50000 BF_"$mum"00.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat



######################################
file="Output/Burst_logL_"$num"21.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$num"21.dat $word 50000 BF_"$num"21.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat


file="Output/Burst_logL_"$num"41.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$num"41.dat $word 50000 BF_"$num"41.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat


file="Output/Burst_logL_"$num"61.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$num"61.dat $word 50000 BF_"$num"61.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat



file="Output/Burst_logL_"$num"81.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$num"81.dat $word 50000 BF_"$num"81.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat


file="Output/Burst_logL_"$mum"01.dat"
read -r firstline<"$file"

for word in $firstline
do  
    srun ./thermo_error Output/Burst_logL_"$mum"01.dat $word 50000 BF_"$mum"01.dat &
    break
done

wait
rm fit*.dat
rm boost.dat
rm chain.dat
rm changes.dat
rm initialfit.dat
rm integrand.dat
rm rawintegrand.dat
rm pchain.dat
rm slopes.dat




