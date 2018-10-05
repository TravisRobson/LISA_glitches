#!/bin/sh
#SBATCH -J enmasse
#SBATCH -o enmasse_output-%j.dat
#SBATCH -e enmasse_error-%j.dat
#SBATCH -p xlarge
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --ntasks=1
#SBATCH --mem 64000
#SBATCH --mail-user travis.robson@montana.edu
#SBATCH --mail-type ALL

make

python CreateInputandBash.py 700 burst 10 15.0e-3 0.011*HOUR priorty ## Generate the Input Files and bash file to run...

sbatch bash_GliBur_700.sh

wait