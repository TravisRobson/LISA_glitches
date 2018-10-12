#!/usr/bin/env

import sys



def main():

    cluster   = sys.argv[1] # xlarge? priority?
    xonlyFlag = sys.argv[2] # how many channels to analyze
    ID        = sys.argv[3] # which input file?
    model     = sys.argv[4]
    data      = sys.argv[5]
    snr       = sys.argv[6]
    f0        = sys.argv[7]
    tau       = sys.argv[8] 
    
    filename = "bash_" + ID + ".sh"
    file = open(filename, "w")
    
    file.write("#!/bin/sh\n")
    file.write("#SBATCH --workdir=\"/mnt/lustrefs/work/travis.robson/CppGlitchBurst\"\n")
    file.write("#SBATCH -J GB_"+ID+"\n")
    file.write("#SBATCH -o GB_"+ID+"_output-%j.dat\n")
    file.write("#SBATCH -e GB_"+ID+"_error-%j.dat\n")
    file.write("#SBATCH -p " + cluster+"\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -t 72:00:00\n")
    file.write("#SBATCH --ntasks=1\n")
    file.write("#SBATCH --mem 64000\n")
    file.write("#SBATCH --mail-user travis.robson@montana.edu\n")
    file.write("#SBATCH --mail-type ALL\n")

    file.write("\n")
    file.write("make\n\n")
    
    file.write("python MakeInput.py " + ID + " " + model + " " + data + " " + snr + " " + f0 + " " + tau + "\n")
                            
    file.write("wait\n\n")
    
    file.write("srun ./g_b --xonly " + xonlyFlag + " --nmcmc 10000 --nburn 3000 --ptmcmc --inputfile Input_" + ID + ".dat")
    
    file.write("\n\nwait\n\n")
    file.write("")
#     file.write("gcc -O2 thermo_error.c -lm -llapack -lblas -o thermo_error \n\n")
    
    file.write("file=\"Output/Burst_logL_" + ID + ".dat\"\n")
    file.write("read -r firstline<\"$file\"\n\n")
    file.write("for word in $firstline\n")
    file.write("do\n")
    file.write("\tsrun ./thermo_error Output/Burst_logL_" + ID + ".dat $word 10000 BF_" + ID + ".dat\n")
    file.write("\tbreak\n")
    file.write("done\n")
    
    file.close()

    return
    
    

if __name__ == "__main__":
    main()



# make
# 
# srun ./g_b --xonly 1 --nmcmc 10000 --nburn 3000 --ptmcmc --inputfile Input_580.dat &
# 
# wait