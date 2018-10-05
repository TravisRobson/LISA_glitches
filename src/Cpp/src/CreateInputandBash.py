#!/usr/bin/env

import sys

# python CreateInputBash.py [idx_start] [model] [snr] [f0] [tau]

# Just some dummy parameters
A         = "1.0e-20"
t0        = "0.5*T"
phi0      = "1.1"

cos_theta = "0.23"
phi       = "2.31"
psi       = "0.45"
ellip     = "0.5"

def print_data_stuff(data_model, file, snr, f0, tau):

    if (data_model=="burst"): # burst
        file.write("data burst" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")
        
        file.write("data_cos_theta " + cos_theta + "\n")
        file.write("data_phi " + phi + "\n")
        file.write("data_psi " + psi + "\n")
        file.write("data_ellip " + ellip + "\n")
                
    elif (data_model=="glitch_OP12"):
        file.write("data glitch_OP12" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")
        
    elif (data_model=="glitch_OP21"):
        file.write("data glitch_OP21" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")

    elif (data_model=="glitch_OP13"):
        file.write("data glitch_OP13" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")
        
    elif (data_model=="glitch_OP31"):
        file.write("data glitch_OP31" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")            
        
    elif (data_model=="glitch_OP23"):
        file.write("data glitch_OP23" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")            
        
    elif (data_model=="glitch_OP32"):
        file.write("data glitch_OP32" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n") 
                   
    #####################################################
              
    elif (data_model=="glitch_AC12"):
        file.write("data glitch_AC12" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")            
                    
    elif (data_model=="glitch_AC21"):
        file.write("data glitch_AC21" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")            
          
    elif (data_model=="glitch_AC13"):
        file.write("model glitch_AC13" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")            

    elif (data_model=="glitch_AC31"):
        file.write("model glitch_AC31" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")          
        
    elif (data_model=="glitch_AC23"):
        file.write("model glitch_AC23" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n")  
                  
    elif (data_model=="glitch_AC32"):
        file.write("model glitch_AC32" + "\n")
        file.write("data_snr " + snr + "\n")
        file.write("data_A " + A + "\n")
        file.write("data_f0 " + f0 + "\n")
        file.write("data_t0 " + t0 + "\n")
        file.write("data_tau " + tau + "\n")
        file.write("data_phi0 " + phi0 + "\n") 
        
    return

def main():
    
    idx_start = int(sys.argv[1]) # the beginning ID integer
    N_models = 13 # number of models and therefore Input_XXX.dat files
    
    file_input_start = "Input_"
    file_input_end   = ".dat"
    
    data_model = sys.argv[2]
    snr        = sys.argv[3]
    f0         = sys.argv[4]
    tau        = sys.argv[5]
    
    # print out the model information
    for i in range(N_models):
        file_name = file_input_start + str(idx_start + i) + file_input_end
        file = open(file_name, "w")
        
        if (i==0): # burst
            file.write("model burst" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")
            
            file.write("model_cos_theta " + cos_theta + "\n")
            file.write("model_phi " + phi + "\n")
            file.write("model_psi " + psi + "\n")
            file.write("model_ellip " + ellip + "\n")
                    
        elif (i==1):
            file.write("model glitch_OP12" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")
            
        elif (i==2):
            file.write("model glitch_OP21" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")

        elif (i==3):
            file.write("model glitch_OP13" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")
            
        elif (i==4):
            file.write("model glitch_OP31" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")            
            
        elif (i==5):
            file.write("model glitch_OP23" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")            
            
        elif (i==6):
            file.write("model glitch_OP32" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n") 
                       
        #####################################################
                  
        elif (i==7):
            file.write("model glitch_AC12" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")            
                        
        elif (i==8):
            file.write("model glitch_AC21" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")            
              
        elif (i==9):
            file.write("model glitch_AC13" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")            

        elif (i==10):
            file.write("model glitch_AC31" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")           
            
        elif (i==11):
            file.write("model glitch_AC23" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")  
                      
        elif (i==12):
            file.write("model glitch_AC32" + "\n")
            file.write("model_snr " + snr + "\n")
            file.write("model_A " + A + "\n")
            file.write("model_f0 " + f0 + "\n")
            file.write("model_t0 " + t0 + "\n")
            file.write("model_tau " + tau + "\n")
            file.write("model_phi0 " + phi0 + "\n")

        # Now print data stuff
        file.write("\n")
        print_data_stuff(data_model, file, snr, f0, tau)
        file.write("\n")    

        file.write("File_true_waveform ")
        file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
        file.write("Burst_True_Waveform_" + str(idx_start + i) + ".dat\n")
           
        file.write("File_cold_chain ")
        file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
        file.write("Burst_Cold_Chain_" + str(idx_start + i) + ".dat\n")

        file.write("File_T1_chain ")
        file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
        file.write("Burst_T1_Chain_" + str(idx_start + i) + ".dat\n")        
        
        file.write("File_Hot_chain ")
        file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
        file.write("Burst_Hot_Chain_" + str(idx_start + i) + ".dat\n")
           
        file.write("File_logL ")
        file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
        file.write("Burst_logL_" + str(idx_start + i) + ".dat\n")

        file.write("File_IDs ")
        file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
        file.write("Burst_IDs_" + str(idx_start + i) + ".dat\n")  
        
        file.write("File_Temps ")
        file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
        file.write("Burst_Temps_" + str(idx_start + i) + ".dat\n")
        
    file.close()
    
    file_name = "bash_GliBur_" + str(idx_start) + ".sh"
    file = open(file_name, "w")

    file.write("#!/bin/sh\n")
    file.write("#SBATCH -J GliBur" + str(idx_start) + "\n")
    file.write("#SBATCH -o GliBur" + str(idx_start) + "_output-%j.dat\n")
    file.write("#SBATCH -e GliBur" + str(idx_start) + "_error-%j.dat\n")
    file.write("#SBATCH -p " + sys.argv[6] + "\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -t 72:00:00\n")
    file.write("#SBATCH --ntasks=13\n")
    file.write("#SBATCH --mem 64000\n")
    file.write("#SBATCH --mail-user travis.robson@montana.edu\n")
    file.write("#SBATCH --mail-type ALL\n")

    file.write("make\n\n")

#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 1000 --ptmcmc --inputfile Input_" + str(idx_start + 0) + ".dat&\n\n")
# 
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 1) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 2) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 3) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 4) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 5) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 6) + ".dat&\n\n")
# 
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 7) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 8) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 9) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 10) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 11) + ".dat&\n")
#     file.write("srun --exclusive ./g_b --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 12) + ".dat&\n\n")


    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 0) + ".dat&\n\n")

    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 1) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 2) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 3) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 4) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 5) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 6) + ".dat&\n\n")

    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 7) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 8) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 9) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 10) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 11) + ".dat&\n")
    file.write("srun --exclusive ./g_b --xonly 1 --nmcmc 50000 --nburn 10000 --ptmcmc --inputfile Input_" + str(idx_start + 12) + ".dat&\n\n")
    
    file.write("wait\n")
    
    file.close()
    
    return

if __name__ == "__main__":
    main()