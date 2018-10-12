#!/usr/bin/env

import sys

# Just some dummy parameters
A         = "1.0e-20"
t0        = "0.5*T"
phi0      = "1.1"

cos_theta = "0.23"
phi       = "2.31"
psi       = "0.45"
ellip     = "0.5"

def print_data_stuff(file, model, snr, f0, tau, which):

    if (model=="burst"): # burst
        file.write(which + " burst" + "\n")
                
    elif (model=="glitch_OP12"):
        file.write(which + " glitch_OP12" + "\n")
        
    elif (model=="glitch_OP21"):
        file.write(which + " glitch_OP21" + "\n")

    elif (model=="glitch_OP13"):
        file.write(which + " glitch_OP13" + "\n")
        
    elif (model=="glitch_OP31"):
        file.write(which + " glitch_OP31" + "\n")         
        
    elif (model=="glitch_OP23"):
        file.write(which + " glitch_OP23" + "\n")         
        
    elif (model=="glitch_OP32"):
        file.write(which + " glitch_OP32" + "\n")
                   
    #####################################################
              
    elif (model=="glitch_AC12"):
        file.write(which + " glitch_AC12" + "\n")
                    
    elif (model=="glitch_AC21"):
        file.write(which + " glitch_AC21" + "\n")          
          
    elif (model=="glitch_AC13"):
        file.write(which + " glitch_AC13" + "\n")        

    elif (model=="glitch_AC31"):
        file.write(which + " glitch_AC31" + "\n")         
        
    elif (model=="glitch_AC23"):
        file.write(which + " glitch_AC23" + "\n")
                  
    elif (model=="glitch_AC32"):
        file.write(which + " glitch_AC32" + "\n")
        
        
    file.write(which + "_snr " + snr + "\n")
    file.write(which + "_A " + A + "\n")
    file.write(which + "_f0 " + f0 + "\n")
    file.write(which + "_t0 " + t0 + "\n")
    file.write(which + "_tau " + tau + "\n")
    file.write(which + "_phi0 " + phi0 + "\n")
    
    if (model=="burst"): # burst
        
        file.write(which + "_cos_theta " + cos_theta + "\n")
        file.write(which + "_phi " + phi + "\n")
        file.write(which + "_psi " + psi + "\n")
        file.write(which + "_ellip " + ellip + "\n")    
        
    return
    
    
    
def main():

    ID = sys.argv[1]
    
    filename = "Input_" + ID + ".dat"
    file  = open(filename, "w")
    
    model = sys.argv[2]
    data  = sys.argv[3]
    snr   = sys.argv[4]
    f0    = sys.argv[5]
    tau   = sys.argv[6]

    print_data_stuff(file, model, snr, f0, tau, "model")
    file.write("\n")  
    
    print_data_stuff(file, data, snr, f0, tau, "data")
    file.write("\n")      

    file.write("File_true_waveform ")
    file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
    file.write("Burst_True_Waveform_" + ID + ".dat\n")

    file.write("File_cold_chain ")
    file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
    file.write("Burst_Cold_Chain_" + ID + ".dat\n")

    file.write("File_T1_chain ")
    file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
    file.write("Burst_T1_Chain_" + ID + ".dat\n")        

    file.write("File_Hot_chain ")
    file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
    file.write("Burst_Hot_Chain_" + ID + ".dat\n")

    file.write("File_logL ")
    file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
    file.write("Burst_logL_" + ID + ".dat\n")

    file.write("File_IDs ")
    file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
    file.write("Burst_IDs_" + ID + ".dat\n")  

    file.write("File_Temps ")
    file.write("/mnt/lustrefs/work/travis.robson/CppGlitchBurst/Output/")
    file.write("Burst_Temps_" + ID + ".dat\n")
        
    file.close()  

    return
    
    

if __name__ == "__main__":
    main()


