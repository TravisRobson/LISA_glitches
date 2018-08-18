import numpy as np

import Wavelet as wv

import LISA as l
import TDI as td

IDX_comp_id = 5

Week = 7.*24*3600
mHz = 1.0e-3

COMP_ID_LS = ['OP12', 'OP21', 'OP13', 'OP31', 'OP23', 'OP32', \
			  'AC12', 'AC21', 'AC13', 'AC31', 'AC23', 'AC32']
		
def set_glitch_properties(self):
    """ determine the disturbance type and the Laser Link on which it afflicts """
    
    if (self.comp_id == 'OP12'):
        instrument_glitch_type = 'Optical Path'
        SC_on    = 1
        SC_point = 2
    elif (self.comp_id == 'OP21'):
        instrument_glitch_type = 'Optical Path'
        SC_on    = 2
        SC_point = 1
    elif (self.comp_id == 'OP13'):
        instrument_glitch_type = 'Optical Path'
        SC_on    = 1
        SC_point = 3
    elif (self.comp_id == 'OP31'):
        instrument_glitch_type = 'Optical Path'
        SC_on    = 3
        SC_point = 1
    elif (self.comp_id == 'OP23'):
        instrument_glitch_type = 'Optical Path'
        SC_on    = 2
        SC_point = 3
    elif (self.comp_id == 'OP32'):
        instrument_glitch_type = 'Optical Path'
        SC_on    = 3
        SC_point = 2    

    elif (self.comp_id == 'AC12'):
        instrument_glitch_type = 'Acceleration'
        SC_on    = 1
        SC_point = 2
    elif (self.comp_id == 'AC21'):
        instrument_glitch_type = 'Acceleration'
        SC_on    = 2
        SC_point = 1
    elif (self.comp_id == 'AC13'):
        instrument_glitch_type = 'Acceleration'
        SC_on    = 1
        SC_point = 3
    elif (self.comp_id == 'AC31'):
        instrument_glitch_type = 'Acceleration'
        SC_on    = 3
        SC_point = 1
    elif (self.comp_id == 'AC23'):
        instrument_glitch_type = 'Acceleration'
        SC_on    = 2
        SC_point = 3
    elif (self.comp_id == 'AC32'):
        instrument_glitch_type = 'Acceleration'
        SC_on    = 3
        SC_point = 2
        
    return instrument_glitch_type, SC_on, SC_point	

def calc_TDI(self):
    """ Generate the TDI channels for an instrumental glitch """

    instrument_glitch_type, SC_on, SC_point = self.set_props()
            
    Wave  = self.Wavelet
    Orbit = self.Orbit
    
    t = np.arange(0.0, Orbit.Tobs, Orbit.dt) # Todo: Don't need Orbit, its in Wavelet
    N = len(t)

    # create empty phase first
    # Todo: be smarter here to save time, i.e. don't twice make phases
    p12 = td.Phase(1,2, t, np.zeros(N,dtype=np.complex_))
    p21 = td.Phase(2,1, t, np.zeros(N,dtype=np.complex_))
    p13 = td.Phase(1,3, t, np.zeros(N,dtype=np.complex_))
    p31 = td.Phase(3,1, t, np.zeros(N,dtype=np.complex_))
    p23 = td.Phase(2,3, t, np.zeros(N,dtype=np.complex_))
    p32 = td.Phase(3,2, t, np.zeros(N,dtype=np.complex_))

	# Handle a Laser Phase glitch
    if (instrument_glitch_type == 'Laser Phase'):
        if (SC_point != None):
            raise ValueError("Lase noise is from the laser on one S/C")
        
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = wv.Wavelet(Wave.A, Wave.f0, Wave.tau, \
                    Wave.t0 + Orbit.L/l.Clight, Wave.phi0, Orbit)
        wave_temp.calc_Psi()
        wave_temp.make_padded_Psi(t) 

        if (SC_on == 1):
            p12 = td.Phase(1,2, t, +wave_temp.Psi_padded)
            p13 = td.Phase(1,3, t, +wave_temp.Psi_padded)

            p21 = td.Phase(2,1, t, -Wave.Psi_padded)
            p31 = td.Phase(3,1, t, -Wave.Psi_padded)

        elif (SC_on == 2):
            p21 = td.Phase(2,1, t, +wave_temp.Psi_padded)
            p23 = td.Phase(2,3, t, +wave_temp.Psi_padded)

            p12 = td.Phase(1,2, t, -Wave.Psi_padded)
            p32 = td.Phase(3,2, t, -Wave.Psi_padded)

        elif (SC_on == 3):
            p31 = td.Phase(3,1, t, +wave_temp.Psi_padded)
            p32 = td.Phase(3,2, t, +wave_temp.Psi_padded)

            p13 = td.Phase(1,3, t, -Wave.Psi_padded)
            p23 = td.Phase(2,3, t, -Wave.Psi_padded)

        else:
            raise ValueError("Invalid SC_on!!!")   


	# Handle an Optical Path glitch
    elif (instrument_glitch_type == 'Optical Path'):
        if (SC_point == None):
            raise ValueError("Optical Path needs an S/C pointed towards.")
        if (SC_point == SC_on):
            raise ValueError("S/C pointed towards must differ from S/C glitch afflicts.")

        # There is only one Laser Link which gets polluted. Do so.
        if (SC_on == 1 and SC_point == 2):
            p12 = td.Phase(1,2, t, Wave.Psi_padded)
            
        elif (SC_on == 2 and SC_point == 1):
            p21 = td.Phase(2,1, t, Wave.Psi_padded)

        elif (SC_on == 1 and SC_point == 3):
            p13 = td.Phase(1,3, t, Wave.Psi_padded)

        elif (SC_on == 3 and SC_point == 1):
            p31 = td.Phase(3,1, t, Wave.Psi_padded)

        elif (SC_on == 2 and SC_point == 3):
            p23 = td.Phase(2,3, t, Wave.Psi_padded)

        elif (SC_on == 3 and SC_point == 2):
            p32 = td.Phase(3,2, t, Wave.Psi_padded)
        else:
            raise ValueError("Invalid SC_on and/or Sc_point!!!")

    # Handle an acceleration noise glitch   
    elif (instrument_glitch_type == 'Acceleration'):
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = wv.Wavelet(Wave.A, Wave.f0, Wave.tau, \
                    Wave.t0 + Orbit.L/l.Clight, Wave.phi0, Orbit)
        wave_temp.calc_Psi()
        wave_temp.make_padded_Psi(t)

        if (SC_on == 1 and SC_point == 2):
            p12 = td.Phase(1,2, t, -Wave.Psi_padded)
            p21 = td.Phase(2,1, t, +wave_temp.Psi_padded)

        elif (SC_on == 1 and SC_point == 3):
            p13 = td.Phase(1,3, t, -Wave.Psi_padded)
            p31 = td.Phase(3,1, t, +wave_temp.Psi_padded)

        elif (SC_on == 2 and SC_point == 1):
            p21 = td.Phase(2,1, t, -Wave.Psi_padded)
            p12 = td.Phase(1,2, t, +wave_temp.Psi_padded)

        elif (SC_on == 3 and SC_point == 1):
            p31 = td.Phase(3,1, t, -Wave.Psi_padded)
            p13 = td.Phase(1,3, t, +wave_temp.Psi_padded)

        elif (SC_on == 2 and SC_point == 3):
            p23 = td.Phase(2,3, t, -Wave.Psi_padded)
            p32 = td.Phase(3,2, t, +wave_temp.Psi_padded)

        elif (SC_on == 3 and SC_point == 2):
            p32 = td.Phase(3,2, t, -Wave.Psi_padded)
            p23 = td.Phase(2,3, t, +wave_temp.Psi_padded)

        else:
            raise ValueError("Invalid SC_on and/or Sc_point!!!")

    else:
        raise ValueError("Unexpected instrument glitch type. Choose from 'Laser Phase', 'Optical Path', or 'Acceleration")
        
    p12.phi = p12.phi.real
    p21.phi = p21.phi.real
    p13.phi = p13.phi.real
    p31.phi = p31.phi.real
    p23.phi = p23.phi.real
    p32.phi = p32.phi.real
    
    # Fourier transform the time-domain phases
    p12.FT_phase(Orbit)
    p21.FT_phase(Orbit)
    p13.FT_phase(Orbit)
    p31.FT_phase(Orbit)
    p23.FT_phase(Orbit)
    p32.FT_phase(Orbit)
    
    self.TDI = td.TDI(p12, p21, p13, p31, p23, p32, Orbit)
    self.TDI.f_min = Wave.f_min
    self.TDI.f_max = Wave.f_max
    
    return

def calc_glitch_snr(self):
    """ Calculate the SNR of the glitch """
    
    #  Todo: create flag such that X channel or AET is an option    
    snr_list = np.sum(self.TDI.get_TDI_snr(self.f_min, self.f_max))

    self.SNR = np.sqrt(snr_list)

    return 
    
def adjust_to_target_snr(self, target):
    """ Adjust the SNR and amplitude to hit target SNR, and its TDI data """
    
    amp = target/self.SNR
    
    self.Wavelet.A *= amp
    self.params[wv.IDX_lnA] = np.log(self.Wavelet.A) # this might already update automatically with python
    self.Wavelet.calc_Psi()
    t = np.arange(0.0, self.Orbit.Tobs, self.Orbit.dt)
    self.Wavelet.make_padded_Psi(t)
    
    self.calc_TDI()
    self.calc_snr()
    
    return
    
def param_vec_to_names(Wavelet, params):
    """ take the parameter vector and update Wavelet parameters (i.e. named ones) """
    
    Wavelet.A    = np.exp(params[wv.IDX_lnA])
    Wavelet.f0   = mHz*params[wv.IDX_f0]
    Wavelet.t0   = Week*params[wv.IDX_t0]
    Wavelet.tau  = Week*params[wv.IDX_tau]
    Wavelet.phi0 = params[wv.IDX_phi0]
    
#     comp_id = COMP_ID_LS[int(params[IDX_comp_id])]
    
    return   
    
    
def calc_Fisher(self):
    """ Calculate the Fisher matrix for the glitch object """
    
    ep = 1.0e-6
    
    for i in range(len(COMP_ID_LS)):
        if (COMP_ID_LS[i] == self.comp_id):
            comp_num = i
            break
    
    t = np.arange(0.0, self.Orbit.Tobs, self.Orbit.dt)  

    Fisher = np.zeros((IDX_comp_id, IDX_comp_id))
    
    # this is not an efficient way to calculate the Fisher matrix
       
    for i in range(IDX_comp_id):
    
        self.params[i] += ep
        
        wave_p_LHS = wv.Wavelet(self.Wavelet.A, self.Wavelet.f0, self.Wavelet.tau, self.Wavelet.t0, self.Wavelet.phi0, self.Orbit)
        param_vec_to_names(wave_p_LHS, self.params)
        
        self.params[i] -= 2*ep
        
        wave_m_LHS = wv.Wavelet(self.Wavelet.A, self.Wavelet.f0, self.Wavelet.tau, self.Wavelet.t0, self.Wavelet.phi0, self.Orbit)
        param_vec_to_names(wave_m_LHS, self.params)
               
        self.params[i] += ep    # to return parameter to initial value           
              
        wave_p_LHS.calc_Psi()
        wave_p_LHS.make_padded_Psi(t)
        wave_m_LHS.calc_Psi()
        wave_m_LHS.make_padded_Psi(t)
        
        glitch_p_LHS = Glitch(wave_p_LHS, comp_num, self.Orbit)
        glitch_m_LHS = Glitch(wave_m_LHS, comp_num, self.Orbit)
        
        glitch_p_LHS.calc_TDI()
        glitch_m_LHS.calc_TDI()
        
        # take the derivatives and stuff in TDI of plus glitch
        glitch_p_LHS.TDI.A = (glitch_p_LHS.TDI.A - glitch_m_LHS.TDI.A)/2/ep
        glitch_p_LHS.TDI.E = (glitch_p_LHS.TDI.E - glitch_m_LHS.TDI.E)/2/ep
        glitch_p_LHS.TDI.T = (glitch_p_LHS.TDI.T - glitch_m_LHS.TDI.T)/2/ep    
        
        for j in range(i, IDX_comp_id):
            
            # RHS 
            self.params[j] += ep
            
            wave_p_RHS = wv.Wavelet(self.Wavelet.A, self.Wavelet.f0, self.Wavelet.tau, self.Wavelet.t0, self.Wavelet.phi0, self.Orbit)
            param_vec_to_names(wave_p_RHS, self.params)

            self.params[j] -= 2*ep

            wave_m_RHS = wv.Wavelet(self.Wavelet.A, self.Wavelet.f0, self.Wavelet.tau, self.Wavelet.t0, self.Wavelet.phi0, self.Orbit)
            param_vec_to_names(wave_m_RHS, self.params)
            
            self.params[j] += ep     
                     
            wave_p_RHS.calc_Psi()
            wave_p_RHS.make_padded_Psi(t)
            wave_m_RHS.calc_Psi()
            wave_m_RHS.make_padded_Psi(t)
            
            glitch_p_RHS = Glitch(wave_p_RHS, comp_num, self.Orbit)
            glitch_m_RHS = Glitch(wave_m_RHS, comp_num, self.Orbit)
            
            glitch_p_RHS.calc_TDI()
            glitch_m_RHS.calc_TDI()
            
            # take the derivatives and stuff in TDI of plus glitch
            glitch_p_RHS.TDI.A = (glitch_p_RHS.TDI.A - glitch_m_RHS.TDI.A)/2/ep
            glitch_p_RHS.TDI.E = (glitch_p_RHS.TDI.E - glitch_m_RHS.TDI.E)/2/ep
            glitch_p_RHS.TDI.T = (glitch_p_RHS.TDI.T - glitch_m_RHS.TDI.T)/2/ep

            Fisher[i][j]   = np.sum(td.get_TDI_overlap(glitch_p_LHS.TDI, glitch_p_RHS.TDI, wave_p_RHS.f_min, wave_p_RHS.f_max))
            
            del glitch_p_RHS
            del glitch_m_RHS
     
        del glitch_p_LHS
        del glitch_m_LHS

    #Fisher[wv.IDX_t0, wv.IDX_phi0] = 0.
        
    # Take advantage of the symmetry of the Fisher matrix
    for i in range(IDX_comp_id):
        for j in range(i+1, IDX_comp_id):
            Fisher[j][i] = Fisher[i][j]
               
    self.Fisher = Fisher
    
    param_vec_to_names(self.Wavelet, self.params) # reset wavelet parameters back to original values
    
    return
    

class Glitch:
    """ Glitch Class """

    def __init__(self, Wavelet, comp_id, Orbit):
        self.Wavelet = Wavelet
        self.comp_id = COMP_ID_LS[comp_id]
        self.Orbit   = Orbit

        self.params = np.zeros(6)

        # create the dimensionless version of the parameters
        self.params[wv.IDX_lnA]   = np.log(Wavelet.A)
        self.params[wv.IDX_f0]    = Wavelet.f0/mHz
        self.params[wv.IDX_t0]    = Wavelet.t0/Week
        self.params[wv.IDX_tau]   = Wavelet.tau/Week
        self.params[wv.IDX_phi0]  = Wavelet.phi0
        self.params[IDX_comp_id]  = comp_id
        
        self.f_min = Wavelet.f_min
        self.f_max = Wavelet.f_max


    # methods
    calc_TDI  = calc_TDI
    set_props = set_glitch_properties
    calc_snr  = calc_glitch_snr
    adjust_snr = adjust_to_target_snr
    calc_Fish = calc_Fisher
    param_vec_to_names = param_vec_to_names
	
	
	
	
	