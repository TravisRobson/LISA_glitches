import numpy as np

import Wavelet as wv

import LISA as l
import TDI as td

from copy import deepcopy

IDX_comp_id = 5

Week = 7.*24*3600
mHz = 1.0e-3

    
def calc_TDI(self):
    """ Analytically Calculate the Glithch FT """
    
    Wave  = self.Wavelet
    Orbit = self.Orbit    

    t = np.arange(0.0, Orbit.Tobs, Orbit.dt) # Todo: Don't need Orbit, its in Wavelet
    N = len(t)
    
    f_min = Wave.f_min
    f_max = Wave.f_max

    N_lo = int(np.floor(f_min*Orbit.Tobs))
    N_hi = int(np.ceil (f_max*Orbit.Tobs))
    
    p12 = td.Phase(1,2, np.zeros(1), np.zeros(1))
    p21 = td.Phase(2,1, np.zeros(1), np.zeros(1))
    p13 = td.Phase(1,3, np.zeros(1), np.zeros(1))
    p31 = td.Phase(3,1, np.zeros(1), np.zeros(1))
    p23 = td.Phase(2,3, np.zeros(1), np.zeros(1))
    p32 = td.Phase(3,2, np.zeros(1), np.zeros(1))
    
    p12.phi_FT = np.zeros(1, dtype=np.complex_)
    p21.phi_FT = np.zeros(1, dtype=np.complex_)
    p13.phi_FT = np.zeros(1, dtype=np.complex_)
    p31.phi_FT = np.zeros(1, dtype=np.complex_)
    p23.phi_FT = np.zeros(1, dtype=np.complex_)
    p32.phi_FT = np.zeros(1, dtype=np.complex_)
        
    if (self.comp_id == 0): # Optical Path glitch, 1->2
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs
        p12.phi_FT = wv.get_Psi_FT(p12.freqs, Wave)
        
    elif (self.comp_id == 1): # Optical Path glitch, 2->1
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs # Yes, p12... (B/C TDI construction)
        p21.phi_FT = wv.get_Psi_FT(p12.freqs, Wave)
        
    elif (self.comp_id == 2): # Optical Path glitch, 1->3
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs # Yes, p12... (B/C TDI construction)
        p13.phi_FT = wv.get_Psi_FT(p12.freqs, Wave)
        
    elif (self.comp_id == 3): # Optical Path glitch, 3->1
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs # Yes, p12... (B/C TDI construction)
        p31.phi_FT = wv.get_Psi_FT(p12.freqs, Wave)
        
    elif (self.comp_id == 4): # Optical Path glitch, 2->3
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs # Yes, p12... (B/C TDI construction)
        p23.phi_FT = wv.get_Psi_FT(p12.freqs, Wave)

    elif (self.comp_id == 5): # Optical Path glitch, 3->2
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs # Yes, p12... (B/C TDI construction)
        p32.phi_FT = wv.get_Psi_FT(p12.freqs, Wave)
        
    elif (self.comp_id == 6): # Acceleration glitch, 1->2
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs
        
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = wv.Wavelet(Wave.A, Wave.f0, Wave.tau, \
                    Wave.t0 + Orbit.L/l.Clight, Wave.phi0 - Wave.f0/l.fstar, Orbit)
                    
        p12.phi_FT = -wv.get_Psi_FT(p12.freqs, Wave) 
        p21.phi_FT =  wv.get_Psi_FT(p12.freqs, wave_temp) 

    elif (self.comp_id == 7): # Acceleration glitch, 2->1
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs
        
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = wv.Wavelet(Wave.A, Wave.f0, Wave.tau, \
                    Wave.t0 + Orbit.L/l.Clight, Wave.phi0 - Wave.f0/l.fstar, Orbit)
                    
        p21.phi_FT = -wv.get_Psi_FT(p12.freqs, Wave) 
        p12.phi_FT =  wv.get_Psi_FT(p12.freqs, wave_temp)         

    elif (self.comp_id == 8): # Acceleration glitch, 1->3
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs
        
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = wv.Wavelet(Wave.A, Wave.f0, Wave.tau, \
                    Wave.t0 + Orbit.L/l.Clight, Wave.phi0 - Wave.f0/l.fstar, Orbit)
                    
        p13.phi_FT = -wv.get_Psi_FT(p12.freqs, Wave) 
        p31.phi_FT =  wv.get_Psi_FT(p12.freqs, wave_temp) 

    elif (self.comp_id == 9): # Acceleration glitch, 3->1
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs
        
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = wv.Wavelet(Wave.A, Wave.f0, Wave.tau, \
                    Wave.t0 + Orbit.L/l.Clight, Wave.phi0 - Wave.f0/l.fstar, Orbit)
                    
        p31.phi_FT = -wv.get_Psi_FT(p12.freqs, Wave) 
        p13.phi_FT =  wv.get_Psi_FT(p12.freqs, wave_temp) 

    elif (self.comp_id == 10): # Acceleration glitch, 2->3
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs
        
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = wv.Wavelet(Wave.A, Wave.f0, Wave.tau, \
                    Wave.t0 + Orbit.L/l.Clight, Wave.phi0 - Wave.f0/l.fstar, Orbit)
                    
        p23.phi_FT = -wv.get_Psi_FT(p12.freqs, Wave) 
        p32.phi_FT =  wv.get_Psi_FT(p12.freqs, wave_temp)  

    elif (self.comp_id == 11): # Acceleration glitch, 3->2
        p12.freqs  = np.arange(N_lo, N_hi, 1)/Orbit.Tobs
        
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = wv.Wavelet(Wave.A, Wave.f0, Wave.tau, \
                    Wave.t0 + Orbit.L/l.Clight, Wave.phi0 - Wave.f0/l.fstar, Orbit)
                    
        p32.phi_FT = -wv.get_Psi_FT(p12.freqs, Wave) 
        p23.phi_FT =  wv.get_Psi_FT(p12.freqs, wave_temp)  
    
    p12.freqs = np.arange(N_lo, N_hi, 1)/Orbit.Tobs
    
    self.TDI = td.TDI(p12, p21, p13, p31, p23, p32, Orbit)
    self.TDI.f_min = Wave.f_min
    self.TDI.f_max = Wave.f_max    
    
    return
    
def calc_glitch_snr(self, X_flag=None):
    """ Calculate the SNR of the glitch """
    
    #  Todo: create flag such that X channel or AET is an option    
    snr_list = np.sum(self.TDI.get_TDI_snr(self.f_min, self.f_max, X_flag))

    self.SNR = np.sqrt(snr_list)

    return 
    
def adjust_to_target_snr(self, target, X_flag=None):
    """ Adjust the SNR and amplitude to hit target SNR, and its TDI data """
    
    amp = target/self.SNR
    
    self.Wavelet.A *= amp
    self.paramsND[wv.IDX_lnA] = np.log(self.Wavelet.A) # this might already update automatically with python
    
    self.calc_TDI()
    self.calc_snr(X_flag)
    
    return
    
def calc_Fisher(self, X_flag=None):
    """ Calculate the Fisher matrix for the glitch object """
    
    ep = 1.0e-9
    
    Fisher = np.zeros((IDX_comp_id, IDX_comp_id))
           
    for i in range(IDX_comp_id):
    
        self.paramsND[i] += ep
        
        glitch_p_LHS = Glitch(deepcopy(self.paramsND), self.Orbit, self.comp_id)
        
        self.paramsND[i] -= 2*ep
        
        glitch_m_LHS = Glitch(deepcopy(self.paramsND), self.Orbit, self.comp_id)
               
        self.paramsND[i] += ep    # to return parameter to initial value           
        
        glitch_p_LHS.calc_TDI()
        glitch_m_LHS.calc_TDI()
        
        # take the derivatives and stuff in TDI of plus glitch
        if (X_flag==None):
            glitch_p_LHS.TDI.A = (glitch_p_LHS.TDI.A - glitch_m_LHS.TDI.A)/(2*ep)
            glitch_p_LHS.TDI.E = (glitch_p_LHS.TDI.E - glitch_m_LHS.TDI.E)/(2*ep)
            glitch_p_LHS.TDI.T = (glitch_p_LHS.TDI.T - glitch_m_LHS.TDI.T)/(2*ep)  
        else:
            glitch_p_LHS.TDI.X = (glitch_p_LHS.TDI.X - glitch_m_LHS.TDI.X)/(2*ep)

        for j in range(i, IDX_comp_id):
            
            # RHS 
            self.paramsND[j] += ep
            
            glitch_p_RHS = Glitch(deepcopy(self.paramsND), self.Orbit, self.comp_id)

            self.paramsND[j] -= 2*ep

            glitch_m_RHS = Glitch(deepcopy(self.paramsND), self.Orbit, self.comp_id)
            
            self.paramsND[j] += ep     
            
            glitch_p_RHS.calc_TDI()
            glitch_m_RHS.calc_TDI()
            
            # take the derivatives and stuff in TDI of plus glitch
            if (X_flag==None):
                glitch_p_RHS.TDI.A = (glitch_p_RHS.TDI.A - glitch_m_RHS.TDI.A)/(2*ep)
                glitch_p_RHS.TDI.E = (glitch_p_RHS.TDI.E - glitch_m_RHS.TDI.E)/(2*ep)
                glitch_p_RHS.TDI.T = (glitch_p_RHS.TDI.T - glitch_m_RHS.TDI.T)/(2*ep)
            else:
                glitch_p_RHS.TDI.X = (glitch_p_RHS.TDI.X - glitch_m_RHS.TDI.X)/(2*ep)
            
            f_min = np.max([glitch_p_RHS.f_min, glitch_m_RHS.f_min, glitch_p_LHS.f_min, glitch_m_LHS.f_min])
            f_max = np.min([glitch_p_RHS.f_max, glitch_m_RHS.f_max, glitch_p_LHS.f_max, glitch_m_LHS.f_max])

            Fisher[i][j]   = np.sum(td.get_TDI_overlap(glitch_p_LHS.TDI, glitch_p_RHS.TDI, glitch_p_RHS.f_min, glitch_p_RHS.f_max, X_flag))
            
            del glitch_p_RHS
            del glitch_m_RHS
     
        del glitch_p_LHS
        del glitch_m_LHS
        
    # Take advantage of the symmetry of the Fisher matrix
    for i in range(IDX_comp_id):
        for j in range(i+1, IDX_comp_id):
            Fisher[j][i] = Fisher[i][j]
               
    self.Fisher = Fisher
        
    return
    

class Glitch:
    """ Glitch Class """

    def __init__(self, paramsND, Orbit, comp_id):
        self.paramsND = paramsND
        self.Orbit   = Orbit
        self.comp_id = int(comp_id)
        
        self.A = np.exp(self.paramsND[wv.IDX_lnA])
        self.f0 = self.paramsND[wv.IDX_f0]*mHz
        self.t0 = self.paramsND[wv.IDX_t0]*Week
        self.tau = self.paramsND[wv.IDX_tau]*Week
        self.phi0 = self.paramsND[wv.IDX_phi0]

        self.Wavelet = wv.Wavelet(self.A, self.f0, self.tau, self.t0, self.phi0, self.Orbit)
        
        self.f_min = self.Wavelet.f_min
        self.f_max = self.Wavelet.f_max


    # methods
    calc_TDI  = calc_TDI
    calc_snr  = calc_glitch_snr
    adjust_snr = adjust_to_target_snr
    calc_Fish = calc_Fisher
	
	
	
	
	