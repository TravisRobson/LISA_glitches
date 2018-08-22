import numpy as np

import LISA as l
import TDI as td

import Wavelet as wv
import copy

IDX_cost  = 5
IDX_phi   = 6
IDX_psi   = 7
IDX_ellip = 8

mHz = 1.0e-3
Week = 3600.*24*7

def calc_xi(self):
	""" Construct the wave variable for times associated with GW """
	
	k_dot_x = self.k[0]*self.x[0,:,:] + self.k[1]*self.x[1,:,:] + self.k[2]*self.x[2,:,:]

	self.xi = self.t.reshape((1,self.N)) - k_dot_x/l.Clight

	return
	
def calc_k_dot_r(self):
	""" Perform the dot product between unit-direction vector and the S/C separation vectors """
	
	self.k_dot_r = self.k[0]*self.rij[0,:,:,:] + self.k[1]*self.rij[1,:,:,:] + self.k[2]*self.rij[2,:,:,:]
	
	return	

def construct_detector_tensor(self):
	""" Caclulate the detector tensor for LISA's response to GWs """
	
	self.x = self.Orbit.get_orbit(self.t) # get the S/C for all relevant times
	
	calc_xi(self) # Calculate the wave variables
	
	self.rij = self.Orbit.get_seps(self.t, self.x) # calculate separation vectors between spacecraft
	
	# calculate the outer product between unit separation vectors
	self.r_outer_r = 0.5*self.Orbit.get_rOr(self.rij)
	 
	calc_k_dot_r(self)

	return
	
def construct_basis_tensors(self):
	""" Calculate the GW basis tensors pieces """
	
	u = np.array([self.cth*self.cphi, self.cth*self.sphi, -self.sth])
	v = np.array([self.sphi, -self.cphi, 0.0])

	ep = np.outer(u,u) - np.outer(v,v)
	ec = np.outer(u,v) + np.outer(v,u)
	
	self.ep = self.c2psi*ep - self.s2psi*ec
	self.ec = self.s2psi*ep + self.c2psi*ec
		
	return
	
def calc_Hcp_ij(self):
    """ Calculate the integrated GW tensor """
    
    hp = self.Wavelet
    ellip = self.ellip

    hp0_delayed = hp.get_Psi(self.xi[0] + self.Orbit.L/l.Clight)
    hp0         = hp.get_Psi(self.xi[0])
    hc0_delayed = ellip*hp0_delayed 
    hc0         = ellip*hp0   

    hp1_delayed = hp.get_Psi(self.xi[1] + self.Orbit.L/l.Clight)
    hp1         = hp.get_Psi(self.xi[1])
    hc1_delayed = ellip*hp1_delayed 
    hc1         = ellip*hp1         

    hp2_delayed = hp.get_Psi(self.xi[2] + self.Orbit.L/l.Clight)
    hp2         = hp.get_Psi(self.xi[2])
    hc2_delayed = ellip*hp2_delayed 
    hc2         = ellip*hp2         

    self.Hpij[0,1] = hp1_delayed - hp0
    self.Hpij[1,0] = hp0_delayed - hp1

    self.Hpij[0,2] = hp2_delayed - hp0
    self.Hpij[2,0] = hp0_delayed - hp2

    self.Hpij[1,2] = hp2_delayed - hp1
    self.Hpij[2,1] = hp1_delayed - hp2

    # cross-polarization
    self.Hcij[0,1] = hc1_delayed - hc0
    self.Hcij[1,0] = hc0_delayed - hc1

    self.Hcij[0,2] = hc2_delayed - hc0
    self.Hcij[2,0] = hc0_delayed - hc2

    self.Hcij[1,2] = hc2_delayed - hc1
    self.Hcij[2,1] = hc1_delayed - hc2

    dft = 2*np.pi*self.q/self.Orbit.Tobs*self.t
    phase_factor = np.exp(-1j*dft)
    self.Hpij *= phase_factor
    self.Hcij *= phase_factor

    return
	
def calc_Hij(self):
	""" Construct the GW tensor """
	
	self.Hij = np.zeros((3,3, 3,3, self.N),dtype=np.complex_)
	
	self.Hij[0,0,:,:,:] = self.Hpij*self.ep[0,0] + self.Hcij*self.ec[0,0]
	self.Hij[0,1,:,:,:] = self.Hpij*self.ep[0,1] + self.Hcij*self.ec[0,1]
	self.Hij[0,2,:,:,:] = self.Hpij*self.ep[0,2] + self.Hcij*self.ec[0,2]
	
	self.Hij[1,0,:,:,:] = self.Hpij*self.ep[1,0] + self.Hcij*self.ec[1,0] 
	self.Hij[1,1,:,:,:] = self.Hpij*self.ep[1,1] + self.Hcij*self.ec[1,1]
	self.Hij[1,2,:,:,:] = self.Hpij*self.ep[1,2] + self.Hcij*self.ec[1,2]
	
	self.Hij[2,0,:,:,:] = self.Hpij*self.ep[2,0] + self.Hcij*self.ec[2,0] 
	self.Hij[2,1,:,:,:] = self.Hpij*self.ep[2,1] + self.Hcij*self.ec[2,1] 
	self.Hij[2,2,:,:,:] = self.Hpij*self.ep[2,2] + self.Hcij*self.ec[2,2]
		
	return
	
def get_l(GW_glitch,i,j):
	""" get the strain induced by the GW """
		   
	temp = np.einsum('nmk,nmk->k', GW_glitch.r_outer_r[:,:,i,j,:], GW_glitch.Hij[:,:,i,j,:])
		   
	return temp
	
def contract_tenors(self):
	""" Contract the detector tensor with the GW tensor """

	self.r_outer_r[:,:,0,1,:] = self.r_outer_r[:,:,0,1,:]/(1. - self.k_dot_r[0,1,:])
	self.r_outer_r[:,:,0,2,:] = self.r_outer_r[:,:,0,2,:]/(1. - self.k_dot_r[0,2,:])
	
	self.r_outer_r[:,:,1,0,:] = self.r_outer_r[:,:,1,0,:]/(1. - self.k_dot_r[1,0,:])
	self.r_outer_r[:,:,1,2,:] = self.r_outer_r[:,:,1,2,:]/(1. - self.k_dot_r[1,2,:])
	
	self.r_outer_r[:,:,2,0,:] = self.r_outer_r[:,:,2,0,:]/(1. - self.k_dot_r[2,0,:])
	self.r_outer_r[:,:,2,1,:] = self.r_outer_r[:,:,2,1,:]/(1. - self.k_dot_r[2,1,:])

	self.delta_l = np.zeros((3,3,self.N),dtype=np.complex_)
 
	self.delta_l[0,1,:] = get_l(self,0,1)
	self.delta_l[1,0,:] = get_l(self,1,0)
	
	self.delta_l[0,2,:] = get_l(self,0,2)
	self.delta_l[2,0,:] = get_l(self,2,0)
	
	self.delta_l[1,2,:] = get_l(self,1,2)
	self.delta_l[2,1,:] = get_l(self,2,1)
  
	return
	
def calculate_strain(self):
	""" Calculate the GW strain due to Sine-Gaussian burst """
	
	self.Wavelet = wv.Wavelet(self.A, self.f0, self.tau, self.t0, self.phi0, self.Orbit)
	
	self.Hpij = np.zeros((3,3,self.N),dtype=np.complex_)
	self.Hcij = np.zeros((3,3,self.N),dtype=np.complex_)
	
	calc_Hcp_ij(self)	

	construct_basis_tensors(self)

	calc_Hij(self)
	
	contract_tenors(self)

	return
	
def calc_k(self):
	""" Calculate the unit-direction vector pointing towards the source """
	
	self.k = -np.array([self.sth*self.cphi, self.sth*self.sphi, self.cth])

	return
	
def set_t(self):
	""" Set the times associated with the TDI footprint """
	
	self.N = 2**9
	dt = self.Orbit.Tobs/self.N
	self.t = np.linspace(0, self.N-1, self.N)*self.Orbit.Tobs/self.N
	
	return
	
def construct_TDI(self, Orbit):
	""" construct the TDI channels for the GW """
	
	#self.make_padded_delta_l(t)

	p12 = td.Phase(1,2, self.t, self.delta_l[0,1])
	p21 = td.Phase(2,1, self.t, self.delta_l[1,0])

	p13 = td.Phase(1,3, self.t, self.delta_l[0,2])
	p31 = td.Phase(3,1, self.t, self.delta_l[2,0])

	p23 = td.Phase(2,3, self.t, self.delta_l[1,2])
	p32 = td.Phase(3,2, self.t, self.delta_l[2,1])
   
	p12.FT_phase_FAST(Orbit)
	p21.FT_phase_FAST(Orbit)
	p13.FT_phase_FAST(Orbit)
	p31.FT_phase_FAST(Orbit)
	p23.FT_phase_FAST(Orbit)
	p32.FT_phase_FAST(Orbit)
	
	p12.phi_FT = np.fft.fftshift(p12.phi_FT)
	p21.phi_FT = np.fft.fftshift(p21.phi_FT)
	p13.phi_FT = np.fft.fftshift(p13.phi_FT)
	p31.phi_FT = np.fft.fftshift(p31.phi_FT)
	p23.phi_FT = np.fft.fftshift(p23.phi_FT)
	p32.phi_FT = np.fft.fftshift(p32.phi_FT)
	
	N_lo = (self.q - self.N/2)
	N_hi = (self.q + self.N/2)
	p12.freqs = np.arange(N_lo, N_hi, 1)/Orbit.Tobs

	tdi_GW = td.TDI(p12, p21, p13, p31, p23, p32, Orbit)
	
	tdi_GW.f_min = p12.freqs[0]
	tdi_GW.f_max = p12.freqs[-1]
	
	return tdi_GW
	
def calc_gw_snr(self, X_flag=None):
    """ Calculate the SNR of the GW burst """
    
    #  Todo: create flag such that X channel or AET is an option    
    snr_list = np.sum(self.TDI.get_TDI_snr(self.TDI.f_min, self.TDI.f_max, X_flag))

    self.SNR = np.sqrt(snr_list)

    return 
    
def adjust_to_target_snr(self, target, X_flag=None):
    """ Adjust the SNR and amplitude to hit target SNR, and its TDI data """
    
    if (self.SNR != target):
        amp = target/self.SNR
    
        self.A *= amp
        self.paramsND[wv.IDX_lnA] = np.log(self.A) # this might already update automatically with python
    
        self.construct_detector_tensor()
        self.calculate_strain()
        self.TDI = self.construct_TDI(self.Orbit)

        self.calc_snr(X_flag)
    
    return
	
def calculate_Fisher(self, X_flag=None):
    """ Calculate the Fisher matrix for a GW burst """

    orb = self.Orbit
    t = np.arange(0.0, orb.Tobs, orb.dt)

    epsilon = 1.0e-6
    NP = IDX_ellip+1 ## wv.IDX_phi0+1
    Fisher = np.zeros((NP, NP))

    for i in range(NP):
        self.paramsND[i] += epsilon
        gw_p_LHS = Burst(copy.deepcopy(self.paramsND), orb)
        gw_p_LHS.construct_detector_tensor()
        gw_p_LHS.calculate_strain()
        gw_p_LHS.TDI = gw_p_LHS.construct_TDI(orb)

        self.paramsND[i] -= 2*epsilon
        gw_m_LHS = Burst(copy.deepcopy(self.paramsND), orb)
        gw_m_LHS.construct_detector_tensor()
        gw_m_LHS.calculate_strain()
        gw_m_LHS.TDI = gw_m_LHS.construct_TDI(orb)

        self.paramsND[i] += epsilon

        if (X_flag==None):
            gw_p_LHS.TDI.A = (gw_p_LHS.TDI.A - gw_m_LHS.TDI.A)/(2*epsilon)
            gw_p_LHS.TDI.E = (gw_p_LHS.TDI.E - gw_m_LHS.TDI.E)/(2*epsilon)
            gw_p_LHS.TDI.T = (gw_p_LHS.TDI.T - gw_m_LHS.TDI.T)/(2*epsilon)
        else:
            gw_p_LHS.TDI.X = (gw_p_LHS.TDI.X - gw_m_LHS.TDI.X)/(2*epsilon)

        for j in range(i, NP):
            self.paramsND[j] += epsilon
            gw_p_RHS = Burst(copy.deepcopy(self.paramsND), orb)
            gw_p_RHS.construct_detector_tensor()
            gw_p_RHS.calculate_strain()
            gw_p_RHS.TDI = gw_p_RHS.construct_TDI(orb)

            self.paramsND[j] -= 2*epsilon
            gw_m_RHS = Burst(copy.deepcopy(self.paramsND), orb)
            gw_m_RHS.construct_detector_tensor()
            gw_m_RHS.calculate_strain()
            gw_m_RHS.TDI = gw_m_RHS.construct_TDI(orb)

            self.paramsND[j] += epsilon    

            if (X_flag==None):
                gw_p_RHS.TDI.A = (gw_p_RHS.TDI.A - gw_m_RHS.TDI.A)/(2*epsilon)
                gw_p_RHS.TDI.E = (gw_p_RHS.TDI.E - gw_m_RHS.TDI.E)/(2*epsilon)
                gw_p_RHS.TDI.T = (gw_p_RHS.TDI.T - gw_m_RHS.TDI.T)/(2*epsilon)          
            else:
                gw_p_RHS.TDI.X = (gw_p_RHS.TDI.X - gw_m_RHS.TDI.X)/(2*epsilon) 

            Fisher[i][j] = np.sum(td.get_TDI_overlap(gw_p_LHS.TDI, gw_p_RHS.TDI, gw_p_LHS.TDI.f_min, gw_p_LHS.TDI.f_max, X_flag))

    del gw_p_RHS
    del gw_m_RHS

    del gw_p_LHS
    del gw_m_LHS

    for i in range(NP):
        for j in range(i+1, NP):
            Fisher[j][i] = Fisher[i][j]

    self.Fisher = Fisher

    return


class Burst:
    """ Gravitational wave Sine-Gaussian """

    def __init__(self, paramsND, Orbit):
        self.paramsND = paramsND
        self.Orbit = Orbit

        # extract the parameters
        self.A    = np.exp(paramsND[wv.IDX_lnA])
        self.f0   = paramsND[wv.IDX_f0]*mHz
        self.t0   = paramsND[wv.IDX_t0]*Week
        self.tau  = paramsND[wv.IDX_tau]*Week
        self.phi0 = paramsND[wv.IDX_phi0]

        # set sky and polarization angles, evaluate respective trig functions
        #self.theta = theta
        #         self.cth = np.cos(self.theta)
        #         self.sth = np.sin(self.theta)
        self.cth = paramsND[IDX_cost]
        self.sth = np.sqrt(1. - self.cth**2)

        #self.phi  = phi
        self.phi  = paramsND[IDX_phi]
        self.sphi = np.sin(self.phi)
        self.cphi = np.cos(self.phi)

        # self.psi   = psi
        self.psi   = paramsND[IDX_psi]
        self.s2psi = np.sin(2.0*self.psi)
        self.c2psi = np.cos(2.0*self.psi)

        self.ellip = paramsND[IDX_ellip]

        calc_k(self) # construct direction-to-source vector

        ########################### This time shit is pretty hacky ATM ###################
        # Todo: Need to make smarter time bounds
        # need to get a common times associated with sampling the wavelets
        # added 3L to lower bound due to temporal footprint of TDI channels

        # First choose them like we have a just a wavelet
        # Todo: smarter choices for these bounds
        #			doesn't seem good when tau ~ 1/f0
        self.t_min = self.t0 - 3.*self.tau
        self.t_max = self.t0 + 3.*self.tau

        # adjust to times LISA actually sampled
        self.t_min = int(self.t_min/Orbit.dt)*Orbit.dt
        self.t_max = int(self.t_max/Orbit.dt)*Orbit.dt

        # incase burst is shorter than detector sampling rate
        if (self.t_min == self.t_max):
            self.t_max = self.t_min + 2*Orbit.dt

        # Todo: how to make sure time doesn't exceed LISA observation times
        if (self.t_min < 0.0): # ensure time is positive
            self.t_min = 0.0

        # now we must adjust for LISA considerations	
        self.t_min = self.t_min - 3.*Orbit.L/l.Clight

        # adjust to make it a sampled time
        self.t_min = int(self.t_min/Orbit.dt)*Orbit.dt
        self.t_max = int(self.t_max/Orbit.dt)*Orbit.dt

        # Todo: how to make sure time doesn't exceed LISA observation times
        if (self.t_min < 0.0): # ensure time is positive
            self.t_min = 0.0
        if (self.t_max > self.Orbit.Tobs-self.Orbit.dt):
            self.t_max = self.Orbit.Tobs-self.Orbit.dt
        ########################### 

        # Todo: smarter choices for these bounds, probably SNR dependent
        self.f_min = self.f0 - 2.5/self.tau
        self.f_max = self.f0 + 2.5/self.tau
        self.f_min = np.floor(self.f_min*self.Orbit.Tobs)/self.Orbit.Tobs
        self.f_min = np.ceil(self.f_max*self.Orbit.Tobs)/self.Orbit.Tobs
        
        # catch to make sure negative frequencies aren't asked for
        if (self.f_min < 0.0): 
            self.f_min = 0.0
        # ensure max frequency is not greater than the Nyquist frequency
        if (self.f_max > Orbit.f_ny):
            self.f_max = Orbit.f_ny

        set_t(self)
        self.q = np.floor(self.f0*self.Orbit.Tobs) # carrier frequency bin

    # Methods
    construct_detector_tensor = construct_detector_tensor	
    calculate_strain = calculate_strain
    construct_TDI = construct_TDI
    calc_Fish = calculate_Fisher
    calc_snr = calc_gw_snr
    adjust_snr = adjust_to_target_snr