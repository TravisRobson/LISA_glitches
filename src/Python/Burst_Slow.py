import numpy as np

import LISA as l
import TDI as td

import Wavelet as wv

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
	
	hp0_delayed = self.hp_wavelet.get_Psi(self.xi[0] + self.Orbit.L/l.Clight)
	hp0         = self.hp_wavelet.get_Psi(self.xi[0])
	hc0_delayed = self.hc_wavelet.get_Psi(self.xi[0] + self.Orbit.L/l.Clight)
	hc0         = self.hc_wavelet.get_Psi(self.xi[0])
	
	hp1_delayed = self.hp_wavelet.get_Psi(self.xi[1] + self.Orbit.L/l.Clight)
	hp1         = self.hp_wavelet.get_Psi(self.xi[1])
	hc1_delayed = self.hc_wavelet.get_Psi(self.xi[1] + self.Orbit.L/l.Clight)
	hc1         = self.hc_wavelet.get_Psi(self.xi[1])
	
	hp2_delayed = self.hp_wavelet.get_Psi(self.xi[2] + self.Orbit.L/l.Clight)
	hp2         = self.hp_wavelet.get_Psi(self.xi[2])
	hc2_delayed = self.hc_wavelet.get_Psi(self.xi[2] + self.Orbit.L/l.Clight)
	hc2         = self.hc_wavelet.get_Psi(self.xi[2])
	
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
	
def set_t(self, Orbit):
	""" Set the times associated with the TDI footprint """
	
	self.t = np.arange(self.t_min, self.t_max, Orbit.dt)
	self.N = len(self.t)
	
	return
	
def set_t_FAST(self):
	""" Set the times associated with the TDI footprint """
	
	self.N = 2**7
	dt = self.Orbit.Tobs/self.N
	self.t = np.linspace(0, self.N-1, self.N)*self.Orbit.Tobs/self.N
	
	return
	
def make_padded_delta_l(self, t):
	""" pad the strain time series if needed """
	
	left_pad  = np.argwhere(t == self.t[0] ).flatten()[0]
	right_pad = len(t) - 1 - np.argwhere(t == self.t[-1]).flatten()[0]

	self.delta_l_padded = np.zeros((3,3,len(t)),dtype=np.complex_)

	self.delta_l_padded[0,1] = np.pad(self.delta_l[0,1], (left_pad,right_pad), 'constant')
	self.delta_l_padded[1,0] = np.pad(self.delta_l[1,0], (left_pad,right_pad), 'constant')
	
	self.delta_l_padded[0,2] = np.pad(self.delta_l[0,2], (left_pad,right_pad), 'constant')
	self.delta_l_padded[2,0] = np.pad(self.delta_l[2,0], (left_pad,right_pad), 'constant')
	
	self.delta_l_padded[1,2] = np.pad(self.delta_l[1,2], (left_pad,right_pad), 'constant')
	self.delta_l_padded[2,1] = np.pad(self.delta_l[2,1], (left_pad,right_pad), 'constant')
	
	return
	
def construct_TDI(self, t, Orbit):
	""" construct the TDI channels for the GW """
	
	self.make_padded_delta_l(t)

	p12 = td.Phase(1,2, t, self.delta_l_padded[0,1,:])
	p21 = td.Phase(2,1, t, self.delta_l_padded[1,0,:])

	p13 = td.Phase(1,3, t, self.delta_l_padded[0,2,:])
	p31 = td.Phase(3,1, t, self.delta_l_padded[2,0,:])

	p23 = td.Phase(2,3, t, self.delta_l_padded[1,2,:])
	p32 = td.Phase(3,2, t, self.delta_l_padded[2,1,:])
   
	p12.FT_phase(Orbit)
	p21.FT_phase(Orbit)
	p13.FT_phase(Orbit)
	p31.FT_phase(Orbit)
	p23.FT_phase(Orbit)
	p32.FT_phase(Orbit)

	tdi_GW = td.TDI(p12, p21, p13, p31, p23, p32, Orbit)
	
	return tdi_GW
	
	
class GW_glitch:
	""" Gravitational wave Sine-Gaussian """
	
	def __init__(self, hp_wavelet, hc_wavelet, theta, phi, psi, Orbit):
		self.hp_wavelet = hp_wavelet
		self.hc_wavelet = hc_wavelet
		self.Orbit = Orbit
		
		# set sky and polarization angles, evaluate respective trig functions
		self.theta = theta
		self.cth = np.cos(self.theta)
		self.sth = np.sin(self.theta)
		
		self.phi  = phi
		self.sphi = np.sin(self.phi)
		self.cphi = np.cos(self.phi)
		
		self.psi   = psi
		self.s2psi = np.sin(2.0*self.psi)
		self.c2psi = np.cos(2.0*self.psi)
		
		calc_k(self) # construct direction-to-source vector
		
		# Todo: Need to make smarter time bounds
		# need to get a common times associated with sampling the wavelets
		# added 3L to lower bound due to temporal footprint of TDI channels
		self.t_min = np.min([self.hp_wavelet.t_min, self.hc_wavelet.t_min]) - 3.*Orbit.L/l.Clight
		self.t_max = np.max([self.hp_wavelet.t_max, self.hc_wavelet.t_max]) 
		
		# adjust to make it a sampled time
		self.t_min = int(self.t_min/Orbit.dt)*Orbit.dt
		self.t_max = int(self.t_max/Orbit.dt)*Orbit.dt

 		# Todo: how to make sure time doesn't exceed LISA observation times
		if (self.t_min < 0.0): # ensure time is positive
			self.t_min = 0.0
		if (self.t_max > self.Orbit.Tobs-self.Orbit.dt):
			self.t_max = self.Orbit.Tobs-self.Orbit.dt
			
		self.f_min = np.min([self.hp_wavelet.f_min, self.hc_wavelet.f_min])
		self.f_max = np.min([self.hp_wavelet.f_max, self.hc_wavelet.f_max])
		
		set_t(self, self.Orbit)

	
	# Methods
	construct_detector_tensor = construct_detector_tensor	
	calculate_strain = calculate_strain
	make_padded_delta_l = make_padded_delta_l
	construct_TDI = construct_TDI
	