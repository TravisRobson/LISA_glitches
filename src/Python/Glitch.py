import numpy as np
import scipy.special as ss

import LISA as l

def evaluate_wavelet(self, t):
	"""
	Return the value of the wavelet at specified t
	"""	
		
	result = np.cos(2.*np.pi*self.f0*(t-self.t0)+self.phi0)*np.exp(-(t-self.t0)**2./(self.tau**2.))*self.A
# 		
# 	if (self.domain == 'frequency'):
# 		raise ValueError("'freqeuency' domain not yet supported. Please use 'time'")
		
	return result
	
	
def set_indexes(self, Orbit):

	self.idx_left  = int(self.t_min/Orbit.dt)
	self.idx_right = int(self.t_max/Orbit.dt)
	
	if (self.idx_left == self.idx_right):
		self.idx_left  = self.idx_left  - 1
		self.idx_right = self.idx_right + 1
	
	# must reset these values to detector times
	self.t_min     = self.idx_left*Orbit.dt
	self.t_max     = self.idx_right*Orbit.dt
	
	return
	
def set_t_h(self, Orbit):
	
	self.t = np.arange(self.t_min, self.t_max, Orbit.dt)
	self.h = self.get_h(self.t)
	
	return
	
def make_padded_h(self, t):
	
	left_pad  = np.argwhere(t == self.t[0] ).flatten()[0]
	right_pad = len(t) - 1 - np.argwhere(t == self.t[-1]).flatten()[0]

	self.h_padded = np.pad(self.h, (left_pad,right_pad), 'constant')
	
	return

def get_integrated_wavelet(self, xi):
	
	a1 = (xi - self.t0)/self.tau
	a2 = np.pi*self.f0*self.tau
	
	arg1 = a1 + 1.0j*a2
	
	temp1 = ss.erf(arg1)*np.exp(-1.0j*self.phi0)

	result = temp1.real
# 	
	result = result.real
	result = result*np.sqrt(np.pi)*0.5*self.A*self.tau*np.exp(-(self.f0*np.pi*self.tau)**2)
	
	return result

class Wavelet:
	kind = 'Morlet-Gabor Wavelet'
	
	def __init__(self, A, f0, tau, t0, phi0):
		self.A      = A
		self.f0     = f0
		self.tau    = tau
		self.t0     = t0
		self.phi0   = phi0
		self.Q  = 2.*np.pi*self.f0*self.tau
		
		# will be modified
		self.t_min = self.t0 - 4.*self.tau
		self.t_max = self.t0 + 4.*self.tau
		
		# TODO: NEEDS CAREFUL CONSIDERATION FOR WHAT THESE BOUNDS SHOULD BE
		self.f_min = self.f0 - 3./self.tau
		self.f_max = self.f0 + 3./self.tau
						
		
	# methods
	get_h = evaluate_wavelet
	set_indexes = set_indexes
	set_t_h = set_t_h
	make_padded_h = make_padded_h
	
	get_integrated_wavelet = get_integrated_wavelet # more for the GW model
	
	
def FT_phase(self, Orbit):
	
	phi = np.fft.rfft(self.phi)
	self.phi_FT = phi[0:len(self.t)/2]*Orbit.Tobs/len(self.t)
	self.freqs  = np.fft.fftfreq(self.t.shape[-1])[0:len(self.t)/2]/Orbit.dt
	
	return
	

class Phase:
	
	def __init__(self, i, j, t, phi):
		self.i      = i 	 # emitting S/C
		self.j      = j	     # receiving S/C
		self.t      = t      # time array
		self.phi    = phi    # phase array
		
	FT_phase = FT_phase				


	
class TDI:
	
	def __init__(self, phi12, phi21, phi13, phi31, phi23, phi32, Orbit):
		"""
		Take the phase comparisons, FFT them, and the construct
			the TDI data channels
		
		TODO: Currently, we are assuming L_{ij} = L!!!!
		"""
		self.freqs = phi12.freqs # TODO: need some logic to handle if another phase is being used instead
		
		fonfs = self.freqs/l.fstar
	
		phase1 = np.cos(fonfs)    - 1.0j*np.sin(fonfs)
		phase2 = np.cos(2.*fonfs) - 1.0j*np.sin(2.*fonfs)
		phase3 = np.cos(3.*fonfs) - 1.0j*np.sin(3.*fonfs)
		
		self.X =  (phi12.phi_FT - phi13.phi_FT)*phase3 \
				+ (phi21.phi_FT - phi31.phi_FT)*phase2 \
		 	    + (phi13.phi_FT - phi12.phi_FT)*phase1 \
		 	    + (phi31.phi_FT - phi21.phi_FT)

		self.Y =  (phi23.phi_FT - phi21.phi_FT)*phase3 \
				+ (phi32.phi_FT - phi12.phi_FT)*phase2 \
		 	    + (phi21.phi_FT - phi23.phi_FT)*phase1 \
		 	    + (phi12.phi_FT - phi32.phi_FT)
		 	    
		self.Z =  (phi31.phi_FT - phi32.phi_FT)*phase3 \
				+ (phi13.phi_FT - phi23.phi_FT)*phase2 \
		 	    + (phi32.phi_FT - phi31.phi_FT)*phase1 \
		 	    + (phi23.phi_FT - phi13.phi_FT)
		
		
		self.A = 1./3.*(2.0*self.X - self.Y - self.Z)
		self.E = 1./np.sqrt(3.)*(self.Z - self.Y)
		self.T = 1./3.*(self.X + self.Y + self.Z)


		
def set_GW_indices(self, Orbit):
	# set the corresponding indices
	self.idx_left  = int(self.t_min/Orbit.dt)
	self.idx_right = int(self.t_max/Orbit.dt)

	if (self.idx_left == self.idx_right):
		self.idx_left  = self.idx_left  - 1
		self.idx_right = self.idx_right + 1	
		
	return
	
def calc_k(self):
	self.k = -np.array([self.sth*self.cphi, self.sth*self.sphi, self.cth])

	return
	
def set_t(self, Orbit):

	self.t = np.arange(self.t_min, self.t_max, Orbit.dt)
	self.N = len(self.t)
	
	return

	
def calc_xi(self):
	
	k_dot_x = self.k[0]*self.x[0,:,:] + self.k[1]*self.x[1,:,:] + self.k[2]*self.x[2,:,:]

	self.xi = np.zeros((3, self.N)) + self.t
	self.xi = self.xi - k_dot_x/l.Clight
	
	
def calc_Hcp_ij(self):
	Hp = np.zeros((3,self.N))
	Hc = np.zeros((3,self.N))
	
	Hp[0] = self.hp_wavelet.get_integrated_wavelet(self.xi[0])
	Hp[1] = self.hp_wavelet.get_integrated_wavelet(self.xi[1])
	Hp[2] = self.hp_wavelet.get_integrated_wavelet(self.xi[2])
	
	Hc[0] = self.hc_wavelet.get_integrated_wavelet(self.xi[0])
	Hc[1] = self.hc_wavelet.get_integrated_wavelet(self.xi[1])
	Hc[2] = self.hc_wavelet.get_integrated_wavelet(self.xi[2])
	
	# plus-polarization
	self.Hpij[0,1] = Hp[1] - Hp[0]
	self.Hpij[1,0] = -self.Hpij[0,1]
	
	self.Hpij[0,2] = Hp[2] - Hp[0]
	self.Hpij[2,0] = -self.Hpij[0,2]
	
	self.Hpij[1,2] = Hp[1] - Hp[2]
	self.Hpij[2,1] = -self.Hpij[1,2]
		
	# cross-polarization
	self.Hcij[0,1] = Hc[1] - Hc[0]
	self.Hcij[1,0] = -self.Hcij[0,1]
	
	self.Hcij[0,2] = Hc[2] - Hc[0]
	self.Hcij[2,0] = -self.Hcij[0,2]
	
	self.Hcij[1,2] = Hc[1] - Hc[2]
	self.Hcij[2,1] = -self.Hcij[1,2]
	
	return
	
def construct_basis_tensors(self):
	u = np.array([self.cth*self.cphi, self.cth*self.sphi, -self.sth])
	v = np.array([self.sphi, -self.cphi, 0.0]).reshape((3,1))
	
	ep = u.reshape((1,3))*u.reshape((3,1)) - v.reshape((1,3))*v.reshape((3,1))
	ec = u.reshape((1,3))*v.reshape((3,1)) + v.reshape((1,3))*u.reshape((3,1))
	
	self.ep = self.c2psi*ep - self.s2psi*ec
	self.ec = self.s2psi*ep + self.c2psi*ec
		
	return
	
def calc_Hij(self):
	self.Hij = np.zeros((3,3, 3,3, self.N))
	
	#self.Hij[0,0,:,:,:] = self.Hpij*self.ep[0,0] + self.Hcij*self.ec[0,0]
	self.Hij[0,1,:,:,:] = self.Hpij*self.ep[0,1] + self.Hcij*self.ec[0,1]
	self.Hij[0,2,:,:,:] = self.Hpij*self.ep[0,2] + self.Hcij*self.ec[0,2]
	
	self.Hij[1,0,:,:,:] = self.Hij[0,1,:,:,:]
	#self.Hij[1,1,:,:,:] = self.Hpij*self.ep[1,1] + self.Hcij*self.ec[1,1]
	self.Hij[1,2,:,:,:] = self.Hpij*self.ep[1,2] + self.Hcij*self.ec[1,2]
	
	self.Hij[2,0,:,:,:] = self.Hij[0,2,:,:,:]
	self.Hij[2,1,:,:,:] = self.Hij[1,2,:,:,:]
	#self.Hij[2,2,:,:,:] = self.Hpij*self.ep[2,2] + self.Hcij*self.ec[2,2]
	
	return

def calc_k_dot_r(self):
	
	self.k_dot_r = np.zeros((3,3,self.N))
	
	self.k_dot_r = self.k[0]*self.rij[0,:,:,:] + self.k[1]*self.rij[1,:,:,:] + self.k[2]*self.rij[2,:,:,:]
	
	return	
	
def get_l(GW_glitch,i,j):

	temp = GW_glitch.r_outer_r[0,1,i,j,:]*GW_glitch.Hij[0,1,i,j,:] + \
		   GW_glitch.r_outer_r[0,2,i,j,:]*GW_glitch.Hij[0,2,i,j,:] + \
		   GW_glitch.r_outer_r[1,0,i,j,:]*GW_glitch.Hij[1,0,i,j,:] + \
		   GW_glitch.r_outer_r[1,2,i,j,:]*GW_glitch.Hij[1,2,i,j,:] + \
		   GW_glitch.r_outer_r[2,0,i,j,:]*GW_glitch.Hij[2,0,i,j,:] + \
		   GW_glitch.r_outer_r[2,1,i,j,:]*GW_glitch.Hij[2,1,i,j,:]	   
	

	return temp
	
def contract_tenors(self):

	#self.r_outer_r[:,:,0,0,:] = self.r_outer_r[:,:,0,0,:]/self.k_dot_r[0,0,:]
	self.r_outer_r[:,:,0,1,:] = self.r_outer_r[:,:,0,1,:]/self.k_dot_r[0,1,:]
	self.r_outer_r[:,:,0,2,:] = self.r_outer_r[:,:,0,2,:]/self.k_dot_r[0,2,:]
	
	self.r_outer_r[:,:,1,0,:] = self.r_outer_r[:,:,1,1,:]
	#self.r_outer_r[:,:,1,1,:] = self.r_outer_r[:,:,1,1,:]/self.k_dot_r[1,1,:]
	self.r_outer_r[:,:,1,2,:] = self.r_outer_r[:,:,1,2,:]/self.k_dot_r[1,2,:]
	
	self.r_outer_r[:,:,2,0,:] = self.r_outer_r[:,:,0,2,:]
	self.r_outer_r[:,:,2,1,:] = self.r_outer_r[:,:,1,2,:]
	#self.r_outer_r[:,:,2,2,:] = self.r_outer_r[:,:,2,2,:]/self.k_dot_r[2,2,:]


	self.delta_l = np.zeros((3,3,self.N))
 

	self.delta_l[0,1,:] = get_l(self,0,1)
	self.delta_l[1,0,:] = get_l(self,1,0)
	
	self.delta_l[0,2,:] = get_l(self,0,2)
	self.delta_l[2,0,:] = get_l(self,2,0)
	
	self.delta_l[1,2,:] = get_l(self,1,2)
	self.delta_l[2,1,:] = get_l(self,2,1)
	

				   
	return
	
def make_padded_delta_l(self, t):
	
	left_pad  = np.argwhere(t == self.t[0] ).flatten()[0]
	right_pad = len(t) - 1 - np.argwhere(t == self.t[-1]).flatten()[0]

	self.delta_l_padded = np.zeros((3,3,len(t)))

	self.delta_l_padded[0,1] = np.pad(self.delta_l[0,1], (left_pad,right_pad), 'constant')
	self.delta_l_padded[1,0] = np.pad(self.delta_l[1,0], (left_pad,right_pad), 'constant')
	self.delta_l_padded[0,2] = np.pad(self.delta_l[0,1], (left_pad,right_pad), 'constant')
	self.delta_l_padded[2,0] = np.pad(self.delta_l[2,0], (left_pad,right_pad), 'constant')
	self.delta_l_padded[1,2] = np.pad(self.delta_l[1,2], (left_pad,right_pad), 'constant')
	self.delta_l_padded[2,1] = np.pad(self.delta_l[2,1], (left_pad,right_pad), 'constant')
	
	return

	
class GW_glitch:

	def __init__(self, hp_wavelet, hc_wavelet, theta, phi, psi, Orbit):
		self.hp_wavelet = hp_wavelet
		self.hc_wavelet = hc_wavelet
		self.Orbit = Orbit
		
		self.theta = theta
		self.cth = np.cos(self.theta)
		self.sth = np.sin(self.theta)
		
		self.phi  = phi
		self.sphi = np.sin(self.phi)
		self.cphi = np.cos(self.phi)
		
		self.psi   = psi
		self.s2psi = np.sin(2.0*self.psi)
		self.c2psi = np.cos(2.0*self.psi)
		
		calc_k(self)
		
		# need to get a common times associated with sampling the wavelets
		self.t_min = np.min([self.hp_wavelet.t_min, self.hc_wavelet.t_min])
		self.t_max = np.max([self.hp_wavelet.t_max, self.hc_wavelet.t_max])
	
		set_GW_indices(self, self.Orbit)
		set_t(self, self.Orbit)
		
		# get the S/C for all relevant times
		self.x = np.zeros((3, 3, self.N))
		self.x = self.Orbit.fill_full_orbit(self.t, self.x)

		calc_xi(self)
		
		self.Hpij = np.zeros((3,3,self.N))
		self.Hcij = np.zeros((3,3,self.N))
		calc_Hcp_ij(self)	
		
		construct_basis_tensors(self)
	
		calc_Hij(self)
		
		# calculate separation vectors between spacecraft
		self.rij = np.zeros((3,3,3,self.N))
		self.rij = self.Orbit.fill_full_seps(self.t, self.x, self.rij)
		
		# calculate the outer product between unit separation vectors
		self.r_outer_r = 0.5*self.rij.reshape((1,3,3,3,self.N))*self.rij.reshape((3,1,3,3,self.N))
		calc_k_dot_r(self)
		
		contract_tenors(self)
		
	make_padded_delta_l = make_padded_delta_l
	
	