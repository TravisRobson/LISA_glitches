import numpy as np
import scipy.special as ss

import LISA as l

def evaluate_wavelet(self, t):
	"""
	Return the value of the Sine-Gaussian wavelet at time(s) t
	
	input:
		self - (Wavelet)
		t    - (float (array)) time(s) at which to evaluate the wavelet
		
	output:
		result - (float (array)) the wavelet evaluated
	"""	
		
	arg1 = 2.*np.pi*self.f0*(t-self.t0) + self.phi0
	arg2 = (t - self.t0)**2./(self.tau**2.)
	
	result = self.A*np.cos(arg1)*np.exp(-arg2)

	return result
	
def calc_Psi(self):
	"""
	In the SSB frame, calculate relevant times and sample the wavelet
	
	input:
		self - (Wavelet)
	"""
	self.t = np.arange(self.t_min, self.t_max, self.Orbit.dt)
	self.Psi = self.get_Psi(self.t)
	
	return
	
def make_padded_Psi(self, t):
	"""
	Pad the wavelet's internal Psi to match shape of the provided time series
	
	input: 
		self - (Wavelet)
		t 	 - (float (Array)) time series for Psi_padded to match
	"""
	
	# find where input time array bounds matches the wavelet's native time array
	left_pad  = np.argwhere(t == self.t[0] ).flatten()[0]
	right_pad = len(t) - 1 - np.argwhere(t == self.t[-1]).flatten()[0]

	self.Psi_padded = np.pad(self.Psi, (left_pad,right_pad), 'constant')
#	self.t_full = np.copy(t)
	
	return

def get_integrated_wavelet(self, t):
	"""
	Evaluate the integral of a Sine-Gaussian from 0 to time(s) t
	
	input:
		self - (Wavelet) 
		t    - (float (array)) time(s) at which to evaluate the wavelet
		
	output:
		result - (float (array)) result of the integration
	"""
	# don't waste evaluating a non-existant wavelet
	if (self.A == 0):
		return np.zeros(len(xi))
	else:
		alpha = (t - self.t0)/self.tau
		beta  = np.pi*self.f0*self.tau
	
		arg = alpha + 1.0j*beta
		
		phase = np.exp(-1.0j*self.phi0)
	
		term1 = phase.real*np.exp(-beta**2.0)
		
		huh = ss.erfcx(arg)
		
		# if the real argument of erfcx is too large, it blows up, 
		#		do the straightforward way
		if (len(np.isnan(huh)) != 0):
			result = ss.erf(arg)*phase
			result = np.sqrt(np.pi)*0.5*self.A*self.tau*result.real
			return result
		
		term2 = huh*np.exp(-arg**2-beta**2)*phase
		term2 = term2.real
		
		result = np.sqrt(np.pi)*0.5*self.A*self.tau*(term1 - term2)
	
	return result


class Wavelet:
	"""
	Sine-Gaussian wavelet class
	"""
	kind = 'Sine-Gaussian'
	
	def __init__(self, A, f0, tau, t0, phi0, Orbit):
			
		self.A      = A
		self.f0     = f0
		self.tau    = tau
		self.t0     = t0
		self.phi0   = phi0
		self.Q      = np.pi*self.f0*self.tau
		self.Orbit  = Orbit
		
		# Todo: smarter choices for these bounds
		#			doesn't seem good when tau ~ 1/f0
		self.t_min = self.t0 - 3.*self.tau
		self.t_max = self.t0 + 3.*self.tau
		
		# adjust to times LISA actually sampled
		self.t_min = int(self.t_min/Orbit.dt)*Orbit.dt
		self.t_max = int(self.t_max/Orbit.dt)*Orbit.dt
		
		# incase burst is shorter than detector sampling rate
		if (self.t_min == self.t_max):
			self.t_max = self.t_min + Orbit.dt
		
		# Todo: how to make sure time doesn't exceed LISA observation times
		if (self.t_min < 0.0): # ensure time is positive
			self.t_min = 0.0
		
		# Todo: smarter choices for these bounds, probably SNR dependent
		self.f_min = self.f0 - 3./self.tau
		self.f_max = self.f0 + 3./self.tau
		
		# catch to make sure negative frequencies aren't asked for
		if (self.f_min < 0.0): 
			self.f_min = 0.0
		# ensure max frequency is not greater than the Nyquist frequency
		if (self.f_max > Orbit.f_ny):
			self.f_max = Orbit.f_ny
						
		
	# methods
	get_Psi                = evaluate_wavelet
	calc_Psi               = calc_Psi
	make_padded_Psi        = make_padded_Psi
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
		
def create_instrument_glitch_TDI(instrument_glitch_type, SC_on, Wave, Orbit, SC_point=None):
    """
    Generate the TDI channels for an instrumental glitch.
    
    input:
    
    output:
    
    """
    t = np.arange(0.0, Orbit.Tobs, Orbit.dt) # Todo: Don't need Orbit, its in Wavelet
    N = len(t)
    
    # create empty phase first
    # Todo: be smarter here to save time, i.e. don't twice make phases
    p12 = Phase(1,2, t, np.zeros(N))
    p21 = Phase(2,1, t, np.zeros(N))
    p13 = Phase(1,3, t, np.zeros(N))
    p31 = Phase(3,1, t, np.zeros(N))
    p23 = Phase(2,3, t, np.zeros(N))
    p32 = Phase(3,2, t, np.zeros(N))
    
    # Handle a Laser Phase glitch
    if (instrument_glitch_type == 'Laser Phase'):
        if (SC_point != None):
            raise ValueError("Lase noise is from the laser on one S/C")
            
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = Wavelet(Wave.A, Wave.f0, Wave.tau, \
                            Wave.t0 - Orbit.L/l.Clight, Wave.phi0, Orbit)
        wave_temp.calc_Psi()
        wave_temp.make_padded_Psi(t) 
                
        if (SC_on == 1):
            p12 = Phase(1,2, t, +wave_temp.Psi_padded)
            p13 = Phase(1,3, t, +wave_temp.Psi_padded)
            
            p21 = Phase(2,1, t, -Wave.Psi_padded)
            p31 = Phase(3,1, t, -Wave.Psi_padded)
            
        elif (SC_on == 2):
            p21 = Phase(2,1, t, +wave_temp.Psi_padded)
            p23 = Phase(2,3, t, +wave_temp.Psi_padded)
            
            p12 = Phase(1,2, t, -Wave.Psi_padded)
            p32 = Phase(3,2, t, -Wave.Psi_padded)

        elif (SC_on == 3):
            p31 = Phase(3,1, t, +wave_temp.Psi_padded)
            p32 = Phase(3,2, t, +wave_temp.Psi_padded)
            
            p13 = Phase(1,3, t, -Wave.Psi_padded)
            p23 = Phase(2,3, t, -Wave.Psi_padded)
            
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
            p12 = Phase(1,2, t, Wave.Psi_padded)
            
        elif (SC_on == 2 and SC_point == 1):
            p21 = Phase(2,1, t, Wave.Psi_padded)
            
        elif (SC_on == 1 and SC_point == 3):
            p13 = Phase(1,3, t, Wave.Psi_padded)
            
        elif (SC_on == 3 and SC_point == 1):
            p31 = Phase(3,1, t, Wave.Psi_padded)
            
        elif (SC_on == 2 and SC_point == 3):
            p23 = Phase(2,3, t, Wave.Psi_padded)
            
        elif (SC_on == 3 and SC_point == 2):
            p32 = Phase(3,2, t, Wave.Psi_padded)
        else:
            raise ValueError("Invalid SC_on and/or Sc_point!!!")
    
    # Handle an acceleration noise glitch
    elif (instrument_glitch_type == 'Acceleration'):
        # construct a wavelet whose central time is shifted to t0-L
        wave_temp = Wavelet(Wave.A, Wave.f0, Wave.tau, \
                            Wave.t0 - Orbit.L/l.Clight, Wave.phi0, Orbit)
        wave_temp.calc_Psi()
        wave_temp.make_padded_Psi(t) 
        
        if (SC_on == 1 and SC_point == 2):
            p12 = Phase(1,2, t, -Wave.Psi_padded)
            p21 = Phase(2,1, t, +wave_temp.Psi_padded)
            
        elif (SC_on == 1 and SC_point == 3):
            p13 = Phase(1,3, t, -Wave.Psi_padded)
            p31 = Phase(3,1, t, +wave_temp.Psi_padded)
            
        elif (SC_on == 2 and SC_point == 1):
            p21 = Phase(2,1, t, -Wave.Psi_padded)
            p12 = Phase(1,2, t, +wave_temp.Psi_padded)
            
        elif (SC_on == 3 and SC_point == 1):
            p31 = Phase(3,1, t, -Wave.Psi_padded)
            p13 = Phase(1,3, t, +wave_temp.Psi_padded)
            
        elif (SC_on == 2 and SC_point == 3):
            p23 = Phase(2,3, t, -Wave.Psi_padded)
            p32 = Phase(3,2, t, +wave_temp.Psi_padded)
            
        elif (SC_on == 3 and SC_point == 2):
            p32 = Phase(3,2, t, -Wave.Psi_padded)
            p23 = Phase(2,3, t, +wave_temp.Psi_padded)
            
        else:
            raise ValueError("Invalid SC_on and/or Sc_point!!!")
    
    else:
        raise ValueError("Unexpected instrument glitch type. Choose from 'Laser Phase', 'Optical Path', or 'Acceleration")
    
    # Fourier transform the time-domain phases
    p12.FT_phase(Orbit)
    p21.FT_phase(Orbit)
    p13.FT_phase(Orbit)
    p31.FT_phase(Orbit)
    p23.FT_phase(Orbit)
    p32.FT_phase(Orbit)
    
    tdi = TDI(p12, p21, p13, p31, p23, p32, Orbit)
    
    return tdi


		
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
		
		# if light travel time from SSB to LISA is longer than tau
		# 		add an extra buffer
		#		Todo: Make this a bit smarter
		if (self.hp_wavelet.tau < 3.*8.0*60. or self.hc_wavelet.tau < 3.*8.0*60.):
			self.t_min = self.t_min - 3.*8.0*60.
			self.t_max = self.t_max + 3.*8.0*60.
	
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
	
	