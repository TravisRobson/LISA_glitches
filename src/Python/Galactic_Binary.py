import numpy as np
import LISA as l


IDX_fT     = 0
IDX_cth    = 1
IDX_phi    = 2
IDX_lnA    = 3
IDX_ci     = 4
IDX_psi    = 5
IDX_phi0   = 6
IDX_fdTsq  = 7
IDX_fddTcb = 8

R_on_C  = 497.

pi = np.pi

def construct_basis_tensors(self):
	""" construct the unrotated GW polarizaiton tensors """
	
	u = np.array([self.cth*self.cphi, self.cth*self.sphi, -self.sth])
	v = np.array([self.sphi, -self.cphi, 0.0])

	self.ep = np.outer(u,u) - np.outer(v,v)
	self.ec = np.outer(u,v) + np.outer(v,u)
		
	return
		
def Calc_BW(self):
	""" Estimate an appropriate bandwdith to simulate signal over """
	
	N_LISA = 0.0 # spreading due to LISA modulation
	N_ss   = 0.0 # sinc (finite time) spreading
	N_RR   = 0.0 # radiation reaction spreading

	SnM = self.Orbit.get_Sn(self.f)/(4.*np.sin(self.f/l.fstar)**2)
	SNR_estimate = self.A*np.sqrt(self.T/SnM[0])

	N_RR += self.fdTsq + 0.5*self.fddTcb

	# result of a fit, i.e. J_{k}(x)^{2} is the power in each harmonic
	#	 how many harmonics to maintain power to a part in 100
	x = 2.*np.pi*self.sth*R_on_C*self.f
	N_LISA = 0.5*(1.01788 + 1.02936*x);

	N_ss = 4.0

	# SNR**2 adjustment based on 1-FF = (D-1)/(2*SNR**2)
	N = N_RR + (SNR_estimate/10.)**2*(N_ss + N_LISA)

	self.N = int( 2**( np.ceil(np.log2(N)))/2)
	
	return
	
def calc_xi(self):
	""" Calculate the wave variable for each S/C """
	
	k_dot_x = self.k[0]*self.x[0,:,:] + self.k[1]*self.x[1,:,:] + self.k[2]*self.x[2,:,:]

	self.xi = self.t.reshape((1, self.N)) - k_dot_x/l.Clight
	
	return
	
def calc_f(self):
	""" Calculate the Doppler-shifted GW frequuency at each S/C """

	self.f_at_SC = self.f + self.fdot*self.xi + 0.5*self.fddot*self.xi**2	
	
	return
	
def calc_k_dot_r(self):
	""" Calculate the dot product between the propagation direction and S/C separation vectors """
	
	self.k_dot_r = self.k[0]*self.rij[0,:,:,:] + self.k[1]*self.rij[1,:,:,:] + self.k[2]*self.rij[2,:,:,:]
	
	return
	
def calc_contracted_detector_tensors(self):
	""" Calculate the S/C unit-separation vector outer product between itself 
		and then contract with the detector tensor """
		
	r_outer_r = self.Orbit.get_rOr(self.rij)
	
	# perform the contraction over the spatial dimensions of the arrays
	self.dp = np.einsum('mnijk,mn', r_outer_r, self.ep)
	self.dc = np.einsum('mnijk,mn', r_outer_r, self.ec) # result dim: (S/C i, S/C j, len(t))
	
	return
	
def calc_trans(self):
	""" Calculate the trasnfer function """
	
	Ap = self.A*(1. + self.ci**2)
	Ac = -2.0*self.A*self.ci
	
	# rotate the tensor by the GW polarization angle
	DPr = Ap*self.c2psi
	DPi = -Ac*self.s2psi
	DCr = -Ap*self.s2psi
	DCi = -Ac*self.c2psi
	
	self.q  = int(self.f*self.T) # the carrier bin (floored)
	df = 2.*pi*self.q/self.T # the distance (in freq) from carrier freq to freq bin

	temp_arg1 = 1. - self.k_dot_r

	arg1 = np.zeros((3,3,self.N))
	temp = 1./l.fstar
	arg1[0,1] = 0.5*self.f_at_SC[0]*temp*temp_arg1[0,1,:]
	arg1[0,2] = 0.5*self.f_at_SC[0]*temp*temp_arg1[0,2,:]
	arg1[1,0] = 0.5*self.f_at_SC[1]*temp*temp_arg1[1,0,:]
	arg1[1,2] = 0.5*self.f_at_SC[1]*temp*temp_arg1[1,2,:]
	arg1[2,0] = 0.5*self.f_at_SC[2]*temp*temp_arg1[2,0,:]
	arg1[2,1] = 0.5*self.f_at_SC[2]*temp*temp_arg1[2,1,:]

	arg2  = 2.*pi*self.f*self.xi + self.phi0 - df*self.t \
	        + pi*self.fdot*self.xi**2 + 1./3.*pi*self.fddot*self.xi**3

	sinc = 0.25*np.sinc(arg1/pi)

	aevol = 1.0 + 2./3.*self.fdot/self.f*self.xi

	tran1r = aevol*(self.dp*DPr + self.dc*DCr) 
	tran1i = aevol*(self.dp*DPi + self.dc*DCi)
	
	arg2 = arg2.reshape((3,1,self.N))
	tran2r = np.cos(arg1 + arg2)
	tran2i = np.sin(arg1 + arg2)
	
	self.TR = sinc*(tran1r*tran2r - tran1i*tran2i)
	self.TI = sinc*(tran1r*tran2i + tran1i*tran2r)
	
	return
	
def calc_TDI(self):
	""" Calculate the TDI channels through phase shifts in frequency domain """
	
	fonfs = self.freqs/l.fstar
	
	phase1 = np.cos(fonfs)    - 1.0j*np.sin(fonfs)
	phase2 = np.cos(2.*fonfs) - 1.0j*np.sin(2.*fonfs)
	phase3 = np.cos(3.*fonfs) - 1.0j*np.sin(3.*fonfs)
	
	self.X =   (self.data_FT[0,1,:] - self.data_FT[0,2,:])*phase3 \
			 + (self.data_FT[1,0,:] - self.data_FT[2,0,:])*phase2 \
			 + (self.data_FT[0,2,:] - self.data_FT[0,1,:])*phase1 \
			 + (self.data_FT[2,0,:] - self.data_FT[1,0,:])
			 
	self.Y =  (self.data_FT[1,2,:] - self.data_FT[1,0,:])*phase3 \
		 	+ (self.data_FT[2,1,:] - self.data_FT[0,1,:])*phase2 \
			+ (self.data_FT[1,0,:] - self.data_FT[1,2,:])*phase1 \
		 	+ (self.data_FT[0,1,:] - self.data_FT[2,1,:])
	
	self.Z =  (self.data_FT[2,0,:] - self.data_FT[2,1,:])*phase3 \
		 	+ (self.data_FT[0,2,:] - self.data_FT[1,2,:])*phase2 \
		 	+ (self.data_FT[2,1,:] - self.data_FT[2,0,:])*phase1 \
		 	+ (self.data_FT[1,2,:] - self.data_FT[0,2,:])
		 	
#	The phase differences between LISA simulator and Neil/Tyson's code
# 	phiLS = 2.*np.pi*self.f*(self.Orbit.dt/2.0-self.Orbit.L/l.Clight);
# 	t1 = np.copy(self.X.real)
# 	t2 = np.copy(self.X.imag)
# 	self.X.real = np.cos(phiLS)*t1 - np.sin(phiLS)*t2
# 	self.X.imag = -(np.sin(phiLS)*t1 + np.cos(phiLS)*t2)
	
	self.A = 1./3.*(2.*self.X - self.Y - self.Z)
	self.E = 1./np.sqrt(3.)*(self.Z - self.Y)
	self.T = 1./3.*(self.X + self.Y + self.Z)
	
	return
	

class GB:
	def __init__(self, params, Orbit):
		self.params = params
		self.Orbit  = Orbit
		self.T      = self.Orbit.Tobs # observation period
		
		# extract parameters, and calculate some respective pieces
		#	to be used in waveform generation
		self.fT         = params[IDX_fT]
		self.f			= self.fT/self.T
		
		self.cth        = params[IDX_cth]
		self.sth        = np.sqrt(1. - self.cth**2)
		
		self.phi        = params[IDX_phi]
		self.cphi       = np.cos(self.phi)
		self.sphi	    = np.sin(self.phi)
		
		self.lnA        = params[IDX_lnA]
		self.A          = np.exp(self.lnA)
		
		self.ci         = params[IDX_ci]
		
		self.psi		= params[IDX_psi]
		self.c2psi      = np.cos(2.*self.psi)
		self.s2psi      = np.sin(2.*self.psi)
		
		self.phi0       = params[IDX_phi0]
		
		self.fdTsq  = params[IDX_fdTsq]
		self.fdot   = self.fdTsq/(self.T**2)
		
		self.fddTcb = params[IDX_fddTcb]
		self.fddot   = self.fddTcb/(self.T**3)
		
		# unit vector pointing towards GB from SSB (Solar System Barycenter)
		self.k = -np.array([self.sth*self.cphi, self.sth*self.sphi, self.cth])
		
		construct_basis_tensors(self)	

		Calc_BW(self)
		#print 'N........ {}'.format(self.N)
		
		# samples in time for slowly evolving piece
		self.t = np.arange(0., Orbit.Tobs, Orbit.Tobs/self.N)
		
		# calculate all spacecraft position
		# dim, S/C, time
		#self.x = np.zeros((3, 3, self.N))
		self.x = self.Orbit.get_orbit(self.t)
		
		calc_xi(self) # calculate the wave variable for all S/C positions
		
		calc_f(self) # calculate frequencies as seen be S/C
		
		# calculate separation vectors between spacecraft
		self.rij = self.Orbit.get_seps(self.t, self.x)

		calc_k_dot_r(self) # calculate k dot rij
		
		# calculate detector matrices i.e. detector tensors contracted with polarization tensors
		calc_contracted_detector_tensors(self)
		
		calc_trans(self) # Calculate the transfer function
		
		self.data = self.TR + 1.0j*self.TI # Calculate the slow time series

		self.data_FT = np.fft.fft(self.data) # FFT the slow time series

		# construct real results
		self.data_FT = np.fft.fftshift(self.data_FT, axes=(2,))*0.5/self.N*np.sqrt(Orbit.Tobs)
		f_lo = (self.q - self.N/2)/self.T
		f_hi = (self.q + self.N/2)/self.T
		self.freqs = np.arange(f_lo, f_hi, 1./self.T)
		
		calc_TDI(self) # Form up the TDI variables





