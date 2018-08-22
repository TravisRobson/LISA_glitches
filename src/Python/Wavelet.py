import numpy as np

IDX_lnA  = 0
IDX_f0   = 1 # in mHz
IDX_t0   = 2 # in weeks
IDX_tau  = 3 # in weeks
IDX_phi0 = 4 # already dimensionless


def evaluate_wavelet(self, t):
    """ Return the value of the Sine-Gaussian wavelet at time(s) t """	

    arg1 = 2.*np.pi*self.f0*(t - self.t0) + self.phi0
    arg2 = (t - self.t0)/self.tau

    result = self.A*(np.cos(arg1) + 1.0j*np.sin(arg1))*np.exp(-arg2**2)

    return result
	
def calc_Psi(self):
	""" In the SSB frame, calculate relevant times and sample the wavelet """

	self.t = np.arange(self.t_min, self.t_max, self.Orbit.dt)
	self.Psi = self.get_Psi(self.t)
	
	return
	
def make_padded_Psi(self, t):
	""" Pad the wavelet's internal Psi to match shape of the provided time series """
	
	# find where input time array bounds matches the wavelet's native time array
	
	left_pad  = np.argwhere(t == self.t[0] ).flatten()[0]
	right_pad = len(t) - 1 - np.argwhere(t == self.t[-1]).flatten()[0]

	self.Psi_padded = np.pad(self.Psi, (left_pad,right_pad), 'constant')
	
	return


class Wavelet:
	""" Sine-Gaussian wavelet class """
	
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
			self.t_max = self.t_min + 2*Orbit.dt
		
		# Todo: how to make sure time doesn't exceed LISA observation times
		if (self.t_min < 0.0): # ensure time is positive
			self.t_min = 0.0
		
		# Todo: smarter choices for these bounds, probably SNR dependent
		self.f_min = self.f0 - 2.5/self.tau
		self.f_max = self.f0 + 2.5/self.tau
		
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
	
	
	
	
	
	
	
	