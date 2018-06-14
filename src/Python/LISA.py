import numpy as np


fm     = 3.168753575e-8   # LISA modulation frequency
YEAR   = 31457280.0 	  # LDC standard for year in seconds
AU     = 1.49597870660e11 # Astronomical unit (meters)
Clight = 299792458.       # speed of light (m/s)
TSUN   = 4.9169e-6        # mass of sun (seconds)
Sps    = 8.321000e-23     # Photon shot power noise
Sacc   = 9.000000e-30     # Acceleration power noise
fstar  = 0.01908538063694777 # transfer frequency


def calc_SC_x(self, i, t):
	"""
	Get the ith (1,2,3) S/C's position at time t
	"""
	
	if (i!=1 and i!=2 and i!=3):
		raise ValueError("There only only 3 S/C! Choose 1, 2, or 3.")
	
	alpha = 2.*np.pi*fm*t + self.kappa 
	
	if (i==1):
		beta = 0.          
	if (i==2):
		beta = 2.*np.pi/3. 
	if (i==3):
		beta = 4.*np.pi/3. 
	
	beta = beta + self.Lambda

	sa = np.sin(alpha)
	ca = np.cos(alpha)
	
	sb = np.sin(beta)
	cb = np.cos(beta)
	
	self.x[i-1][0] = AU*ca + AU*self.ecc*(sa*ca*sb - (1. + sa*sa)*cb);
	self.x[i-1][1] = AU*sa + AU*self.ecc*(sa*ca*cb - (1. + ca*ca)*sb);
	self.x[i-1][2] = -np.sqrt(3.)*AU*self.ecc*(ca*cb + sa*sb);
	
	return
	

def full_SC_orbit(self, t, x):
	N = len(t)
	
	alpha = 2.*np.pi*fm*t + self.kappa 
	sa = np.sin(alpha).reshape((1,N))
	ca = np.cos(alpha).reshape((1,N))
	
	beta = np.array([0.0, 2.*np.pi/3., 4.*np.pi/3.]) + self.Lambda
	sb = np.sin(beta).reshape((3,1))
	cb = np.cos(beta).reshape((3,1))
	
	# dim, S/C, time
	# x = np.zeros((3, 3, N))

	x[0] = x[0] + AU*ca + AU*self.ecc*(sa*ca*sb - (1. + sa*sa)*cb)
	x[1] = x[1] + AU*sa + AU*self.ecc*(sa*ca*cb - (1. + ca*ca)*sb)
	x[2] = x[2] - np.sqrt(3.)*AU*self.ecc*(ca*cb + sa*sb)
	
	return x
	
def get_SC_x(self, i):
	"""
	Return the desired vector
	"""
		
	return self.x[i-1]
	
def calc_SC_sep(self, i, j, t):
	"""
	Calculate the separation vectors for the specified SC
	"""
	
	# which separation vector
	if (i==1 and j==2):
		k = 0
	elif (i==2 and j==1):
		k = 1
	elif (i==1 and j==3):
		k = 2
	elif (i==3 and j==1):
		k = 3
	elif (i==2 and j==3):
		k = 4
	elif (i==3 and j==2):
		k = 5
			
	# calculate the S/C positions
	self.calc_x(i, t)
	self.calc_x(j, t)
		
	# take the difference
	self.rij[k] = self.get_x(j) - self.get_x(i)

	return
	
def full_SC_seps(self, t, x, rij):
 	N = len(t)

	rij = (x.reshape((3,1,3,N)) - x.reshape((3,3,1,N)))/self.L
	
	return rij
	
def get_SC_sep(self, i, j):
	"""
	Return the desired separation vector
	"""
	
	if (i==1 and j==2):
		k = 0
	elif (i==2 and j==1):
		k = 1
	elif (i==1 and j==3):
		k = 2
	elif (i==3 and j==1):
		k = 3
	elif (i==2 and j==3):
		k = 4
	elif (i==3 and j==2):
		k = 5
		
	return self.rij[k]
	
def get_instrument_noise(self, f):
	"""
	Get the instrument noise power spectral density
	"""
	
	f_star = Clight/(2.*np.pi*self.L) # TODO: will change when we go to second gen TDI
	fonfs  = f/f_star

	red  = 16.0*( (2.0e-5/f)**10.0 + (1.0e-4/f)**2. )
	Sloc = 2.89e-24;

	trans = pow(np.sin(fonfs), 2.0); # transfer function

	#     SnAE = 16.0/3.0*trans*( (2.0+np.cos(fonfs))*(Sps + Sloc) 
	#     					    +2.0*( 3.0 + 2.0*cos(fonfs) + cos(2.0*fonfs) )
	#     					        *( Sloc/2.0 + Sacc/(2.0*np.pi*f)**4.0*(1.0+red) ) )
	#     					  /(2.0*self.L)**2.0;

	SXYZ = 4.0*trans*(4.0*(Sps+Sloc) + 8.0*(1.0+np.cos(fonfs)**2.0)*(Sloc/2.0 + Sacc/(2.0*np.pi*f)**4.*(1.0+red)))/(2.0*self.L)**2.0

	return SXYZ
	
class Orbit:
	# 

	def __init__(self, Tobs, type='analytic', kappa=0.0, Lambda=0.0, dt=15.0, ecc=0.0048241852, L=2.5e9):
		self.Tobs = Tobs     # observation period
		self.type = type     # 'analytic' or 'numeric'
		
		# check to make sure that a supported orbit type has been selections
		if (self.type != 'analytic'):
			raise ValueError("Orbit type incorrect or not supported! We currently only support 'analytic'.")
			
		self.kappa = kappa   # initial azimuthal position of guiding center
		self.Lambda = Lambda # initial orientation of LISA constellation
		self.dt = dt         # cadence (sec)
		
		# throw an if number of samples is not a power of 2
		num = np.int(self.Tobs/self.dt)
		if (num != 0 and ((num & (num - 1)) == 0)):
			pass
		else:
			raise ValueError("Number of data samples is not a power of 2!")
		
		self.ecc = ecc       # eccentricity of S/C orbits
		self.L = L           # (meters) average S/C separation
		
		# Cache quantities for easy use of object
		self.x   = np.zeros((3,3)) # position vectors, 3 S/C's, and 3 spatial dimensions
		self.rij = np.zeros((6,3)) # separation vectors (1,2), (2,1), (1,3), (3,1), (2,3), and (3,2)
		self.Lij = np.zeros(6)     # separation distances
		
	# Methods
	calc_x = calc_SC_x
	get_x  = get_SC_x
	
	calc_r = calc_SC_sep
	get_r  = get_SC_sep
	get_Sn = get_instrument_noise
	
	fill_full_orbit = full_SC_orbit
	fill_full_seps = full_SC_seps
	
	
	

	
	
	
	
	
	
	
