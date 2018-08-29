import numpy as np


fm     = 3.168753575e-8   # LISA modulation frequency
YEAR   = 31457280.0 	  # LDC standard for year in seconds
AU     = 1.49597870660e11 # Astronomical unit (meters)
Clight = 299792458.       # speed of light (m/s)
TSUN   = 4.9169e-6        # mass of sun (seconds)
Sps    = 8.321000e-23     # Photon shot power noise
Sacc   = 9.000000e-30     # Acceleration power noise
fstar  = 0.01908538063694777 # transfer frequency

pi = np.pi

def full_SC_orbit(self, t):
	""" Calculate the analytic (leading order in eccentricity) LISA orbits """
	
	N = len(t)
	ecc = self.ecc
	
	alpha = (2.*pi*fm*t + self.kappa).reshape((1,N))
	sa = np.sin(alpha) # to later go to (S/C, len(t))
	ca = np.cos(alpha)
	
	beta = (np.array([0.0, 2.*pi/3., 4.*pi/3.]) + self.Lambda).reshape((3,1))
	sb = np.sin(beta)
	cb = np.cos(beta) # (S/C, len(t))
	
	# dim, S/C, time
	x = np.zeros((3, 3, N))
	
	x[0] = AU*ca + AU*ecc*(sa*ca*sb - (1. + sa*sa)*cb)
	x[1] = AU*sa + AU*ecc*(sa*ca*cb - (1. + ca*ca)*sb)
	x[2] = -np.sqrt(3.)*AU*ecc*(ca*cb + sa*sb)


	return x
	
def full_SC_seps(self, t, x):
	""" calculate S/C separation vectors """
	
	N = len(t)
		
	rij = np.zeros((3,3,3,N))

	rij[:,0,1,:] = x[:,1,:] - x[:,0,:]
	rij[:,1,0,:] = -rij[:,0,1,:]

	rij[:,0,2,:] = x[:,2,:] - x[:,0,:]
	rij[:,2,0,:] = -rij[:,0,2,:]

	rij[:,1,2,:] = x[:,2,:] - x[:,1,:]
	rij[:,2,1,:] = -rij[:,1,2,:]
	
	return rij/self.L
	
def get_instrument_noise(self, f, X_flag=None):
    """
    Get the instrument noise power spectral density
    """

    f_star = Clight/(2.*np.pi*self.L) # TODO: will change when we go to second gen TDI
    fonfs  = f/f_star
    
#     huh = np.min(f)
#     if (huh<1.0e-4):
#         print('Go fuck yourself.... {}', huh)
    
    red = 16.0*( (2.0e-5/f)**10.0 + (1.0e-4/f)**2. )
    red[np.isnan(red)] = 1.
        
    Sloc = 2.89e-24;

    trans = np.sin(fonfs)**2.0 # transfer function

    if (X_flag==None):

        SnAE = 16./3.*trans*( (2.0+np.cos(fonfs))*(Sps + Sloc)  + \
                        2.0*(3.0 + 2.0*np.cos(fonfs) + np.cos(2.0*fonfs)) \
                           *(Sloc/2.0 + Sacc/(2.0*np.pi*f)**4.0*(1.0+red)))/(2.0*self.L)**2.0

        SnT = 16./3.*trans*(0.5*(Sps + Sloc)*(1. - np.cos(fonfs)) + \
                           (1. - 2.*np.cos(fonfs) + np.cos(fonfs)**2) \
                          *(Sloc/2.0 + Sacc/(2.0*np.pi*f)**4.0*(1.0+red)))/(2.0*self.L)**2.0 
            
        

        return SnAE, SnT

    else:
        SXYZ = 4.0*trans*(4.0*(Sps+Sloc) + 8.0*(1.0+np.cos(fonfs)**2.0)*(Sloc/2.0 + Sacc/(2.0*np.pi*f)**4.*(1.0+red)))/(2.0*self.L)**2.0

    return SXYZ

def get_Sn(f):
    """
    Get the instrument noise power spectral density
    """

    f_star = Clight/(2.*np.pi*2.5*10**9) # TODO: will change when we go to second gen TDI
    fonfs  = f/fstar

    red  = 16.0*( (2.0e-5/f)**10.0 + (1.0e-4/f)**2. )
    Sloc = 2.89e-24;

    SXYZ = 4.0*(4.0*(Sps+Sloc) + 8.0*(1.0+np.cos(fonfs)**2.0)*(Sloc/2.0 + Sacc/(2.0*np.pi*f)**4.*(1.0+red)))/(2.0*2.5*10**9)**2.0

    return SXYZ


def construct_r_outer_r(self, rij):
	""" Take the outer product of \hat{r}_{ij} shaped as Dim, S/C i, S/C j, len(t) """

	# find number of time samples
	N = np.shape(rij)[-1]
	
	r_outer_r = rij.reshape((3,1, 3,3, N))*rij.reshape((1,3, 3,3, N))

	return r_outer_r

class Orbit:

	def __init__(self, Tobs, type='analytic', kappa=0.0, Lambda=0.0, dt=15.0, ecc=0.0048241852, L=2.5e9):
		self.Tobs = Tobs     # observation period
		self.type = type     # 'analytic' or 'numeric'
		
		# check to make sure that a supported orbit type has been selections
		if (self.type != 'analytic'):
			raise ValueError("Orbit type incorrect or not supported! We currently only support 'analytic'.")
			
		self.kappa = kappa          # initial azimuthal position of guiding center
		self.Lambda = Lambda        # initial orientation of LISA constellation
		self.dt = dt                # cadence (sec)
		self.f_ny = 1./(2.*self.dt) # Nyquist frequency (Hz)
		
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
	get_Sn = get_instrument_noise
	
	get_orbit = full_SC_orbit
	get_seps = full_SC_seps
	get_rOr = construct_r_outer_r
	
	
	

	
	
	
	
	
	
	
