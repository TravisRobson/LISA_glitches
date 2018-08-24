import numpy as np

import LISA as l

############## LISA Phase Measurement Class Implementation #################

def FT_phase(self, Orbit):
	""" Take the Fourier transform of the Phase object """
	
	phi = np.fft.fft(self.phi)
	
	self.phi_FT = phi[:len(self.t)//2]*Orbit.Tobs/len(self.t)*0.5 # just need the positive frequencies

	self.freqs  = np.fft.fftfreq(self.t.shape[-1], Orbit.dt)[:len(self.t)//2] 

	return
	
def FT_phase_FAST(self, Orbit):
	""" Take the Fourier transform of the Phase object """
	
	phi = np.fft.fft(self.phi)
	
	self.phi_FT = phi*Orbit.Tobs/len(self.t)*0.5

	return
	

class Phase:
	""" LISA laser phase comparison """
	
	def __init__(self, i, j, t, phi):
		self.i      = i 	 # emitting S/C
		self.j      = j	     # receiving S/C
		self.t      = t      # time array
		self.phi    = phi    # phase array
		
	FT_phase = FT_phase		
	FT_phase_FAST = FT_phase_FAST
	
	
################### LISA TDI Class Implementation ###################
	
def get_TDI_snr(self, f_min, f_max, X_flag=None):
    """ Calculate the SNR in A,E, and T channels """
    
    if (len(self.freqs) < 2): # i.e. the signals do not overlap at all..
        return 0.0

    mask = (self.freqs>f_min) & (self.freqs<f_max)
    df = self.freqs[1] - self.freqs[0]
    
    if (X_flag == None):
        SnAE, SnT = self.Orbit.get_Sn(self.freqs[mask])

        SNR_A = 4.0*np.sum( self.A[mask]*np.conjugate(self.A[mask])/SnAE).real*df
        SNR_E = 4.0*np.sum( self.E[mask]*np.conjugate(self.E[mask])/SnAE).real*df

        SNR_T = 4.0*np.sum( self.T[mask]*np.conjugate(self.T[mask])/SnT).real*df

        return SNR_A, SNR_E, SNR_T	
    else:
        SnX = self.Orbit.get_Sn(self.freqs[mask], X_flag)

        SNR_X = 4.0*np.sum( self.X[mask]*np.conjugate(self.X[mask])/SnX).real*df

        return SNR_X
        	
def get_TDI_overlap(tdi1, tdi2, f_min, f_max, X_flag=None):
    """ Calculate the overlap between TDI data """

    mask1 = (tdi1.freqs >= f_min) & (tdi1.freqs <= f_max)
    mask2 = (tdi2.freqs >= f_min) & (tdi2.freqs <= f_max) 

    freqs = tdi1.freqs[mask1] 
    
    if (len(freqs) < 2): # i.e. the signals do not overlap at all..
        return 0.0

    df = freqs[1] - freqs[0]

    if (X_flag == None):
        d1_A = tdi1.A[mask1]
        d1_E = tdi1.E[mask1]
        d1_T = tdi1.T[mask1]

        d2_A = tdi2.A[mask2]
        d2_E = tdi2.E[mask2]
        d2_T = tdi2.T[mask2]
        
        SnAE, SnT = tdi1.Orbit.get_Sn(freqs)

        overlap_A = 0
        overlap_E = 0
        overlap_T = 0

        overlap_A += np.sum(d1_A*np.conjugate(d2_A)/SnAE)
        overlap_E += np.sum(d1_E*np.conjugate(d2_E)/SnAE)
        overlap_T += np.sum(d1_T*np.conjugate(d2_T)/SnT )

        return 4*np.array([overlap_A, overlap_E, overlap_T]).real*df

    else:
        d1_X = tdi1.X[mask1]
        d2_X = tdi2.X[mask2]

        SnX  = tdi1.Orbit.get_Sn(freqs, X_flag)

        overlap_X  = 0
        overlap_X += np.sum(d1_X*np.conjugate(d2_X)/SnX)

        return 4*overlap_X.real*df

	
class TDI:
	""" Time Delay Inferometery class X, Y, Z, A, E, and T data channels """
	
	def __init__(self, phi12, phi21, phi13, phi31, phi23, phi32, Orbit):
		""" Take the phase comparisons and construct the TDI data channels """

		self.freqs = phi12.freqs # TODO: need some logic to handle if another phase is being used instead
		self.Orbit = Orbit
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
		 	
		 			
	# Methods
	get_TDI_snr = get_TDI_snr
	
	
	
	
	
	
	
	
	
	