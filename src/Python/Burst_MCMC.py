import numpy as np
import scipy.stats as ss
import copy

import Wavelet as wv
import Glitch as gl
import Burst as bu

import TDI as td

def f0_sample(self):
    """ Generate a sample from the uniform/Gaussian hybrid """
    
    x = ss.uniform.rvs(0,1,1)
    
    if (x<self.uniform_frac):
        samp = self.uni_dist.rvs()
    
    else:
        samp = self.norm_dist.rvs()
    
    return samp
    
def f0_pdf(self, f0):
	""" find the value of the probability distribution function at f0 """
	
	prob  = self.uniform_frac*self.uni_dist.pdf(f0)
	prob += self.gaussian_frac*self.norm_dist.pdf(f0)
	
	return np.log(prob)


class f0_proposal:
    """ 'Cheat' propsal distribution for wavelet carrier frequency """
    
    def __init__(self, f0, f_min, f_max, sigma):
            self.f0    = f0
            self.f_min = f_min
            self.f_max = f_max
            self.sigma = sigma
            
            self.uni_dist  = ss.uniform(self.f_min, self.f_max-self.f_min)
            self.norm_dist = ss.norm(self.f0, self.sigma)
            
            self.uniform_frac  = 0.2
            self.gaussian_frac = 0.8
            
    sample = f0_sample
    pdf    = f0_pdf

    
def t0_sample(self):
    """ Generate a sample from the uniform/Gaussian hybrid """
    
    x = ss.uniform.rvs(0,1,1)
    
    if (x<self.uniform_frac):
        samp = ss.uniform.rvs(self.t_min, self.t_max - self.t_min,1)[0]
    
    else:
        samp = ss.norm.rvs(self.t0, self.sigma,1)[0]

    return samp
    
def t0_pdf(self, t0):
	""" find the value of the probability distribution fucntion at f0 """
	
	prob  = self.uniform_frac*ss.uniform.pdf(t0, self.t_min, self.t_max-self.t_min)
	prob += self.gaussian_frac*ss.norm.pdf(t0, self.t0, self.sigma)
	
	return np.log(prob)


class t0_proposal:
    """ 'Cheat' propsal distribution for wavelet carrier frequency """
    
    def __init__(self, t0, t_min, t_max, sigma):
            self.t0    = t0
            self.t_min = t_min
            self.t_max = t_max
            self.sigma = sigma
            
            self.uniform_frac  = 0.2
            self.gaussian_frac = 0.8
            
    sample = t0_sample
    pdf    = t0_pdf

def GW_sample(self, x):
	""" generate a sample for the GW proposal distribution """
	
	# Extract current parameters that aren't built into t0, f0 proposal distributions
	lnA_x  = x[0]
	phi0_x = x[1]
	tau_x  = x[2]
	ep_x   = x[3]
	
	# Todo: will need to screw around with variance choices I'm sure
	lnA_y  = ss.norm.rvs(lnA_x,  self.sgm_lnA)
	phi0_y = ss.norm.rvs(phi0_x, self.sgm_phi0)
	tau_y  = ss.norm.rvs(tau_x,  self.sgm_tau)
	ep_y   = ss.norm.rvs(ep_x,   self.sgm_ep)
	
	t0_y = self.t0_prop.sample()
	f0_y = self.f0_prop.sample()
	
	return np.array([lnA_y, phi0_y, tau_y, ep_y, t0_y, f0_y])
	
def GW_pdf(self, x, y):
	""" Caclulate the pdf value of the GW proposal distribution """
	
	# Extract current parameters that aren't built into t0, f0 proposal distributions
	lnA_x  = x[0]
	phi0_x = x[1]
	tau_x  = x[2]
	ep_x   = x[3]
# 	t0_x   = x[4]
# 	f0_x   = x[5]
	
	lnA_y  = y[0]
	phi0_y = y[1]
	tau_y  = y[2]
	ep_y   = y[3]
	t0_y   = y[4]
	f0_y   = y[5]
	
	# Todo: I don't like how I would have to change the variance in two places
	prob  = ss.norm.logpdf(lnA_y, lnA_x, self.sgm_lnA)
	prob += ss.norm.logpdf(phi0_y, phi0_x, self.sgm_phi0)
	prob += ss.norm.logpdf(tau_y, tau_x, self.sgm_tau)
	prob += ss.norm.logpdf(ep_y, ep_x, self.sgm_ep)
	prob += self.t0_prop.pdf(t0_y)
	prob += self.f0_prop.pdf(f0_y)
	
	return prob
	
    
class GW_proposal:
	""" Proposal distribution for GW (lnA, phi0, tau, ep, t0, f0) """    
	
	def __init__(self, sgm_lnA, sgm_phi0, sgm_tau, sgm_ep, t0_prop, f0_prop,\
						lnA_true, phi0_true, tau_true, ep_true):
						
		# the variances for gaussian proposal distributions
		self.sgm_lnA  = sgm_lnA
		self.sgm_phi0 = sgm_phi0
		self.sgm_tau  = sgm_tau
		self.sgm_ep   = sgm_ep
		
		# the true parameters
		self.lnA_true  = lnA_true
		self.phi0_true = phi0_true
		self.tau_true  = tau_true
		self.ep_true   = ep_true
		
		# these proposals will have to be setup first as its currently implemented
		self.t0_prop = t0_prop
		self.f0_prop = f0_prop
		
	
	# methods
	sample = GW_sample
	pdf    = GW_pdf
	
	
# Will also need some prior class
# and then I should be able to start the MCMC
	
	
def overlap(tdi_a, tdi_b):

	return	
	
def GW_MCMC(N, Orbit, TDI, GW, seed, flag=0, Flag_fast=0):
	""" Perform the MCMC on the GW burst """
	
	t = np.arange(0.0, Orbit.Tobs, Orbit.dt) # set up the time of this orbit

	data_snr = np.sum(td.get_TDI_overlap(TDI, TDI, GW.hp_wavelet.f_min, GW.hp_wavelet.f_max))
# 	data_snr = gl.get_TDI_overlap(TDI, TDI, GW.hp_wavelet.f_min, GW.hp_wavelet.f_max)[0]
	SNR = np.sqrt(data_snr)
	print("Data SNR: {}\n".format(SNR))
	np.random.seed(seed)
	
	# Extract the true parameters
	lnA_true   = np.log(GW.hp_wavelet.A)
	phi0_true  = GW.hp_wavelet.phi0
	tau_true   = GW.hp_wavelet.tau
	epsil_true = GW.hc_wavelet.A/GW.hp_wavelet.A
	t0_true    = GW.hp_wavelet.t0
	f0_true    = GW.hp_wavelet.f0
	
	# setup the proposal distribution
	df = TDI.freqs[1] - TDI.freqs[0]
	
	temp = (np.pi*tau_true*f0_true)**2
	prop_sigma_lnA = 1/SNR #np.sqrt( (1 + np.exp(-2*temp))  )/SNR
	prop_sigma_phi0 = 1/SNR #np.sqrt( (1 - np.exp(-2*temp)) )/SNR
	prop_sigma_tau = 2./np.sqrt(3.)*tau_true/SNR#/np.sqrt(1 + (16/3*temp**2-8*temp)*np.exp(-temp))
	prop_sigma_ep = 1/SNR
	prop_sigma_f0 = 1/(SNR*np.pi*tau_true)#/np.sqrt(1 + 4*temp*np.exp(-2*temp))
	prop_sigma_t0 = 1/(2*np.pi*f0_true*SNR)
	
	f0_prop = f0_proposal(f0_true, GW.hp_wavelet.f_min, GW.hp_wavelet.f_max, prop_sigma_f0) 

	t0_prop = t0_proposal(t0_true, t0_true - 30.0*Orbit.dt/SNR, t0_true + 30.*Orbit.dt/SNR, prop_sigma_t0)

	gw_prop = GW_proposal(prop_sigma_lnA, prop_sigma_phi0, prop_sigma_tau, prop_sigma_ep, t0_prop, f0_prop,\
                          lnA_true, phi0_true, tau_true, epsil_true)	
	
	# setup state and its likelihood
	x = np.array([lnA_true, phi0_true, tau_true, epsil_true, t0_true, f0_true])
	tdi_x = copy.deepcopy(TDI)
	
	# calculate the likelihood
	if (flag==1):
		lkl_x = 1.0
	else:
		lkl_x = data_snr - 2.*np.sum(td.get_TDI_overlap(TDI,   tdi_x, GW.hp_wavelet.f_min, GW.hp_wavelet.f_max)) \
							+ np.sum(td.get_TDI_overlap(tdi_x, tdi_x, GW.hp_wavelet.f_min, GW.hp_wavelet.f_max))
		lkl_x *= -0.5

	# the output chain
	out = np.zeros( (N, len(x)+1) )
	props = np.zeros((N,len(x)))
	accept = 0
	
	tdi_y = copy.deepcopy(TDI) # Todo: prolly doesn't need to happen every time

	for i in range(N):
		# Lets sort out the proposed state
		y = gw_prop.sample(x)
		#y[0] = lnA_true
		#y[1] = phi0_true
		#y[2] = tau_true
		#y[3] = epsil_true
		#y[4] = t0_true
		#y[5] = f0_true
		props[i] = y
		
		# check the prior range
		if not (np.log(1.0e-25) < y[0] < np.log(1.0e-15)):
			out[i,0]  = lkl_x
			out[i,1:] = np.copy(x)
			continue	
			
		elif not (-np.pi < y[1] < np.pi):
			out[i,0]  = lkl_x
			out[i,1:] = np.copy(x)
			continue
			
		elif not (0.1 < y[2] < Orbit.Tobs):
			out[i,0]  = lkl_x
			out[i,1:] = np.copy(x)
			continue
				
		elif not (-1 < y[3] < 1):
			out[i,0]  = lkl_x
			out[i,1:] = np.copy(x)
			continue
		
		w1_y = wv.Wavelet(     np.exp(y[0]), y[5], y[2], y[4], y[1], Orbit)
		w2_y = wv.Wavelet(y[3]*np.exp(y[0]), y[5], y[2], y[4], y[1], Orbit)
						 ### (epsil*A, f0, tau, t0, phi0, orb)
	
		# Todo: I am not marginalizing over sky angles yet...
		phi   = GW.phi
		theta = GW.theta
		psi   = GW.psi
	
		gw_y = bu.GW_glitch(w1_y, w2_y, theta, phi, psi, Orbit, Flag_fast)
		
		# calculate the likelihood
		if (flag==1):
			lkl_y = 1.0
		else:
			gw_y.construct_detector_tensor() # Doesn't need to be calculated if done before
			
			gw_y.calculate_strain()
			if (Flag_fast == 0):
				tdi_y = gw_y.construct_TDI(t, Orbit)
			else:
				tdi_y = gw_y.construct_TDI_FAST(Orbit)
								
			f_min = np.max([np.min(TDI.freqs), np.min(tdi_y.freqs)])
			f_max = np.min([np.max(TDI.freqs), np.max(tdi_y.freqs)])
			
			lkl_y = data_snr - 2.*np.sum(td.get_TDI_overlap(TDI, tdi_y, f_min, f_max)) \
					+ np.sum(td.get_TDI_overlap(tdi_y, tdi_y, gw_y.hp_wavelet.f_min, gw_y.hp_wavelet.f_max))	
					
			lkl_y *= -0.5
		
		q_yx = gw_prop.pdf(x,y)
		q_xy = gw_prop.pdf(y,x)
	
		lnH = (lkl_y - lkl_x) + (q_xy - q_yx) # Todo: priors have yet to be incorporated
		
		# Make a decision
		u = np.log(ss.uniform.rvs(0,1))
		if (u < np.sum(lnH)):
			lkl_x = np.copy(lkl_y)
			#tdi_x = tdi_y #copy.deepcopy(tdi_y)
			x = y
			accept += 1

		out[i,0]  = np.copy(lkl_x)
		out[i,1:] = np.copy(x)
			
	print("Acceptance rate: {}%\n".format(accept/N*100))
	
	return out,props   
    
