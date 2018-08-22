import numpy as np
import scipy.stats as ss
import scipy.linalg as slin

from copy import deepcopy
import functools
import time
from tqdm import tqdm


import LISA as l

import Wavelet as wv
import Burst as bu
import TDI as td

import MCMC_tools as mct

# constants
mHz = 1.0e-3
Week = 3600.*24*7

def calculate_log_likelihood(self, TDI_data, X_flag=None):
    """ 
        Calculate the log-likelihood of the GW burst model 
        Assumes that SNR has been calculated 
    """
    
    # Find the bounds of the NWIP integrals
    f_min = np.max([TDI_data.f_min, self.Burst.TDI.f_min])
    f_max = np.min([TDI_data.f_max, self.Burst.TDI.f_max])
    
    overlap = np.sum(td.get_TDI_overlap(TDI_data, self.Burst.TDI, f_min, f_max, X_flag))
    
    self.Burst.calc_snr(X_flag)
    
    self.loglkl = overlap - 0.5*self.Burst.SNR**2
    
    return
    

class Model:
    """ GW burst model class """
    
    def __init__(self, Burst, Orbit):
        self.Burst = Burst
        self.Orbit = Orbit
        
    # Methods
    get_loglkl = calculate_log_likelihood
    
    
def burst_prior_logpdf(paramsND):
    """ GW Burst prior log PDF on parameters """
    
    result = 0.0
    
    A     = np.exp(paramsND[wv.IDX_lnA]) 
    f0    = paramsND[wv.IDX_f0]*mHz
    t0    = paramsND[wv.IDX_t0]*Week
    tau   = paramsND[wv.IDX_tau]*Week
    phi0  = paramsND[wv.IDX_phi0]
    
    cost  = paramsND[bu.IDX_cost]
    phi   = paramsND[bu.IDX_phi]
    psi   = paramsND[bu.IDX_psi]
    ellip = paramsND[bu.IDX_ellip]
    
    # wavelet priors
    if (A < 0.0): # TODO: gotta develop and SNR prior, but for now enforce positivity
        result -= np.inf
    if not (1.0e-5 < f0 < 1.0e-1):
        result -= np.inf
    if not (0.0 < t0 < 1.5*Week):
        result -= np.inf
    if (tau < 0.0):
        result -= np.inf
    if not (0.0 < phi0 < 2*np.pi):
        result -= np.inf
       
    # burst specific priors
    if not (-1 < cost < 1):
        result -= np.inf
    if not (0 < phi < 2*np.pi):
        result -= np.inf
    if not (0 < psi < np.pi):
        result -= np.inf
    if not (0 < ellip < 1):
        result -= np.inf
    
    return result


class Prior:
    """ GW burst prior class """
    
    def __init__(self, func_logpdf, func_rvs=None):
        self.logpdf = func_logpdf
        

class Proposal:
    """ GW burst Proposal class """
    
    def __init__(self, weight):
        self.weight = weight       
    

def Fisher_logpdf(self):
    """ This distribution is symmetric in location so return a constant """
    
    return 1.0
    
def Fisher_rvs(self, paramsND_X, T=1.0):
    """ Given fisher matrix eigen-bs provide a random sample from its guassian approximation """
    
    e_vals = self.e_vals
    e_vecs = self.e_vecs
    
    e_dir = np.random.choice( range(len(e_vals)) ) # choose and eigen-direction
    e_val = e_vals[e_dir]
    if (e_val < 0):
        e_val = 0.1
    jump = ss.norm.rvs(0,1)/np.sqrt(e_val)*e_vecs[:,e_dir]*np.sqrt(T)

    return paramsND_X + jump
    
        
class Proposal_Fisher(Proposal):
    """ GW burst Fisher proposal class """
    
    def __init__(self, weight, e_vals, e_vecs):
        self.weight = weight
        self.e_vals = e_vals
        self.e_vecs = e_vecs
        
    logpdf = Fisher_logpdf
    rvs    = Fisher_rvs
    
    
def DE_logpdf(self):
    """ This distribution is symmetric in location so return a constant """
    
    return 1.0
    
def DE_rvs(self, paramsND_X, T=1.0):
    """ Propose a sample along history generated vector direction """
    
    hist = self.history
    N_hist = len(hist)
    itr  = int(self.itr)%N_hist
   
    
    if (itr < 2):
        return paramsND_X
    else:
        j = int(itr*ss.uniform.rvs(0,1))
        k = deepcopy(j) 
        while (k==j):
            k = int(itr*ss.uniform.rvs(0,1))
        
        alpha = 1.0
        beta  = ss.uniform.rvs(0,1)
        
        if (beta < 0.9):
            alpha = ss.norm.rvs(0,1)
        
        direction = (hist[j] - hist[k])

    return paramsND_X + alpha*direction*np.sqrt(T)
    
    
class Proposal_DE(Proposal):
    """ Differential Evolution Proposal Class """
    
    def __init__(self, weight, history, itr):
        self.weight = weight
        self.history = history
        self.itr = itr
        
    logpdf = DE_logpdf
    rvs    = DE_rvs
    
    
def SkyRef_logpdf(self):
    """ This distribution is symmetric in location so return a constant """
    
    return 1.0
    
def SkyRef_rvs(self, paramsND_X, T=1.0):
    """ Given fisher matrix eigen-bs provide a random sample from its guassian approximation """
    
    theta = np.arccos(paramsND_X[IDX_cost])
    phi   = paramsND_X[idx_phi]
    
    n_hat = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    n_hatp = n_hat - 2*np.dot(n_hat, self.p_hat)*self.p_hat
    
    phip_mean   = np.arctan2(n_hatp[1]/np.sin(thetap), n_hatp[0]/np.sin(thetap))

    costp = ss.norm.rvs(n_hatp[-1], 0.1)
    phip = ss.norm.rvs(phip_mean, 0.1)
    
    paramsND_Y = deepcopy(paramND_X)
    paramsND_Y[IDX_cost] = costp
    paramsND_Y[IDX_phi] = phip

    return paramsND_Y

    
class Proposal_SkyRef(Proposal):
    """ Sky reflection: across plane of detector """
    
    def __init__(self, weight, p_hat):
        self.weight = weight
        self.p_hat  = p_hat # vector defining plane of LISA
        
    logpdf = SkyRef_logpdf
    rvs    = SkyRef_rvs
   
   
def Propose_Parameters(proposal_list, paramsND_X, T):
    """ Take a list of proposals and the current state. Return a new set of parameters """
    
    N_props = len(proposal_list)
    
    # gather the weights
    weights = np.zeros(N_props)
    for i in range(N_props):
        weights[i] = proposal_list[i].weight
    
    # convert weights to a distribution
    cumsum = 0
    ls = np.zeros(N_props+1)
    for i in range(N_props):
        cumsum += weights[i]
        ls[i+1] = cumsum

    # find which distribution was chosen
    u = ss.uniform.rvs(0,1)
    for i in range(N_props):
        if (ls[0] < u < ls[1]):
            break

    return proposal_list[i].rvs(paramsND_X, T), i
    

class Flag:
    """ GW burst flag class """
    
    def __init__(self, PT, X):
        self.PT = PT
        self.X  = X # X=1--use X-channel only

    
def Burst_MCMC(TDI_data, Orbit, modelX0, seed, N, Flag):
    """ Burst MCMC """
    
    X_flag = Flag.X
    
    under_sample = 50
    
    # Setup temperatures for MCMC
    N_temps = 1 # will be changing later
    T_ladder = np.array([1.0])
    if (Flag.PT): 
        data_snr = np.sum(td.get_TDI_overlap(TDI_data, TDI_data, TDI_data.f_min, TDI_data.f_max, X_flag))
        Tmax = data_snr/5**2 # so that hottest chain can see SNR 5 features
        N_temps = 4
        T_ladder = np.geomspace(1., Tmax, N_temps) # Todo: fix this..., just a choice for now...
        
    who = np.arange(0, N_temps, 1)
    
    # keep track of the ML solution
    max_logL = 1.0e-30 
    max_model = deepcopy(modelX0)
    
    np.random.seed(seed) # set Numpy's seed for reproducibility 
    
    # Storage of acceptance rates
    acc_cnt = 0
    acc = np.zeros(N_temps) 
    swap_cnt = np.zeros(N_temps-1)
    swap = np.zeros(N_temps-1)
    
    # MCMC data storage
    chain = np.zeros((len(modelX0.Burst.paramsND), N, N_temps))
    lkl   = np.zeros((N, N_temps))
    #props = np.zeros((len(modelX0.Burst.paramsND), N, N_temps)) # Todo: should be a flag for this, not really needed
    
    # create the modelX list for each temperature
    modelX_list = []
    for i in range(N_temps):
        modelX_list.append(deepcopy(modelX0))
        
        
    # setup the propsals
    proposal_list = []
    
    # calculate the fisher matrix and eigen-BS off initial states
    modelX0.Burst.calc_Fish(X_flag)
#     e_vals1, e_vecs1 = mct.get_Fisher_Eignen_BS(modelX0.Burst.Fisher)
#     mask = (e_vals1 > 0) # remove pesky eigen-vectors
#     e_vals1 = e_vals1[mask]
#     e_vecs1 = e_vecs1[:,mask]  
#     
#     Prop_Fish1 = Proposal_Fisher(0.4, e_vals1, e_vecs1)  
#     proposal_list.append(Prop_Fish1) 
    
    # Perform SVD on Fisher matrix to eliminate bad directions
    M,NN = modelX0.Burst.Fisher.shape
    U,s,Vh = slin.svd(modelX0.Burst.Fisher)
    Sig = slin.diagsvd(s,M,NN)

    for i in range(len(Sig)):
        if (Sig[i,i] < 0.01):
            Sig[i,i] = 0.

    Fish = U.dot(Sig.dot(Vh))
    e_vals2, e_vecs2 = mct.get_Fisher_Eignen_BS(Fish)
    mask = (e_vals2 > 1.0e-10)
    e_vals2 = e_vals2[mask]
    e_vecs2 = e_vecs2[:,mask]

    Prop_Fish2 = Proposal_Fisher(0.6, e_vals2, e_vecs2)
    proposal_list.append(Prop_Fish2) 
    
    
    # DE proposal
    N_hist = 2000
    history = np.zeros((N_hist, len(modelX0.Burst.paramsND)))
    Prop_Hist = Proposal_DE(0.2, history, itr=0)
    proposal_list.append(Prop_Hist) 
    
    # sky reflection proposal
    ti = 0
    t = np.arange(0, Orbit.Tobs, Orbit.dt)
    Orbit.x = Orbit.get_orbit(t)
    Orbit.rij = Orbit.get_seps(t, Orbit.x)
#    print(
    r12 = Orbit.rij[:,0,1,ti]
    r13 = Orbit.rij[:,0,2,ti]

    plane_vec = np.cross(r12, r13)
    plane_vec /= np.sqrt(np.dot(plane_vec, plane_vec))
    
    PropSky = Proposal_SkyRef(0.2, plane_vec)
    
    # set up the priors
    prior = Prior(burst_prior_logpdf)
    
    if (N_temps != 1):
        T_ladder[-1] = 1.0e6 # for adaptive T scheme, make largest T inf
        t0 = 10000
        nu = t0/10
        
    for i in tqdm(range(N)):
        if (Flag.PT):
            # decide whether to swap temperatures or not...
            u = ss.uniform.rvs(0,1)
            
            # adjust the temperature ladder
            if (i > 100 and N_temps != 1):
                T_old = deepcopy(T_ladder)
                A = swap/swap_cnt
                mask = np.isnan(A)
                A[mask] = 0
                
                for j in range(1, N_temps-1):
                    kappa = t0/(t0+i)/nu
                    T_ladder[j] = T_ladder[j-1] + (T_old[j] - T_old[j-1])*np.exp(kappa*(A[j-1] - A[j]))
                    
                    if (T_ladder[j]/T_ladder[j-1] < 1.1):
                        T_ladder[j] = 1.1*T_ladder[j-1]
                                    
        else:
            u = 1 # i.e. force the use of normal MCMC jumps
            
        if (u<0.5):
            # propose a swap in temperatures
            k = int((N_temps-1)*ss.uniform.rvs(0,1))
            swap_cnt[k] += 1
            whoA = deepcopy(who[k])
            whoB = deepcopy(who[k+1])
            
            beta  = (modelX_list[whoA].loglkl - modelX_list[whoB].loglkl)/T_ladder[k+1]
            beta -= (modelX_list[whoA].loglkl - modelX_list[whoB].loglkl)/T_ladder[k]
            
            u = np.log(ss.uniform.rvs(0,1))
            if (u < beta):
                swap[k] += 1
                hold = deepcopy(whoA)
                whoA = deepcopy(whoB)
                whob = deepcopy(hold)
                
            for j in range(N_temps):
                ID = deepcopy(who[j])
                chain[:,i,j] = deepcopy(modelX_list[ID].Burst.paramsND)
                lkl[i,j]     = deepcopy(modelX_list[ID].loglkl)

            Prop_Hist.itr += 1 # add on an iteration
            Prop_Hist.history[i%N_hist,:] = modelX_list[who[0]].Burst.paramsND
            
        else:
            for j in range(N_temps):
                acc_cnt += 1
            
                ID = who[j] 
                T  = T_ladder[j]

                paramsND_Y, propID = Propose_Parameters(proposal_list, modelX_list[ID].Burst.paramsND, T)
                
                if (prior.logpdf(paramsND_Y) != -np.inf): # check that parameters are acceptable
                
                    # generate Y burst and calculate its signal
                    Burst_Y = bu.Burst(deepcopy(paramsND_Y), Orbit)
                    Burst_Y.construct_detector_tensor()
                    Burst_Y.calculate_strain()
                    Burst_Y.TDI = Burst_Y.construct_TDI(Orbit)

                    modelY = Model(Burst_Y, Orbit) 
                    modelY.get_loglkl(TDI_data, X_flag)
                
                    # calculate the log PDF of the priors
                    logH  = prior.logpdf(modelY.Burst.paramsND) - prior.logpdf(modelX_list[ID].Burst.paramsND)
                    logH += proposal_list[propID].logpdf() - proposal_list[propID].logpdf() # doesn't matter
                    logH += (modelY.loglkl - modelX_list[ID].loglkl)/T
                
                    u = np.log(ss.uniform.rvs(0,1))
                
                    if (u < logH): # accept the move
                        acc[j] += 1
                        modelX_list[ID] = deepcopy(modelY)
                    
                        if (j == 0 and modelX_list[ID].loglkl > max_logL):
                            max_logL  = deepcopy(modelX_list[ID].loglkl)
                            max_model = deepcopy(modelX_list[ID])
                        
                chain[:,i,j] = deepcopy(modelX_list[ID].Burst.paramsND)
                lkl[i,j]     = deepcopy(modelX_list[ID].loglkl)
            
                if (j==0):
                    Prop_Hist.itr += 1 # add on an iteration
                    Prop_Hist.history[i%N_hist,:] = modelX_list[ID].Burst.paramsND

    print("Temperature Ladder...... {}".format(T_ladder))
    print("acceptance rate......... {}%".format(acc/acc_cnt*100))
    print("swap acceptance rate.... {}%".format(swap/swap_cnt*100))
    print("max logL................ {}" .format(max_logL))          
    
    return lkl[::under_sample,:], chain[:,::under_sample,:], max_model
    
















       