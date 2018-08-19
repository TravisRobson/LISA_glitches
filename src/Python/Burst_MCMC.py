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

def calculate_log_likelihood(self, TDI_data):
    """ 
        Calculate the log-likelihood of the GW burst model 
        Assumes that SNR has been calculated 
    """
    
    # Find the bounds of the NWIP integrals
    f_min = np.max([TDI_data.f_min, self.Burst.TDI.f_min])
    f_max = np.min([TDI_data.f_max, self.Burst.TDI.f_max])
    
    overlap = np.sum(td.get_TDI_overlap(TDI_data, self.Burst.TDI, f_min, f_max))
    
    self.Burst.calc_snr()
    
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
    if not (0 < psi < np.pi/2):
        result -= np.inf
    if not (-1 < ellip < 1):
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
    
#     # methods
#     logpdf = func_logpdf
#     rvs    = func_rvs
    

def Fisher_logpdf(self):
    """ This distribution is symmetric in location so return a constant """
    
    return 1.0
    
def Fisher_rvs(self, paramsND_X):
    """ Given fisher matrix eigen-bs provide a random sample from its guassian approximation """
    
    e_vals = self.e_vals
    e_vecs = self.e_vecs
    
    e_dir = np.random.choice( range(len(e_vals)) ) # choose and eigen-direction
    e_val = e_vals[e_dir]
    if (e_val < 0):
        e_val = 0.1
    jump = ss.norm.rvs(0,1)/np.sqrt(e_val)*e_vecs[:,e_dir]

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
    
def DE_rvs(self, paramsND_X):
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
        
    return paramsND_X + alpha*direction
    
    
class Proposal_DE(Proposal):
    """ Differential Evolution Proposal Class """
    
    def __init__(self, weight, history, itr):
        self.weight = weight
        self.history = history
        self.itr = itr
        
    logpdf = DE_logpdf
    rvs    = DE_rvs
   
   
def Propose_Parameters(proposal_list, paramsND_X):
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

    return proposal_list[i].rvs(paramsND_X), i
    

class Flag:
    """ GW burst flag class """
    
    def __init__(self, PT):
        self.PT = PT

    
def Burst_MCMC(TDI_data, Orbit, modelX0, seed, N, Flag):
    """ Burst MCMC """
    
    # Setup temperatures for MCMC
    N_temps = 1 # will be changing later
    T_ladder = np.array([1.0])
    if (Flag.PT): 
        data_snr = np.sum(td.get_TDI_overlap(TDI_data, TDI_data, TDI_data.f_min, TDI_data.f_max))
        Tmax = data_snr/1 # so that hottest chain can see SNR 5 features
        N_temps = 2
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
    props = np.zeros((len(modelX0.Burst.paramsND), N, N_temps)) # Todo: should be a flag for this, not really needed
    
    # create the modelX list for each temperature
    modelX_list = []
    for i in range(N_temps):
        modelX_list.append(deepcopy(modelX0))
        
    # calculate the fisher matrix and eigen-BS off initial states
    modelX0.Burst.calc_Fish()
    e_vals, e_vecs = mct.get_Fisher_Eignen_BS(modelX0.Burst.Fisher)

    # setup the propsals
    proposal_list = []
    
    Prop_Fish = Proposal_Fisher(0.7, e_vals, e_vecs)
    
    proposal_list.append(Prop_Fish) 
    
    # DE proposal
    N_hist = 1000
    history = np.zeros((N_hist, len(modelX0.Burst.paramsND)))
    Prop_Hist = Proposal_DE(0.3, history, itr=0)
    proposal_list.append(Prop_Hist) 
    
    # set up the priors
    prior = Prior(burst_prior_logpdf)
    
    #for i in range(N):
    for i in tqdm(range(N)):

        if (Flag.PT):
            # decide whether to swap temperatures or not...
            u = ss.uniform.rvs(0,1)
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
            
                paramsND_Y, propID = Propose_Parameters(proposal_list, modelX_list[ID].Burst.paramsND)

                if (prior.logpdf(paramsND_Y) != -np.inf): # check that parameters are acceptable
                    # generate Y burst and calculate its signal
                    Burst_Y = bu.Burst(deepcopy(paramsND_Y), Orbit)
                    Burst_Y.construct_detector_tensor()
                    Burst_Y.calculate_strain()
                    Burst_Y.TDI = Burst_Y.construct_TDI(Orbit)
                
                    modelY = Model(Burst_Y, Orbit) 
                    modelY.get_loglkl(TDI_data)
                
                    # calculate the log PDF of the priors
                    logH  = prior.logpdf(modelY.Burst.paramsND) - prior.logpdf(modelX_list[ID].Burst.paramsND)
                    logH += proposal_list[propID].logpdf() - proposal_list[propID].logpdf() # doesn't matter
                    logH += modelY.loglkl - modelX_list[ID].loglkl
                
                    u = np.log(ss.uniform.rvs(0,1))
                
                    if (u < logH): # accept the move
                        acc += 1
                        modelX_list[ID] = deepcopy(modelY)
                    
                        if (j == 0 and modelX_list[ID].loglkl > max_logL):
                            max_logL  = deepcopy(modelX_list[ID].loglkl)
                            max_model = deepcopy(modelX_list[ID])
                        
                chain[:,i,j] = deepcopy(modelX_list[ID].Burst.paramsND)
                lkl[i,j]     = deepcopy(modelX_list[ID].loglkl)
            
                if (j==0):
                    Prop_Hist.itr += 1 # add on an iteration
                    Prop_Hist.history[i%N_hist,:] = modelX_list[ID].Burst.paramsND

    print("acceptance rate......... {}%".format(acc/acc_cnt*100))
    print("swap acceptance rate.... {}%".format(swap/swap_cnt*100))
    print("max logL................ {}" .format(max_logL))          
    
    return lkl, chain, props, max_model
    
















       