import time
from tqdm import tqdm
from copy import deepcopy

import numpy as np
import scipy.linalg as slin
import scipy.stats as ss

import MCMC_tools as mct

import LISA as l
import TDI as td
import Wavelet as wv
import Glitch as gl

# constants
mHz = 1.0e-3
Week = 7*24*3600 

def get_log_likelihood(self, TDI_data, X_flag=None):
    """ Get the log likelihood for the glitch given the data TDI """
    
    f_min = np.max([TDI_data.f_min, self.Glitch.TDI.f_min])
    f_max = np.min([TDI_data.f_max, self.Glitch.TDI.f_max])
    
    overlap = np.sum(td.get_TDI_overlap(TDI_data, self.Glitch.TDI, f_min, f_max, X_flag))
    
    self.Glitch.calc_snr(X_flag)
    
    self.loglkl = overlap - 0.5*self.Glitch.SNR**2
    
    return
    
    
class Model:
    """ Model class """
    
    def __init__(self, Glitch, Orbit):
        self.Glitch = Glitch
        self.Orbit  = Orbit 
            
    # Methods
    get_loglkl = get_log_likelihood
    
    
def Propose_Parameters(proposal_list, paramsND_X, T):
    """ Take a list of proposals and the current state. Return a new set of parameters """
    
    N_props = len(proposal_list)
    
    # gather the weights
    weights = np.zeros(N_props)
    for i in range(N_props):
        weights[i] = proposal_list[i].weight
    
    # convert weights to a distrigltion
    cumsum = 0
    ls = np.zeros(N_props+1)
    for i in range(N_props):
        cumsum += weights[i]
        ls[i+1] = cumsum

    # find which distrigltion was chosen
    u = ss.uniform.rvs(0,1)
    for i in range(N_props):
        if (ls[0] < u < ls[1]):
            break

    return proposal_list[i].rvs(paramsND_X, T), i
    
    
class Proposal:
    """ GW Glitch Proposal class """
    
    def __init__(self, name, weight):
        self.name = name
        self.weight = weight  


def Fisher_logpdf(self):
    """ This distrigltion is symmetric in location so return a constant """
    
    return 1.0
    
def Fisher_rvs(self, paramsND_X, T=1.0):
    """ Given fisher matrix eigen-bs provide a random sample from its guassian approximation """
    
    e_vals = self.e_vals
    e_vecs = self.e_vecs
    
    e_dir = np.random.choice( range(len(e_vals)) ) # choose and eigen-direction
    e_val = e_vals[e_dir]
        
    jump = ss.norm.rvs(0,1)/np.sqrt(e_val)*e_vecs[:,e_dir]*np.sqrt(T)

    return paramsND_X + jump
    
        
class Proposal_Fisher(Proposal):
    """ GW Glitch Fisher proposal class """
    
    def __init__(self, name, weight, e_vals, e_vecs):
        self.name = name
        self.weight = weight
        self.e_vals = e_vals
        self.e_vecs = e_vecs
        
    logpdf = Fisher_logpdf
    rvs    = Fisher_rvs

   
def DE_logpdf(self):
    """ This distrigltion is symmetric in location so return a constant """
    
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
    
    def __init__(self, name, weight, history, itr):
        self.name = name
        self.weight = weight
        self.history = history
        self.itr = itr
        
    logpdf = DE_logpdf
    rvs    = DE_rvs


def glitch_prior_logpdf(paramsND):
    """ GW Glitch prior log PDF on parameters """
    
    result = 0.0
    
    A     = np.exp(paramsND[wv.IDX_lnA]) 
    f0    = paramsND[wv.IDX_f0]*mHz
    t0    = paramsND[wv.IDX_t0]*Week
    tau   = paramsND[wv.IDX_tau]*Week
    phi0  = paramsND[wv.IDX_phi0]
    
    # wavelet priors
    if (A < 1.0e-25): # TODO: gotta develop and SNR prior, glt for now enforce positivity
        result -= np.inf
    if not (1.0e-4 < f0 < 1.0e-1):
        result -= np.inf
    if not (0.0 < t0 < 1.5*Week):
        result -= np.inf
    if not (0.0 < tau < 1.5*Week):
        result -= np.inf
    if not (0.0 < phi0 < 2*np.pi):
        result -= np.inf
    
    return result


class Prior:
    """ GW glrst prior class """
    
    def __init__(self, func_logpdf, func_rvs=None):
        self.logpdf = func_logpdf  
            
    
class Flag:
    """ Glitch MCMC flag class """
    
    def __init__(self, comp_switch, PT, X):
        self.comp_switch = comp_switch
        self.PT = PT
        self.X = X
        
        
def MCMC_glitch(TDI_data, Orbit, modelX0, seed, N, Burn, Flag):
    """ Perform an MCMC in the glitch model """
    
    comp_id = modelX0.Glitch.comp_id
    
    np.random.seed(seed) 
    X_flag = Flag.X
    
    under_sample = 10
    
    # Setup temperatures for MCMC
    N_temps = 1 
    T_ladder = np.array([1.0])
    if (Flag.PT): 
        data_snr = np.sum(td.get_TDI_overlap(TDI_data, TDI_data, TDI_data.f_min, TDI_data.f_max, X_flag))
        Tmax = data_snr/5**2 # so that hottest chain can see SNR 5 features
        N_temps = 16
        T_ladder = np.geomspace(1., Tmax, N_temps)
            
    who = np.arange(0, N_temps, 1) # ID tags for each chain
    
    # For tracking the ML solution
    max_logL = modelX0.loglkl
    max_model = deepcopy(modelX0)
    
    # Storage of acceptance rates
    acc_cnt = np.zeros(N_temps) 
    acc = np.zeros(N_temps) 
    swap_cnt = np.zeros(N_temps-1)
    swap = np.zeros(N_temps-1)
        
    # MCMC data storage
    chain = np.zeros((len(modelX0.Glitch.paramsND), N, N_temps))
    lkl   = np.zeros((N, N_temps))  
    comp_id_ls   = np.zeros((N, N_temps)) 
    
    # create the modelX list for each temperature
    modelX_list = []
    for i in range(N_temps):
        modelX_list.append(deepcopy(modelX0))
        
    ######## setup the propsals ########
    proposal_list = []    
    modelX0.Glitch.calc_Fish(X_flag)

    # Perform SVD on Fisher matrix to eliminate bad directions
    M,NN = modelX0.Glitch.Fisher.shape
    U,s,Vh = slin.svd(modelX0.Glitch.Fisher)
    Sig = slin.diagsvd(s,M,NN)

    for i in range(len(Sig)):
        if (Sig[i,i] < 0.01):
            Sig[i,i] = 0.

    Fish = U.dot(Sig.dot(Vh))
    e_vals2, e_vecs2 = mct.get_Fisher_Eignen_BS(Fish)
    mask = (e_vals2 > 1.0e-10)
    e_vals2 = e_vals2[mask]
    e_vecs2 = e_vecs2[:,mask]
    #e_vals2, e_vecs2 = mct.get_Fisher_Eignen_BS(modelX0.Glitch.Fisher)
    Prop_Fish2 = Proposal_Fisher("Fisher", 0.8, e_vals2, e_vecs2)
    proposal_list.append(Prop_Fish2) 
    
    # DE proposal
    N_hist = 1000
    history = np.zeros((N_hist, len(modelX0.Glitch.paramsND)))
    Prop_Hist = Proposal_DE("DiffEv", 0.2, history, itr=0)
    proposal_list.append(Prop_Hist)    
    
    # set up the priors
    prior = Prior(glitch_prior_logpdf) 
    
#     if (N_temps != 1):
#         T_ladder[-1] = 1.0e6 # for adaptive T scheme, make largest T inf
#         t0 = 1000
#         nu = t0/10

    for i in tqdm(range(-Burn, N)):
        if (Flag.PT):
            # decide whether to swap temperatures or not...
            u = ss.uniform.rvs(0,1)
            
            # adjust the temperature ladder
#             if (i > 100 and N_temps != 1):
#                 T_old = deepcopy(T_ladder)
#                 A = swap/swap_cnt
#                 mask = np.isnan(A)
#                 A[mask] = 0
#                 
#                 for j in range(1, N_temps-1):
#                     kappa = t0/(t0+i)/nu
#                     T_ladder[j] = T_ladder[j-1] + (T_old[j] - T_old[j-1])*np.exp(kappa*(A[j-1] - A[j]))
#                     
#                     if (T_ladder[j]/T_ladder[j-1] < 1.1):
#                         T_ladder[j] = 1.1*T_ladder[j-1]
                        
        else:
            u = 1 # i.e. force the use of normal MCMC jumps
            
            
        if (u<0.5):
            # propose a swap in temperatures
            k = int((N_temps-1)*ss.uniform.rvs(0,1))
            swap_cnt[k] += 1
            whoA = np.copy(who[k])
            whoB = np.copy(who[k+1])
            
            beta  = (modelX_list[whoA].loglkl - modelX_list[whoB].loglkl)/T_ladder[k+1]
            beta -= (modelX_list[whoA].loglkl - modelX_list[whoB].loglkl)/T_ladder[k]
            
            u = np.log(ss.uniform.rvs(0,1))
            if (u < beta):
                swap[k] += 1
                hold = deepcopy(whoA)
                who[k]   = deepcopy(whoB)
                who[k+1] = deepcopy(hold)
                
            for j in range(N_temps):
                ID = np.copy(who[j])
                if (i>=0):
                    chain[:,i,j] = deepcopy(modelX_list[ID].Glitch.paramsND)
                    lkl[i,j]     = deepcopy(modelX_list[ID].loglkl)
                    comp_id_ls[i,j] = modelX_list[ID].Glitch.comp_id

            Prop_Hist.itr += 1 # add on an iteration
            Prop_Hist.history[i%N_hist,:] = deepcopy(modelX_list[who[0]].Glitch.paramsND)
            
        else:
            acc_cnt += 1
            for j in range(N_temps):

                ID = np.copy(who[j])
                T  = np.copy(T_ladder[j])
            
                paramsND_Y, propID = Propose_Parameters(proposal_list, deepcopy(modelX_list[ID].Glitch.paramsND), T)
                # wrap the phase
                if (paramsND_Y[wv.IDX_phi0] > 2*np.pi):
                    paramsND_Y[wv.IDX_phi0] -= 2*np.pi
                if (paramsND_Y[wv.IDX_phi0] < 0):
                    paramsND_Y[wv.IDX_phi0] += 2*np.pi
                    
                if (Flag.comp_switch == 1):#and i%4==0):
                    comp_id_Y = np.random.choice(range(12))
                    
                    # random phi0 flick for the perfect degeneracies in acceleration glitches
                    u = ss.uniform.rvs(0,1)
                    if (u<0.1):
                        paramsND_Y[wv.IDX_phi0] += np.pi + ss.norm.rvs(0,0.05)
                    elif (0.1 < u < 0.2):
                        paramsND_Y[wv.IDX_phi0] -= np.pi + ss.norm.rvs(0,0.05)
                    
                else:
                    comp_id_Y = np.copy(modelX_list[ID].Glitch.comp_id)
                
                if (prior.logpdf(paramsND_Y) != -np.inf): # check that parameters are acceptable
                
                    # generate Y Glitch and calculate its signal
                    Glitch_Y = gl.Glitch(deepcopy(paramsND_Y), Orbit, np.copy(comp_id_Y))
                    Glitch_Y.calc_TDI()

                    modelY = Model(Glitch_Y, Orbit) 
                    modelY.get_loglkl(TDI_data, X_flag)
                
                    # calculate the log PDF of the priors
                    logH  = prior.logpdf(modelY.Glitch.paramsND) - prior.logpdf(modelX_list[ID].Glitch.paramsND)
                    logH += proposal_list[propID].logpdf() - proposal_list[propID].logpdf() # doesn't matter
                    logH += (modelY.loglkl - modelX_list[ID].loglkl)/T
                
                    u = np.log(ss.uniform.rvs(0,1))
                
                    if (u < logH): # accept the move
                        acc[j] += 1
                        modelX_list[ID] = deepcopy(modelY)
                        modelX_list[ID].Glitch.comp_id = np.copy(comp_id_Y)
                        
                        if (j == 0 and modelX_list[ID].loglkl > max_logL):
                            max_logL  = deepcopy(modelX_list[ID].loglkl)
                            max_model = deepcopy(modelX_list[ID])
                
                if (i>=0):        
                    chain[:,i,j] = deepcopy(modelX_list[ID].Glitch.paramsND)
                    lkl[i,j]     = deepcopy(modelX_list[ID].loglkl)
                    comp_id_ls[i,j] = deepcopy(modelX_list[ID].Glitch.comp_id)
                  
                # update the fisher matrix  
#                 if (j == 0  and i%20 == 0):
#                     name = 'huh'
#                     k = -1
#                     while (name != 'Fisher'):
#                         k += 1
#                         name = proposal_list[k].name
#                     modelX_list[who[j]].Glitch.calc_Fish(X_flag)
#                     M,NN = modelX0.Glitch.Fisher.shape
#                     U,s,Vh = slin.svd(modelX0.Glitch.Fisher)
#                     Sig = slin.diagsvd(s,M,NN)
# 
#                     for i in range(len(Sig)):
#                         if (Sig[i,i] < 0.01):
#                             Sig[i,i] = 0.    
# 
#                     Fish = U.dot(Sig.dot(Vh))
#                     e_vals2, e_vecs2 = mct.get_Fisher_Eignen_BS(Fish)
#                     mask = (e_vals2 > 1.0e-10)
#                     e_vals2 = e_vals2[mask]
#                     e_vecs2 = e_vecs2[:,mask]
#                     proposal_list[k] = Proposal_Fisher("Fisher", 0.8, e_vals2, e_vecs2)                   
            
                if (j==0):
                    Prop_Hist.itr += 1 # add on an iteration
                    Prop_Hist.history[i%N_hist,:] = deepcopy(modelX_list[ID].Glitch.paramsND)

    print("Temperature Ladder...... {}".format(T_ladder))
    print("acceptance rate......... {}%".format(acc/acc_cnt*100))
    print("swap acceptance rate.... {}%".format(swap/swap_cnt*100))
    print("max logL................ {}" .format(max_logL))          
    
    return lkl[::under_sample,:], chain[:,::under_sample,:], comp_id_ls[::under_sample,:], max_model    
    