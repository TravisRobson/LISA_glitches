import numpy as np
import scipy.stats as ss
import scipy.linalg as slin

from copy import deepcopy
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
    """ Calculate the log-likelihood of the GW burst model """
    
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
    
    
def burst_prior_logpdf(paramsND, T=1.0):
    """ GW Burst prior log PDF on parameters """
    
    result = 0.0
        
    lnA   = paramsND[wv.IDX_lnA]
    f0    = paramsND[wv.IDX_f0]*mHz
    t0    = paramsND[wv.IDX_t0]*Week
    tau   = paramsND[wv.IDX_tau]*Week
    phi0  = paramsND[wv.IDX_phi0]
    
    cost  = paramsND[bu.IDX_cost]
    phi   = paramsND[bu.IDX_phi]
    psi   = paramsND[bu.IDX_psi]
    ellip = paramsND[bu.IDX_ellip]
    
    # wavelet priors
    if not (np.log(1.0e-25) < lnA < np.log(1.0e-15)): # TODO: gotta develop and SNR prior, but for now enforce positivity
        result -= np.inf
    if not (0.1e-3 < f0 < 40.e-3):
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
    

def Fisher_logpdf(self, paramsND_X=None):
    """ This distribution is symmetric in location so return a constant """
    
    return 1.0
    
def Fisher_rvs(self, paramsND_X, T=1.0):
    """ Given fisher matrix eigen-bs provide a random sample from its guassian approximation """
    
    T_idx = np.where(T==self.T_ls)[0][0]

    e_vals = self.e_vals[T_idx]
    e_vecs = self.e_vecs[T_idx]
    
    e_dir = np.random.choice( range(len(e_vals)) ) # choose and eigen-direction
    e_val = e_vals[e_dir]

    jump = ss.norm.rvs(0,1)/np.sqrt(e_val)*e_vecs[:,e_dir]*np.sqrt(T)

    return paramsND_X + jump
    
        
class Proposal_Fisher(Proposal):
    """ GW burst Fisher proposal class """
    
    def __init__(self, name, weight, e_vals_ls, e_vecs_ls, T_ls):
        self.name = name
        self.weight = weight
        self.e_vals = e_vals_ls
        self.e_vecs = e_vecs_ls
        self.N_temps = len(e_vals_ls)
        self.T_ls = T_ls
        
    logpdf = Fisher_logpdf
    rvs    = Fisher_rvs
    
    
def DE_logpdf(self, paramsND_X=None):
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
    
    def __init__(self, name,weight, history, itr):
        self.name = name
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

    theta = np.arccos(paramsND_X[bu.IDX_cost])
    phi   = paramsND_X[bu.IDX_phi]
    
    n_hat = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    n_hatp = n_hat - 2*np.dot(n_hat, self.p_hat)*self.p_hat
    
    thetap = np.arccos(n_hatp[-1])
    phip_mean   = np.arctan2(n_hatp[1]/np.sin(thetap), n_hatp[0]/np.sin(thetap))

    costp = ss.norm.rvs(n_hatp[-1], 0.05)
    phip = ss.norm.rvs(phip_mean, 0.05)
        
    paramsND_Y = deepcopy(paramsND_X)
    paramsND_Y[bu.IDX_cost] = costp
    paramsND_Y[bu.IDX_phi] = phip

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

    
def Burst_MCMC(TDI_data, Orbit, modelX0, seed, N, N_Burn, Flag):
    """ Burst MCMC """
    
    X_flag = Flag.X
    
    under_sample = 10
    
    # Setup temperatures for MCMC
    N_temps = 1 # will be changing later
    T_ladder = np.array([1.0])
    if (Flag.PT): 
        data_snr = np.sum(td.get_TDI_overlap(TDI_data, TDI_data, TDI_data.f_min, TDI_data.f_max, X_flag))
        Tmax = data_snr/5**2 # so that hottest chain can see SNR 5 features
        N_temps = 15 #30
        T_ladder = np.geomspace(1., Tmax, N_temps) # Todo: fix this..., just a choice for now...
        
    who = np.arange(0, N_temps, 1)
    
    # keep track of the ML solution
    max_logL = 1.0e-30 
    max_model = deepcopy(modelX0)
    
    np.random.seed(seed) # set Numpy's seed for reproducibility 
    
    # Storage of acceptance rates
    acc_cnt = np.zeros(N_temps) 
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
#     Prop_Fish1 = Proposal_Fisher(0.6, e_vals1, e_vecs1)  
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
    
    e_vals_ls = [e_vals2]*N_temps
    e_vecs_ls = [e_vecs2]*N_temps
    #Prop_Fish2 = Proposal_Fisher("Fisher", 0.8, e_vals2, e_vecs2)
    Prop_Fish2 = Proposal_Fisher("Fisher", 0.8, e_vals_ls, e_vecs_ls, T_ladder)
    proposal_list.append(Prop_Fish2) 
    
    
    # DE proposal
    N_hist = 1000
    history = np.zeros((N_hist, len(modelX0.Burst.paramsND)))
    Prop_Hist = Proposal_DE("DiffEv", 0.2, history, itr=0)
    proposal_list.append(Prop_Hist) 
    
    # sky reflection proposal
#     t = np.arange(0, Orbit.Tobs, Orbit.dt)
#     mask = (t-modelX0.Burst.t0 < Orbit.dt) & (modelX0.Burst.t0-t < Orbit.dt)
#     ti = t[mask]
#     ti_idx = np.where(t == ti)[0][0]
#     Orbit.x = Orbit.get_orbit(t)
#     Orbit.rij = Orbit.get_seps(t, Orbit.x)
#     r12 = Orbit.rij[:,0,1,ti_idx]
#     r13 = Orbit.rij[:,0,2,ti_idx]
# 
#     plane_vec = np.cross(r12, r13)
#     plane_vec /= np.sqrt(np.dot(plane_vec, plane_vec))
#     
#     Prop_Sky = Proposal_SkyRef(0.0, plane_vec)
#     proposal_list.append(Prop_Sky) 
    
    # set up the priors
    prior = Prior(burst_prior_logpdf)
    
#     if (N_temps != 1):
#         T_ladder[-1] = 1.0e6 # for adaptive T scheme, make largest T inf
#         t0 = 10000
#         nu = 10
        
    T_chain = np.zeros((N, N_temps))    
    
    for i in tqdm(range(-N_Burn, N)):
        if (Flag.PT):
            # decide whether to swap temperatures or not...
            u = ss.uniform.rvs(0,1)
            
            # adjust the temperature ladder
#             if (i > 0 and N_temps != 1):
#                 T_old = np.copy(T_ladder)
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
#                                     
        else:
            u = 1 # i.e. force the use of normal MCMC jumps
            
        if (u<0.3):
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
                who[k]   = deepcopy(whoB)
                who[k+1] = deepcopy(hold)
                
            for j in range(N_temps):
                ID = deepcopy(who[j])
                if (i>=0):
                    chain[:,i,j] = deepcopy(modelX_list[ID].Burst.paramsND)
                    lkl[i,j]     = deepcopy(modelX_list[ID].loglkl)

            Prop_Hist.itr += 1 # add on an iteration
            Prop_Hist.history[i%N_hist,:] = deepcopy(modelX_list[who[0]].Burst.paramsND)
            
        else:
            acc_cnt += 1
            for j in range(N_temps):
            
                ID = np.copy(who[j])
                T  = np.copy(T_ladder[j])

                paramsND_Y, propID = Propose_Parameters(proposal_list, deepcopy(modelX_list[ID].Burst.paramsND), T)
                
                u = ss.uniform.rvs(0,1)
                if (u<0.05):
#                     paramsND_Y[bu.IDX_phi]   = ss.uniform.rvs(0,2*np.pi)  
#                     paramsND_Y[bu.IDX_cost]  = ss.uniform.rvs(-1,2)  
#                     paramsND_Y[bu.IDX_psi]   = ss.uniform.rvs(0,np.pi)  
#                     paramsND_Y[bu.IDX_ellip] = ss.uniform.rvs(0,1) 
                    paramsND_Y[bu.IDX_phi]   += ss.norm.rvs(0,2*np.pi*0.05)  
                    paramsND_Y[bu.IDX_cost]  += ss.norm.rvs(0, 0.05*2)  
                    paramsND_Y[bu.IDX_psi]   += ss.norm.rvs(0, 0.05*np.pi)  
                    paramsND_Y[bu.IDX_ellip] += ss.norm.rvs(0, 0.05) 
                    
                u = ss.uniform.rvs(0,1)
                if (u<0.05):
                    paramsND_Y[wv.IDX_phi0] += np.pi + ss.norm.rvs(0,0.1)
                if (0.05<u<0.10):
                    paramsND_Y[wv.IDX_phi0] -= np.pi + ss.norm.rvs(0,0.1)                      
                    
                # wrap the phase
                while (paramsND_Y[wv.IDX_phi0] > 2*np.pi):
                    paramsND_Y[wv.IDX_phi0] -= 2*np.pi
                while (paramsND_Y[wv.IDX_phi0] < 0):
                    paramsND_Y[wv.IDX_phi0] += 2*np.pi       
                #paramsND_Y[bu.IDX_ellip] = 0.3 # fix this as a test of degeneracy
                if (prior.logpdf(paramsND_Y) > -np.inf): # check that parameters are acceptable

                    # generate Y burst and calculate its signal
                    Burst_Y = bu.Burst(deepcopy(paramsND_Y), Orbit)
                    Burst_Y.construct_detector_tensor()
                    Burst_Y.calculate_strain()
                    Burst_Y.TDI = Burst_Y.construct_TDI(Orbit)

                    modelY = Model(Burst_Y, Orbit) 
                    modelY.get_loglkl(TDI_data, X_flag)
                    if not (np.isnan(modelY.loglkl)):
                
                        # calculate the log PDF of the priors
                        logH  = prior.logpdf(modelY.Burst.paramsND) - prior.logpdf(modelX_list[ID].Burst.paramsND)
                        logH += proposal_list[propID].logpdf(paramsND_Y) - proposal_list[propID].logpdf(modelX_list[ID].Burst.paramsND) # doesn't matter
                        logH += (modelY.loglkl - modelX_list[ID].loglkl)/T
                
                        u = np.log(ss.uniform.rvs(0,1))
                
                        if (u < logH): # accept the move
                            acc[j] += 1
                            modelX_list[ID] = deepcopy(modelY)
                    
                            if (j == 0 and modelX_list[ID].loglkl > max_logL):
                                max_logL  = deepcopy(modelX_list[ID].loglkl)
                                max_model = deepcopy(modelX_list[ID])
                 
                if (i>=0):       
                    chain[:,i,j] = deepcopy(modelX_list[ID].Burst.paramsND)
                    lkl[i,j]     = deepcopy(modelX_list[ID].loglkl)
            
                if (j==0):
                    Prop_Hist.itr += 1 # add on an iteration
                    Prop_Hist.history[i%N_hist,:] = deepcopy(modelX_list[ID].Burst.paramsND)
                    
                # update the fisher matrix  
                if (i%100 == 0):
                    name = 'huh'
                    k = -1
                    while (name != 'Fisher'):
                        k += 1
                        name = proposal_list[k].name
                    modelX_list[who[j]].Burst.calc_Fish(X_flag)
                    M,NN = modelX_list[who[j]].Burst.Fisher.shape
                    U,s,Vh = slin.svd(modelX_list[who[j]].Burst.Fisher)
                    Sig = slin.diagsvd(s,M,NN)

                    for i in range(len(Sig)):
                        if (Sig[i,i] < 0.01):
                            Sig[i,i] = 0.    

                    Fish = U.dot(Sig.dot(Vh))
                    e_vals2, e_vecs2 = mct.get_Fisher_Eignen_BS(Fish)
                    mask = (e_vals2 > 1.0e-10)
                    e_vals2 = e_vals2[mask]
                    e_vecs2 = e_vecs2[:,mask]
                    
                    proposal_list[k].e_vals[j] = e_vals2
                    proposal_list[k].e_vecs[j] = e_vecs2
                    
                    #proposal_list[k] = Proposal_Fisher("Fisher", 0.8, e_vals2, e_vecs2)    
                    
        if (i>=0):
            T_chain[i,:] = np.copy(T_ladder)

    print("Temperature Ladder...... {}".format(T_ladder))
    print("acceptance rate......... {}%".format(acc/acc_cnt*100))
    print("swap acceptance rate.... {}%".format(swap/swap_cnt*100))
    print("max logL................ {}" .format(max_logL))          
    
    return lkl[::under_sample,:], chain[:,::under_sample,:], max_model, T_chain
    



def set_T_ladder(T_max, T_ladder):
    """ Use geometric spacing to fix Temperature ladder """
    
    gamma = 1.3
    T = T_ladder[0]
    while (T < T_max):
        T *= gamma
        T_ladder = np.append(T_ladder, T)
        
    return T_ladder
    
def condition_matrix(matrix):
    """ Use the SVD to condition a matrix """
    
    M,NN = matrix.shape
    U,s,Vh = slin.svd(matrix)
    Sig = slin.diagsvd(s,M,NN)

    for i in range(len(Sig)):
        if (Sig[i,i] < 0.01):
            Sig[i,i] = 0.

    matrix = U.dot(Sig.dot(Vh))    
    
    return matrix
    
def prior_logpdf(self, paramsND, T=1.0):
    """ GW Burst prior log PDF on parameters """
    
    result = 0.0
        
    lnA   = paramsND[wv.IDX_lnA]
    f0    = paramsND[wv.IDX_f0]*mHz
    t0    = paramsND[wv.IDX_t0]*Week
    tau   = paramsND[wv.IDX_tau]*Week
    phi0  = paramsND[wv.IDX_phi0]
    
    cost  = paramsND[bu.IDX_cost]
    phi   = paramsND[bu.IDX_phi]
    psi   = paramsND[bu.IDX_psi]
    ellip = paramsND[bu.IDX_ellip]
    
    # wavelet priors
    if not (np.log(1.0e-25) < lnA < np.log(1.0e-15)): # TODO: gotta develop and SNR prior, but for now enforce positivity
        result -= np.inf
    if not (0.1e-3 < f0 < 40.e-3):
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
   
def Prior_rvs(self, paramsND_X, T=1.0):
    """ Randomly draw from prior distribution """
    
    lnA   = ss.uniform.rvs(np.log(1.0e-25), np.log(1.0e-15)-np.log(1.0e-25))
    f0    = ss.uniform.rvs(0.1e-3, 40e-3-0.1e-3)
    t0    = ss.uniform.rvs(0, 1.5*Week)
    tau   = ss.uniform.rvs(0, 1.5*Week)
    phi0  = ss.uniform.rvs(0, 2*np.pi)
    
    cost  = ss.uniform.rvs(-1,2)
    phi   = ss.uniform.rvs(0, 2*np.pi)
    psi   = ss.uniform.rvs(0, np.pi)
    ellip = ss.uniform.rvs(0,1)
    
    paramsND = np.array([lnA, f0/mHz, t0/Week, tau/Week, phi0, cost, phi, psi, ellip])
    
    return paramsND   
   
        
class Proposal_Prior(Proposal):
    """ GW burst prior proposal class """
    
    def __init__(self, name, weight):
        self.name = name
        self.weight = weight
        
    logpdf = prior_logpdf
    rvs    = Prior_rvs   

    
def wrap_phase(phi0):
    
    # wrap the phase
    while (phi0 > 2*np.pi):
        phi0 -= 2*np.pi
    while (phi0 < 0):
        phi0 += 2*np.pi  
        
    return phi0
    
def target_logpdf(self, paramsND, T=1.0):

    result = 0
    
    f0  = paramsND[wv.IDX_f0]*mHz
    tau = paramsND[wv.IDX_fau]*Week
    
    result += ss.norm.logpdf(f0,  self.f0,  self.f0_sigma)
    result += ss.norm.lodpdf(tau, self.tau, self.tau_sigma)
    
    return result
    
def target_rvs(self, paramsND, T=1.0):
    
    paramsND[wv.IDX_f0]  = ss.norm.rvs(self.f0,  self.f0_sigma)/mHz
    paramsND[wv.IDX_tau] = ss.norm.rvs(self.tau, self.tau_sigma)/Week
    
    return paramsND
    
    
class Proposal_target(Proposal):
    """ Targeting proposal """
    
    def __init__(self, name, weight, f0, f0_sigma, tau, tau_sigma):
        self.name      = name
        self.weight    = weight
        self.f0        = f0
        self.f0_sigma  = f0_sigma
        self.tau       = tau
        self.tau_sigma = tau_sigma

    # methods
    logpdf = target_logpdf
    rvs    = target_rvs  
    
    
def Burst_MCMC2(TDI_data, Orbit, modelX0, seed, N, N_Burn, Flag, model_True):
    """ Maybe a cleaner version of the MCMC """
    
    X_flag = Flag.X
    under_sample = 50
    
    # Setup temperatures for MCMC
    N_temps = 1 # will be changing later
    T_ladder = np.array([1.0])
    
    if (Flag.PT):
        data_snr = np.sum(td.get_TDI_overlap(TDI_data, TDI_data, TDI_data.f_min, TDI_data.f_max, X_flag))
        T_max = data_snr/5**2
        T_ladder = set_T_ladder(T_max, T_ladder)
        N_temps = len(T_ladder)
    print("Num of Temps....... {}".format(N_temps))
    who = np.arange(0, N_temps, 1)
    
    # keep track of the ML solution
    max_logL = 1.0e-30 
    max_model = deepcopy(modelX0)
    
    np.random.seed(seed) # set Numpy's seed for reproducibility 
    
    # Storage of acceptance rates
    acc_cnt  = np.zeros(N_temps) 
    acc      = np.zeros(N_temps) 
    swap_cnt = np.zeros(N_temps-1)
    swap     = np.zeros(N_temps-1)
    
    # MCMC data storage
    Chains = np.zeros((len(modelX0.Burst.paramsND), N, N_temps))
    lkl   = np.zeros((N, N_temps))
    
    ##### set up the priors #####
    prior = Prior(burst_prior_logpdf)
    
    # create the modelX list for each temperature
    modelX_list = []
    modelX0.logprior = prior.logpdf(modelX0.Burst.paramsND)
    for i in range(N_temps):
        modelX_list.append(deepcopy(modelX0))
          
    ########################### setup the propsals ##############################

    proposal_list = []
    
    ##### Fisher #####
    modelX0.Burst.calc_Fish(X_flag)
    Fisher = condition_matrix(modelX0.Burst.Fisher)

    e_vals, e_vecs = mct.get_Fisher_Eignen_BS(Fisher)
    
    # remove the bad eigenBS
    mask = (e_vals > 1.0e-10)
    e_vals = e_vals[mask]
    e_vecs = e_vecs[:,mask]
    
    # multiple the list by number of temperatures such that each temp has its own fisher
    e_vals_ls = [e_vals]*N_temps
    e_vecs_ls = [e_vecs]*N_temps
    
    # append Fisher proposal to list of distributions
    Prop_Fish2 = Proposal_Fisher("Fisher", 0.7, e_vals_ls, e_vecs_ls, T_ladder)
    proposal_list.append(Prop_Fish2) 
    
    ##### Differential Evolution #####  
    N_hist = 10
    history = np.zeros((N_hist, len(modelX0.Burst.paramsND)))
    
    Prop_Hist = Proposal_DE("DiffEv", 0.15, history, itr=0)
    proposal_list.append(Prop_Hist)  
    
    ##### Prior #####     
    Prop_prior = Proposal_Prior("Prior", 0.05)
    proposal_list.append(Prop_prior)
    
    ##### Targeting #####
    model_True.Burst.calc_Fish(X_flag)
    Fisher = condition_matrix(model_True.Burst.Fisher)
    Cov = np.linalg.inv(Fisher)
    Prop_target = Proposal_target("target", 0.1, model_True.Burst.paramsND[wv.IDX_f0]*mHz, 
                                                 2*np.sqrt(Cov[wv.IDX_f0, wv.IDX_f0]),
                                                 model_True.Burst.paramsND[wv.IDX_tau]*Week, 
                                                 2*np.sqrt(Cov[wv.IDX_tau, wv.IDX_tau]))
    
    #############################################################################

    
    ############################### MCMC Loop ##################################    
    
    for i in tqdm(range(-N_Burn, N)):
    
        ##### Perform the normal MCMC #####
        for chain in range(N_temps):
            # Extract chain ID and temperature
            ID = int(np.copy(who[chain]))
            T  = np.copy(T_ladder[chain])

            for step in range(under_sample): # iterate over each chain some number of times (steps)
                acc_cnt[chain] += 1

                # propose a new set of parameters
                paramsND_Y, propID = Propose_Parameters(proposal_list, np.copy(modelX_list[ID].Burst.paramsND), T)
                paramsND_Y[wv.IDX_phi0] = wrap_phase(paramsND_Y[wv.IDX_phi0])
                log_priorY = prior.logpdf(paramsND_Y)
                
                if (log_priorY > -np.inf): # check that parameters lie within priors before proceeding
                
                    # generate Y burst and calculate its signal
                    Burst_Y = bu.Burst(np.copy(paramsND_Y), Orbit)
                    Burst_Y.construct_detector_tensor()
                    Burst_Y.calculate_strain()
                    Burst_Y.TDI = Burst_Y.construct_TDI(Orbit)                    
                    
                    modelY = Model(Burst_Y, Orbit) 
                    modelY.get_loglkl(TDI_data, X_flag)
                    modelY.logprior = log_priorY
                    
                    if not (np.isnan(modelY.loglkl)): # check that log-lkl is an acceptable number
                        logH  = modelY.logprior - modelX_list[ID].logprior
                        logH += proposal_list[propID].logpdf(paramsND_Y) \
                                - proposal_list[propID].logpdf(modelX_list[ID].Burst.paramsND)
                        logH += (modelY.loglkl - modelX_list[ID].loglkl)/T     
                        
                        u = np.log(ss.uniform.rvs(0,1))                              
                        if (u < logH): # accept the move
                            acc[chain] += 1
                            modelX_list[ID] = deepcopy(modelY)
                    
                            if (chain == 0 and modelX_list[ID].loglkl > max_logL):
                                #max_logL  = deepcopy(modelX_list[ID].loglkl)
                                max_model = deepcopy(modelX_list[ID])  
                                
        
        ##### Perform PTMCMC #####
        swap_cnt += 1
        for chain in np.arange(N_temps-1,0,-1): # starting at hottest chain
            whoA = np.copy(who[chain-1])
            whoB = np.copy(who[chain])
            
            beta  = (modelX_list[whoA].loglkl - modelX_list[whoB].loglkl)/T_ladder[chain]
            beta -= (modelX_list[whoA].loglkl - modelX_list[whoB].loglkl)/T_ladder[chain-1]   
            
            u = np.log(ss.uniform.rvs(0,1))
            if (u < beta):
                swap[chain-1] += 1
                who[chain-1]   = np.copy(whoB)
                who[chain]     = np.copy(whoA)         
        
        # store the chains and append to DE proposal history matrix
        if (i>=0):   
            for chain in range(N_temps):
                ID = int(np.copy(who[chain]))
    
                Chains[:,i,chain] = np.copy(modelX_list[ID].Burst.paramsND)
                lkl[i,chain]     = np.copy(modelX_list[ID].loglkl)

                if (chain==0 and i%N_hist==0):
                    Prop_Hist.itr += 1 # add on an iteration
                    Prop_Hist.history[i%N_hist,:] = np.copy(modelX_list[ID].Burst.paramsND) 
                    
        
        # Identify where in the proposal list the Fisher matrix lives
        name = 'huh'
        k = -1
        while (name != 'Fisher'):
            k += 1
            name = proposal_list[k].name  
             
        ##### update each of the Fisher matrices #####     
        for chain in range(N_temps):
            ID = int(np.copy(who[chain]))
            modelX_list[ID].Burst.calc_Fish(X_flag)
            Fisher = condition_matrix(modelX_list[ID].Burst.Fisher)
            e_vals, e_vecs = mct.get_Fisher_Eignen_BS(Fisher)
            # remove the bad eigenBS
            mask = (e_vals > 1.0e-10)
            e_vals = e_vals[mask]
            e_vecs = e_vecs[:,mask]
            
            proposal_list[k].e_vals[chain] = e_vals
            proposal_list[k].e_vecs[chain] = e_vecs

    return lkl, Chains, max_model





       