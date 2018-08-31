import time
from tqdm import tqdm
from copy import deepcopy

import numpy as np
import scipy.linalg as slin
import scipy.stats as ss

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
    
    # convert weights to a distribution
    cumsum = 0
    ls = np.zeros(N_props+1)
    for i in range(N_props):
        cumsum += weights[i]
        ls[i+1] = cumsum
                
    # find which distribution was chosen
    u = ss.uniform.rvs(0,1)
    for i in range(N_props):
        if (ls[i] < u < ls[i+1]):
            break

    return proposal_list[i].rvs(paramsND_X, T), i
    
    
def set_T_ladder(T_max, T_ladder):
    """ Use geometric spacing to fix Temperature ladder """
    
    gamma = 1.2
    T = T_ladder[0]
    while (T < T_max):
        T *= gamma
        T_ladder = np.append(T_ladder, T)
        
    return T_ladder
    
    
class Proposal:
    """ GW Glitch Proposal class """
    
    def __init__(self, name, weight):
        self.name = name
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

    jump = ss.norm.rvs(0,1)/np.sqrt(e_val)*e_vecs[:,e_dir]#*np.sqrt(T)

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
    
    T_idx = np.where(T==self.T_ls)[0][0]
    
    hist = self.history[:,T_idx,:]
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

    return paramsND_X + alpha*direction#*np.sqrt(T)
    
    
class Proposal_DE(Proposal):
    """ Differential Evolution Proposal Class """
    
    def __init__(self, name,weight, history, itr, T_ls):
        self.name    = name
        self.weight  = weight
        self.history = history
        self.itr     = itr
        self.T_ls    = T_ls
        
    logpdf = DE_logpdf
    rvs    = DE_rvs


def glitch_prior_logpdf(paramsND, t0_True):
    """ GW Glitch prior log PDF on parameters """
    
    result = 0.0
    
    lnA   = paramsND[wv.IDX_lnA]
    f0    = paramsND[wv.IDX_f0]*mHz
    t0    = paramsND[wv.IDX_t0]*Week
    tau   = paramsND[wv.IDX_tau]*Week
    phi0  = paramsND[wv.IDX_phi0]
    
    # wavelet priors
    if not (np.log(1.0e-22) < lnA < np.log(1.0e-18)):
        result -= np.inf
    if not (0.1e-3 < f0 < 40.e-3):
        result -= np.inf
    if not (t0_True - 5*(3+2*np.pi)*2.5e9/l.Clight < t0 < t0_True + 5*(3+2*np.pi)*2.5e9/l.Clight):
        result -= np.inf
    if not (60 < tau < Week):
        result -= np.inf
    if not (0.0 < phi0 < 2*np.pi):
        result -= np.inf
    
    return result
    
def prior_logpdf(self, paramsND, T=1.0):
    """ GW Glitch prior log PDF on parameters """
    
    result = 0.0
        
    lnA   = paramsND[wv.IDX_lnA]
    f0    = paramsND[wv.IDX_f0]*mHz
    t0    = paramsND[wv.IDX_t0]*Week
    tau   = paramsND[wv.IDX_tau]*Week
    phi0  = paramsND[wv.IDX_phi0]
    
    # wavelet priors
    if not (np.log(1.0e-22) < lnA < np.log(1.0e-18)): 
        result -= np.inf
    if not (0.1e-3 < f0 < 40.e-3):
        result -= np.inf
    if not (self.t0_True - 5*(3+2*np.pi)*2.5e9/l.Clight < t0 < self.t0_True + 5*(3+2*np.pi)*2.5e9/l.Clight):
        result -= np.inf
    if not (60 < tau < Week):
        result -= np.inf
    if not (0.0 < phi0 < 2*np.pi):
        result -= np.inf
    
    return result
    
    
def Prior_rvs(self, paramsND_X, T=1.0):
    """ Randomly draw from prior distribution """
    
    lnA   = ss.uniform.rvs(np.log(1.0e-22), np.log(1.0e-18)-np.log(1.0e-22))
    f0    = ss.uniform.rvs(0.1e-3, 40e-3-0.1e-3)
    t0    = ss.uniform.rvs(self.t0_True - 5*(3+2*np.pi)*2.5e9/l.Clight, 2*5*(3+2*np.pi)*2.5e9/l.Clight)
    tau   = ss.uniform.rvs(60, Week-60)
    phi0  = ss.uniform.rvs(0, 2*np.pi)
    
    paramsND = np.array([lnA, f0/mHz, t0/Week, tau/Week, phi0])
    
    return paramsND   
   
        
class Proposal_Prior(Proposal):
    """ GW burst prior proposal class """
    
    def __init__(self, name, weight, t0_True):
        self.name = name
        self.weight = weight
        self.t0_True = t0_True
        
    logpdf = prior_logpdf
    rvs    = Prior_rvs   

    
def wrap_phase(phi0):
    
    # wrap the phase
    while (phi0 > 2*np.pi):
        phi0 -= 2*np.pi
    while (phi0 < 0):
        phi0 += 2*np.pi  
        
    return phi0
    

class Prior:
    """ GW glrst prior class """
    
    def __init__(self, func_logpdf, func_rvs=None):
        self.logpdf = func_logpdf  
        
        
def target_logpdf(self, paramsND, T=1.0):

    result = 0
    t1 = 0
    t2 = 0
    
    f0  = paramsND[wv.IDX_f0]
    tau = paramsND[wv.IDX_tau]
    t0  = paramsND[wv.IDX_t0] 
        
    t1 += self.norm_weight*ss.norm.pdf(f0,  self.f0,  self.f0_sigma)
    t1 += self.unif_weight*ss.uniform.pdf(f0, self.f0-self.unif_width/2*self.f0_sigma, self.unif_width*self.f0_sigma)
    
    t2 += self.norm_weight*ss.norm.pdf(tau, self.tau, self.tau_sigma)
    t2 += self.unif_weight*ss.uniform.pdf(tau, self.tau-self.unif_width/2*self.tau_sigma, self.unif_width*self.tau_sigma)
    
#     result += self.norm_weight*ss.norm.pdf(t0, self.t0, self.t0_sigma)
#     result += self.unif_weight*ss.uniform.pdf(t0, self.t0-self.unif_width/2*self.t0_sigma, self.unif_width*self.t0_sigma)
    
    result = t1*t2
    
#     result += self.norm_weight*ss.norm.pdf(t0, self.t0, self.t0_sigma)
#     result += self.unif_weight*ss.uniform.pdf(t0, self.t0-self.unif_width/2*self.t0_sigma, self.unif_width*self.t0_sigma)
    
    if (result == 0):
        return -np.inf
    else:
        return np.log(result)
    
def target_rvs(self, paramsND, T=1.0):
    
    u = ss.uniform.rvs(0,1,3)

    if (u[0] < self.unif_weight):
        paramsND[wv.IDX_f0]  = ss.uniform.rvs(self.f0-self.unif_width/2*self.f0_sigma, self.unif_width*self.f0_sigma)
    else:
        paramsND[wv.IDX_f0]  = ss.norm.rvs(self.f0,  self.f0_sigma)
        
    if (u[1] < self.unif_weight):
        paramsND[wv.IDX_tau] = ss.uniform.rvs(self.tau-self.unif_width/2*self.tau_sigma, self.unif_width*self.tau_sigma)
    else:
        paramsND[wv.IDX_tau] = ss.norm.rvs(self.tau, self.tau_sigma)
        
#     if (u[2] < self.unif_weight):
#         paramsND[wv.IDX_t0]  = ss.uniform.rvs(self.t0-self.unif_width/2*self.t0_sigma, self.unif_width*self.t0_sigma)/Week
#     else:
#         paramsND[wv.IDX_t0]  = ss.norm.rvs(self.t0, self.t0_sigma)/Week   
            
    return paramsND
    
    
class Proposal_target(Proposal):
    """ Targeting proposal """
    
    def __init__(self, name, weight, f0, f0_sigma, tau, tau_sigma, t0):
        self.name      = name
        self.weight    = weight
        self.f0        = f0
        self.f0_sigma  = f0_sigma
        self.tau       = tau
        self.tau_sigma = tau_sigma
        self.t0        = t0
        
        self.unif_weight = 0.2
        self.norm_weight = 0.8
        self.unif_width  = 10
        self.t0_sigma    = 3*2.5e9/l.Clight
        
    # methods
    logpdf = target_logpdf
    rvs    = target_rvs  
            
    
class Flag:
    """ Glitch MCMC flag class """
    
    def __init__(self, comp_switch, PT, X):
        self.comp_switch = comp_switch
        self.PT = PT
        self.X = X
        
        
def get_eigenBS(Fisher):
    """ Perform a modified Cholesky Decompisition to ensure positive-definitness of Fisher """
    
    a = np.copy(Fisher)
    N = len(Fisher)
    
    u = np.finfo(float).eps # machine precision
    gamma = 0
    for i in range(N):
        if (a[i,i] > gamma):
            gamma = a[i,i]
            
    xi = 0
    for i in range(N):
        for j in range(i+1,N):
            if (a[i,j] > xi):
                xi = np.abs(a[i,j])
    delta = u*np.max([gamma + xi, 1])
    #print("delta.... {}".format(delta))
    beta = np.sqrt( np.max([u, gamma, xi/np.sqrt(N**2 - 1)]) )
    #print("beta.... {}".format(beta))
    #print()

    c = np.zeros((N,N))
    d = np.zeros((N,N))
    l = np.zeros((N,N))

    for i in range(N):
        l[i,i] = 1

    for j in range(N):
        c[j,j] = np.copy(a[j,j])
        for s in range(j-1):
            c[j,j] -= d[s,s]*l[j,s]**2
            
        #d[j,j] = np.copy(c[j,j]) # This is the step to be modified

        if (j==0 or j==N-1):
            thetaj = 0
        else:
            thetaj = np.max(np.abs(c[j+1:,j]))

        d[j,j] = np.max([delta, np.abs(c[j,j]), (thetaj/beta)**2])


        for i in range(j+1, N):
            c[i,j] = np.copy(a[i,j])
            for s in range(j-1):
                c[i,j] -= d[s,s]*l[i,s]*l[j,s]
                
            l[i,j] = c[i,j]/d[j,j]

    result = np.dot(l, np.dot(d, l.transpose()))
    
    try:
        e_vals, e_vecs = np.linalg.eigh(result)
    except np.linalg.LinAlgError:
        e_vals = np.ones(N)*1
        e_vecs = np.zeros((N,N))
        for i in range(N):
            e_vecs[i,i] = 1
    
    return e_vals, e_vecs, result

def adapt_T_ladder(T_ladder, swap, swap_cnt, i):
    """ Adapt the temperature ladder """
    
    t0 = 10000
    nu = 10
    
    T_old = np.copy(T_ladder)
    N_temps = len(T_ladder)
    
    A = swap/swap_cnt
    mask = np.isnan(A)
    A[mask] = 0
    
    for j in range(1, N_temps-1):
        kappa = t0/(t0+i)/nu
        T_ladder[j] = T_ladder[j-1] + (T_old[j] - T_old[j-1])*np.exp(kappa*(A[j-1] - A[j]))
    
        if (T_ladder[j]/T_ladder[j-1] < 1.1):
            T_ladder[j] = 1.1*T_ladder[j-1] 
                
    return T_ladder   
        
def MCMC_glitch2(TDI_data, Orbit, modelX0, seed, N, N_Burn, Flag, model_True):
    """ A Clearner/Better glitch MCMC """
    
    X_flag = Flag.X
    under_sample = 10
    
    # Setup temperatures for MCMC
    N_temps = 1 # will be changing later
    T_ladder = np.array([1.0])
    
    if (Flag.PT):
        data_snr = np.sum(td.get_TDI_overlap(TDI_data, TDI_data, TDI_data.f_min, TDI_data.f_max, X_flag))
        T_max = data_snr/2**2
        T_ladder = set_T_ladder(T_max, T_ladder)
        N_temps = len(T_ladder)
    print("Num of Temps....... {}".format(N_temps))
    who = np.arange(0, N_temps, 1)
    T_ladder[-1] = 1.0e6

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
    Chains = np.zeros((len(modelX0.Glitch.paramsND), N, N_temps))
    lkl    = np.zeros((N, N_temps))
    comp_id_ls   = np.zeros((N, N_temps)) 
    
    ##### set up the priors #####
    prior = Prior(glitch_prior_logpdf, model_True.Glitch.paramsND[wv.IDX_t0]*Week)
    
    # create the modelX list for each temperature
    modelX_list = []
    t0_True = model_True.Glitch.paramsND[wv.IDX_t0]*Week
    modelX0.logprior = prior.logpdf(modelX0.Glitch.paramsND,  t0_True)
    for i in range(N_temps):
        modelX_list.append(deepcopy(modelX0))    
        
    ########################### setup the propsals ##############################

    proposal_list = []
    
    ##### Fisher #####
    model_True.Glitch.calc_Fish(X_flag)
    e_vals, e_vecs, Fisher = get_eigenBS(model_True.Glitch.Fisher)
    
    # multiple the list by number of temperatures such that each temp has its own fisher
    e_vals_ls = [e_vals]*N_temps
    e_vecs_ls = [e_vecs]*N_temps
    
    # append Fisher proposal to list of distributions
    Prop_Fish = Proposal_Fisher("Fisher", 0.65, e_vals_ls, e_vecs_ls, T_ladder)
    proposal_list.append(Prop_Fish) 
    
    ##### Differential Evolution #####  
    N_hist = 1
    history = np.zeros((N_hist, N_temps, len(modelX0.Glitch.paramsND)))
    
    Prop_Hist = Proposal_DE("DiffEv", 0.20, history, 0, T_ladder)
    proposal_list.append(Prop_Hist)  
    
    ##### Prior #####     
    Prop_prior = Proposal_Prior("Prior", 0.05, t0_True)
    proposal_list.append(Prop_prior)
    
    ##### Targeting #####
    model_True.Glitch.calc_Fish(X_flag)
    Cov = np.linalg.inv(model_True.Glitch.Fisher)
    Prop_target = Proposal_target("target", 0.1, model_True.Glitch.paramsND[wv.IDX_f0], 
                                                 2*np.sqrt(Cov[wv.IDX_f0, wv.IDX_f0]),
                                                 model_True.Glitch.paramsND[wv.IDX_tau], 
                                                 2*np.sqrt(Cov[wv.IDX_tau, wv.IDX_tau]),
                                                 model_True.Glitch.paramsND[wv.IDX_t0])
    proposal_list.append(Prop_target)        
 
    #############################################################################

    ############################### MCMC Loop ##################################    
    
    for i in tqdm(range(-N_Burn, N)): #range(-N_Burn, N): #tqdm(range(-N_Burn, N)):

        ##### Perform the normal MCMC #####
        for chain in range(N_temps):
            # Extract chain ID and temperature
            ID = int(np.copy(who[chain]))
            T  = np.copy(T_ladder[chain])
            
            for step in range(under_sample): # iterate over each chain some number of times (steps)
                acc_cnt[chain] += 1
                # propose a new set of parameters
                paramsND_Y, propID = Propose_Parameters(proposal_list, np.copy(modelX_list[ID].Glitch.paramsND), T)
                comp_id_Y = int(model_True.Glitch.comp_id)
                
                paramsND_Y[wv.IDX_phi0] = wrap_phase(paramsND_Y[wv.IDX_phi0])
                logpriorY = prior.logpdf(paramsND_Y,  t0_True)

                if (logpriorY > -np.inf): # check that parameters lie within priors before proceeding

                    # generate Y burst and calculate its signal
                    Glitch_Y = gl.Glitch(np.copy(paramsND_Y), Orbit, np.copy(comp_id_Y))
                    Glitch_Y.calc_TDI()                  
                    
                    modelY = Model(Glitch_Y, Orbit) 
                    modelY.get_loglkl(TDI_data, X_flag)
                    modelY.logprior = logpriorY
                    
                    if not (np.isnan(modelY.loglkl)): # check that log-lkl is an acceptable number
                        logH  = modelY.logprior - modelX_list[ID].logprior
                        logH -= proposal_list[propID].logpdf(paramsND_Y) \
                                - proposal_list[propID].logpdf(modelX_list[ID].Glitch.paramsND)
                        logH += (modelY.loglkl - modelX_list[ID].loglkl)/T  
                        
                        u = np.log(ss.uniform.rvs(0,1))                              
                        if (u < logH): # accept the move
                            acc[chain] += 1
                            modelX_list[ID] = deepcopy(modelY)
                    
                            if (chain == 0 and modelX_list[ID].loglkl > max_logL):
                                max_logL  = np.copy(modelX_list[ID].loglkl)
                                max_model = deepcopy(modelX_list[ID])  
                                
        
        ##### Perform PTMCMC #####
        swap_cnt += 1
        for chain in np.arange(N_temps-1,0,-1): # starting at hottest chain
            whoA = int(np.copy(who[chain-1]))
            whoB = int(np.copy(who[chain]))
            
            beta  = (modelX_list[whoA].loglkl - modelX_list[whoB].loglkl)/T_ladder[chain]
            beta -= (modelX_list[whoA].loglkl - modelX_list[whoB].loglkl)/T_ladder[chain-1]   
            
            u = np.log(ss.uniform.rvs(0,1))
            if (u < beta):
                swap[chain-1] += 1
                who[chain-1]   = int(np.copy(whoB))
                who[chain]     = int(np.copy(whoA))
                
        T_ladder = adapt_T_ladder(np.copy(T_ladder), swap, swap_cnt, i+N_Burn)
        
        ####### store the chains  ######
        if (i>=0):   
            for chain in range(N_temps):
                ID = int(np.copy(who[chain]))
    
                Chains[:,i,chain] = np.copy(modelX_list[ID].Glitch.paramsND)
                lkl[i,chain]      = np.copy(modelX_list[ID].loglkl)
                comp_id_ls[i,chain] = int(np.copy(modelX_list[ID].Glitch.comp_id))
                
        ##### Append to DE proposal history matrix ######
        if (i%N_hist==0):
            Prop_Hist.itr += 1 # add on an iteration
            for chain in range(N_temps):
                ID = int(np.copy(who[chain]))
                Prop_Hist.history[i%N_hist,chain,:] = np.copy(modelX_list[ID].Glitch.paramsND) 
                    
        
        # Identify where in the proposal list the Fisher matrix lives
        name = 'huh'
        k = -1
        while (name != 'Fisher'):
            k += 1
            name = proposal_list[k].name  
        proposal_list[k].T_ls = np.copy(T_ladder) # need to update Temperature ladder
        ##### update each of the Fisher matrices #####     
        if (i%20 == 0):
            for chain in range(N_temps):
                ID = int(np.copy(who[chain]))
                modelX_list[ID].Glitch.calc_Fish(X_flag)
                e_vals, e_vecs, Fisher = get_eigenBS(modelX_list[ID].Glitch.Fisher)
            
                proposal_list[k].e_vals[chain] = np.copy(e_vals)
                proposal_list[k].e_vecs[chain] = np.copy(e_vecs)
                
        name = 'huh'
        k = -1
        while (name != 'DiffEv'):  # need to update Temperature ladder
            k += 1
            name = proposal_list[k].name  
        proposal_list[k].T_ls = np.copy(T_ladder) 
    
    ##### Print some summary results #####
    print("Temperature Ladder...... {}".format(T_ladder))
    print("acceptance rate......... {}%".format(acc/acc_cnt*100))
    print("swap acceptance rate.... {}%".format(swap/swap_cnt*100))
    print("max logL................ {}" .format(max_logL))        

    return lkl, Chains, comp_id_ls, max_model



    

        
