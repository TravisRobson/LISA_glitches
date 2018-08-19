import numpy as np
import scipy.linalg as slin

import Glitch as gl
import TDI as td
import scipy.stats as ss
import Wavelet as wv
from copy import deepcopy
import MCMC_tools as mct
import LISA as l

mHz = 1.0e-3
Week = 7*24*3600 

def get_log_lkl(self, TDI_data):
    """ Get the log likelihood for the glitch given the data TDI """
    
    self.Glitch.calc_snr()
    snr = self.Glitch.SNR
    
    f_min = np.max([TDI_data.f_min, self.Glitch.TDI.f_min])
    f_max = np.min([TDI_data.f_max, self.Glitch.TDI.f_max])
    
    overlap = td.get_TDI_overlap(TDI_data, self.Glitch.TDI, f_min, f_max)
    
    return np.sum(overlap) - 0.5*snr**2

def Fisher_rvs(e_vals, e_vecs, x):
    """ Given fisher matrix eigen-bs provide a random sample from its guassian approximation """
    
    e_dir = np.random.choice(range(len(e_vals)))
    jump = ss.norm.rvs(0,1)/np.sqrt(e_vals[e_dir])*e_vecs[:,e_dir]

    return x + np.append(jump, 0)
    
def Fisher_logpdf():
    """ This distribution is symmetric in location so return a constant """
    
    return 1.0
   
    
class Model:
    """ Model class """
    
    def __init__(self, Glitch, Orbit, prior, proposal):
        self.Glitch = Glitch
        self.Orbit  = Orbit 
        self.prior = prior
        self.proposal = proposal
            
    # Methods
    log_lkl = get_log_lkl
    

class Proposal:
    """ Glitch Proposal class """
    
    def __init__(self, weight, rvs_func, logpdf_func):
        self.weight = weight
        self.rvs = rvs_func
        self.logpdf = logpdf_func
 
 
def glitch_prior_logpdf(params):
    
    result = 0

    phi0 = params[wv.IDX_phi0]
    
    result += ss.uniform.logpdf(phi0, 0, 2*np.pi) # domain for phi0
    
    Sn = l.get_Sn(params[wv.IDX_f0]*mHz)
    A_lo = np.sqrt( np.sqrt(2/np.pi)*Sn/params[wv.IDX_tau]/Week )

    if (np.exp(params[wv.IDX_lnA]) < A_lo):
        result += -np.inf
        
    # Todo: Need the observation period built in here somehow...
    if (params[wv.IDX_t0]*Week < 0 or params[wv.IDX_t0]*Week > 1.5*Week):
        result += -np.inf
    if (params[wv.IDX_tau] < 0):
        result -= np.inf
    
    return result 
        
        
class Prior:
    """ Glitch Prior class """
    
    def __init__(self, logpdf_func, rvs_func=None):
        self.logpdf = logpdf_func    

    
def decide_on_proposal(self, TDI_data, T):
    """ Decide whether to take jump or not """
    
    acc = 0
    
    logH = self.modelY.prior.logpdf(self.modelY.Glitch.params) - self.modelX.prior.logpdf(self.modelX.Glitch.params)
    
    if (logH > -np.inf):
        logH += self.modelY.proposal.logpdf() - self.modelX.proposal.logpdf()
        
        logH += (self.modelY.loglkl - self.modelX.loglkl)/T
        
        logU = np.log(ss.uniform.rvs(0,1))
        
        if (logU < logH): # accept the proposal
            acc = 1
            self.modelX = deepcopy(self.modelY)

    return acc
    

class Link:
    """ Link class """
    
    def __init__(self, modelX, modelY, T):
        self.modelX = modelX
        self.modelY = modelY
        self.T = T # temperature
        
    # Methods
    decide = decide_on_proposal
    
    
class Flag:
    """ Glitch MCMC flag class """
    
    def __init__(self, comp_switch, PT):
        self.comp_switch = comp_switch
        self.PT = PT
        
    
def glitch_MCMC(TDI_data, Orbit, modelX, seed, N, Flag):
    """ Perform an MCMC in the glitch model """
    
    t = np.arange(0.0, Orbit.Tobs, Orbit.dt) # set time samples for observation
    
    max_logL = 1.0e-30
    max_model = deepcopy(modelX)
    
    np.random.seed(seed)
    
    # find max temperature and set up ladder
    N_temps = 1
    T_ladder = np.array([1.0])
    if (Flag.PT): 
        gamma = 1.3
        data_snr = np.sum(td.get_TDI_overlap(TDI_data, TDI_data, TDI_data.f_min, TDI_data.f_max))
        Tmax = data_snr/25 # so that hottest chain can see SNR 5 features
        N_temps = 5
        T_ladder = np.geomspace(1., Tmax, N_temps) # Todo: fix this..., just a choice for now...
       
    acc_cnt = 0
    acc = np.zeros(N_temps) # to store acceptance rates
    swap_cnt = np.zeros(N_temps-1)
    swap = np.zeros(N_temps-1)
    who = np.arange(0, N_temps, 1) # to keep track of which chain is at which temperature
    
    # store the accepted states and likelihoods
    chain = np.zeros((len(modelX.Glitch.params), N, N_temps))
    lkl  = np.zeros((N, N_temps))
    props = np.zeros((len(modelX.Glitch.params), N, N_temps))
    link = Link(deepcopy(modelX), deepcopy(modelX), 1.0)
        
    # create a link list at various temperatures (NOT a linked list!)
    link_ls = []
    for i in range(N_temps):
        link_ls.append(deepcopy(link))
        link_ls[i].T = T_ladder[i] # set the appropriate temperature      
    
    # calculate the fisher matrix and eigen-BS off initial states
    modelX.Glitch.calc_Fish()
     
#     M,N = modelX.Glitch.Fisher.shape
#     U,s,Vh = slin.svd(modelX.Glitch.Fisher)
#     Sig = slin.diagsvd(s,M,N)
# 
#     for i in range(len(Sig)):
#         if (Sig[i,i] < 0.1):
#             Sig[i,i] = 0.
#         
#     modelX.Glitch.Fisher = U.dot(Sig.dot(Vh))    

    e_vals, e_vecs = mct.get_Fisher_Eignen_BS(modelX.Glitch.Fisher)

    for i in range(N): 
        if (Flag.PT):
            # decide whether to swap temperatures or not...
            u = ss.uniform.rvs(0,1)
        else:
            u = 1 # i.e. force the use of normal MCMC jumps
            
        if (u<0.5):
            # propose a swap in temperatures
            k = int((N_temps-1)*ss.uniform.rvs(0,1))
            swap_cnt[k] += 1
            whoA = who[k]
            whoB = who[k+1]
            
            beta  = (link_ls[whoA].modelX.loglkl - link_ls[whoB].modelX.loglkl)/T_ladder[k+1]
            beta -= (link_ls[whoA].modelX.loglkl - link_ls[whoB].modelX.loglkl)/T_ladder[k]
            
            u = np.log(ss.uniform.rvs(0,1))
            if (u < beta):
                swap[k] += 1
                hold = deepcopy(whoA)
                whoA = deepcopy(whoB)
                whob = deepcopy(hold)
                
            for j in range(N_temps):
                ID = who[j]
                chain[:,i,j] = deepcopy(link_ls[ID].modelX.Glitch.params)
                lkl[i,j] = deepcopy(link_ls[ID].modelX.loglkl)
            
        else:
            for j in range(N_temps):
                acc_cnt += 1
                ID = who[j]
                T = T_ladder[j]
                
                # propose a new model
                paramsY = link_ls[ID].modelX.proposal.rvs(e_vals, e_vecs, \
                                                         deepcopy(link_ls[ID].modelX.Glitch.params))
                if (Flag.comp_switch): # if we are vary the component in which glitch appears... do so
                    comp_id = np.random.choice(range(len(gl.COMP_ID_LS)))
                    paramsY[-1] = comp_id
                    #paramsY[wv.IDX_t0] -= ss.norm.rvs(0, 1/(2*np.pi*l.fstar)/Week)
                props[:,i,j] = paramsY
                
                A = np.exp(paramsY[wv.IDX_lnA])
                f0 = paramsY[wv.IDX_f0]*mHz
                t0 = paramsY[wv.IDX_t0]*Week
                tau = paramsY[wv.IDX_tau]*Week
                phi0 = paramsY[wv.IDX_phi0]
                
                waveY = wv.Wavelet(A, f0, tau, t0, phi0, Orbit)
                
                if (modelX.prior.logpdf(paramsY) != -np.inf):     
                    waveY.calc_Psi()
                    waveY.make_padded_Psi(t)
            
                    glitchY = gl.Glitch(waveY, int(paramsY[-1]), Orbit) 
                    glitchY.params = paramsY
                    glitchY.calc_TDI()
                    glitchY.calc_snr()
    
                    link_ls[ID].modelY = Model(deepcopy(glitchY), Orbit, \
                                               deepcopy(modelX.prior), deepcopy(modelX.proposal))
                    link_ls[ID].modelY.loglkl = link_ls[ID].modelY.log_lkl(TDI_data)
    
                    link_ls[ID] = Link(deepcopy(link_ls[ID].modelX), link_ls[ID].modelY, 1.0)
                
                dec = link_ls[ID].decide(TDI_data, T)
                acc[j] += dec
                if (T == 1.0 and dec == 1):
                    if (link_ls[ID].modelX.loglkl > max_logL):
                        max_logL = link_ls[ID].modelX.loglkl
                        max_model = deepcopy(link_ls[ID].modelX)
    
                chain[:,i,j] = deepcopy(link_ls[ID].modelX.Glitch.params)
                lkl[i,j] = deepcopy(link_ls[ID].modelX.loglkl)
       
    print("acceptance rate......... {}%".format(acc/acc_cnt*100))
    print("swap acceptance rate.... {}%".format(swap/swap_cnt*100))
    print("max logL................ {}" .format(max_logL))
     
    return lkl, chain, props, max_model
        
    
    
    
    
    
