import numpy as np
import matplotlib.pyplot as plt

# let Latex be used, and choose a font
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def get_Fihser_Eignen_BS(Fisher):
    """ Obtain eigenvalues and eigenvectors of a Fisher Matrix """
    results  = np.linalg.eigh(Fisher)
    
    e_vals = results[0]
    e_vecs = results[1]

    return e_vals, e_vecs

def get_autocorr_length(chain):
	""" Find the autocorrelation length of an MCMC chain """
	
	term = chain - np.mean(chain)
	corr = np.correlate(term, term, mode='full')
	corr = corr[len(chain)-1:]/np.max(corr) # normalize and take the positive stride half
	len_corr = np.where(corr < 0.01)[0][0]  # find first element where corr ~ 0
	
	return len_corr
	
def plot_chain(chain, y_axis_label, n_bins, true_value=None):
	""" Plot a MCMC chain """
	
	fig, ax = plt.subplots(nrows=1,ncols=2, figsize=(8*2, 6))
	
	
	# set axes labels
	ax[0].set_ylabel(y_axis_label, fontsize=20)
	ax[0].set_xlabel('Chain Itr.', fontsize=20)
	
	# make tick look nice
	ax[0].tick_params(axis='both', which='major', labelsize=20)
	
	# plot true--or ML--value if provided
	if (true_value != None):
		ax[0].axhline(true_value, color='k', linestyle='-.', linewidth=2)
		
	ax[0].plot(chain)
	
	# now plot the histogram
	lo_bin = np.min(chain)
	hi_bin = np.max(chain)
	bins = np.linspace(lo_bin, hi_bin, n_bins)
	
	ax[1].set_xlabel(y_axis_label, fontsize=20)
	ax[1].hist(chain, bins=bins, histtype='step', density=True)
	# plot true--or ML--value if provided
	if (true_value != None):
		ax[1].axvline(true_value, color='k', linestyle='-.', linewidth=2)
	
	
	return 