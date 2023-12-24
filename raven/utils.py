import numpy as np
from pathlib import Path
from itertools import combinations
from scipy.stats import norm
from scipy.interpolate import interp1d
from figaro.utils import make_gaussian_mixture
from figaro.montecarlo import MC_integral

class gaussian_product:
    """
    Wrapper to add the .pdf method to a product of gaussians
    """
    def __init__(self, gaussians):
        self.gaussians = gaussians
    
    def pdf(self, x):
        return np.exp(self.logpdf(x))
        
    def logpdf(self, x):
        return np.sum([g[0].logpdf(x) for g in self.gaussians], axis = 0)

def make_mixtures(folder, out_folder = '.', rel_error = 0.1, error = False):
    """
    Create gaussian mixture for stars observed in globular clusters.
    If no associated uncertainty is available, a relative error is used.
    
    Arguments:
        str or Path folder:     path to stars
        str or path out_folder: path for individual star reconstruction
        double rel_error:       relative uncertainty, default 3%
        bool error:             whether to return a no-errors flag or not
    
    Returns:
        list:       mixtures
        np.ndarray: bounds
        list:       names
        bool:       error flag (optional)
    """
    parent_folder = Path(folder).parent
    means    = []
    covs     = []
    names    = []
    err      = []
    mixtures = []
    flag_glob_err = False
    for star in Path(folder).glob('*.[tc][xs][tv]'):
        star_name = star.parts[-1].split('.')[0]
        data = np.atleast_2d(np.genfromtxt(star))
        # Check for measurements with no uncertainty
        flag_no_errors = False
        if (data[:,1] == 0.).any():
            flag_no_errors = True
            idx = np.where(data[:,1] == 0)[0]
            data[:,1][idx] = np.abs(data[:,0][idx]*rel_error)
        if flag_no_errors == True:
            flag_glob_err = True
            err.append('{0}: star with no uncertainty. {1:.0f}% by hand.'.format(star_name, rel_error*100.))
        # Store means and covariances
        means.append(data[:,0])
        covs.append(np.atleast_2d(data[:,1])**2)
        names.append(star_name)
    with open(Path(parent_folder, 'no_errors.txt'), 'w') as f:
        f.write('\n'.join(err))
    # Build mixtures
    bounds   = np.array([np.min([m.min() for m in means])-3*np.sqrt(np.max([c.max() for c in covs])), np.max([m.max() for m in means])+3*np.sqrt(np.max([c.max() for c in covs]))])
    for mi, ci, name in zip(means, covs, names):
        star_folder = Path(out_folder, name)
        if not star_folder.exists():
            star_folder.mkdir(parents = True)
        # One gaussian mixture per epoch
        mixtures.append(make_gaussian_mixture(mi.T, ci.T, bounds = bounds, names = ['epoch_{}'.format(i+1) for i in range(len(mi))], out_folder = star_folder, save = True, save_samples = True, probit = False))
    if not error:
        return mixtures, bounds, names
    else:
        return mixtures, bounds, names, not(flag_glob_err)
    
def find_mean_weight(draws, vel_disp = 5.):
    """
    Find the mean (maximum of the median) and relative weight of the 5 km/s std Gaussian distribution of single stars
    
    Arguments:
        list draws:      list of figaro.mixture instances
        double vel_disp: velocity dispersion for single stars
    
    Returns:
        double: mean of the Gaussian distribution
        double: relative weight of the specific feature in the mixture
    """
    bounds = draws[0].bounds
    x      = np.linspace(bounds[0,0], bounds[0,1], 10000)
    prob   = np.median([d.pdf(x) for d in draws], axis = 0)
    prob   = prob/np.sum(prob*(x[1]-x[0]))
    mu     = x[np.argmax(prob)]
    max_p  = np.max(prob)
    return mu, max_p/norm(mu, vel_disp).pdf(mu)

def probability_single_star(epochs, rv_star, rv_dist, mu, weight, bounds, vel_disp = 5., n_pts = 1e4):
    """
    Probability for an object to be a single star.
    
    Arguments:
        mixture epochs:    mixture representing the (potentially multi-epoch) observations
        mixture rv_star:   (H)DPGMM reconstructions of the radial velocity distribution of the star
        mixture rv_dist:   (H)DPGMM reconstructions of the radial velocity distribution of the cluster
        double mu:         mean of the single star distribution. If None, is inferred from draws directly.
        double weight:     relative weight of the single star component in the overall radial velocity distribution. If None, is inferred from draws directly.
        np.ndarray bounds: bounds for integration
        double vel_disp:   velocity dispersion for single stars
        int n_pts:         number of points for the integral calculation
    
    Returns:
        double: probability for the object to be a single star
    """
    # Integrals
    bounds = np.atleast_2d(bounds)
    with np.errstate(divide = 'ignore'): # Suppress warning for underflows
        logP_single = np.log(MC_integral(gaussian_product(epochs), norm(mu, vel_disp), n_draws = int(n_pts), error = False))
        logP_binary = np.sum([np.log(MC_integral(rv_star, e, n_draws = int(n_pts), error = False)) for e in epochs])
    # Numerical stability
    norm_val = np.max([logP_binary, logP_single])
    return np.exp(logP_single - norm_val), np.exp(logP_binary - norm_val)

def variability_test_single(means, covs, threshold = 4):
    """
    Variability test introduced in Eq. 1 of Sana et al. (2012) (https://arxiv.org/pdf/1209.4638.pdf)
    """
    sigmas   = np.array([np.abs(mc[0] - mc[1])/np.sqrt(sc[0] + sc[1]) for mc, sc in zip(combinations(means, 2), combinations(covs, 2))])
    if len(sigmas) == 0:
        return None
    else:
        s_detect = np.max(sigmas)
    if s_detect > threshold:
        # Binary star
        return False
    else:
        # Single star
        return True
