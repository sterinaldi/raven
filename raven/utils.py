import numpy as np
from pathlib import Path
from itertools import combinations
from scipy.stats import norm
from scipy.interpolate import interp1d
from scipy.spatial.distance import jensenshannon as js
from numba import njit
from figaro.utils import make_gaussian_mixture
from figaro.montecarlo import MC_integral

@njit
def norm_pdf(x, mu, sigma):
    return np.exp(-0.5*((x-mu)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)

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
        mixtures.append(make_gaussian_mixture(mi.T, ci.T, bounds = bounds, names = ['epoch_{}'.format(i+1) for i in range(len(mi))], out_folder = star_folder, save = True, save_samples = True, probit = False, make_comp = False))
    if not error:
        return mixtures, bounds, names
    else:
        return mixtures, bounds, names, not(flag_glob_err)
    
def find_mean_weight(draws, n_populations = 1, vel_disp = None):
    """
    Find the mean (maximum of the median) and relative weight of the 5 km/s std Gaussian distribution of single stars
    
    Arguments:
        list draws:      list of figaro.mixture instances
        double vel_disp: velocity dispersion for single stars
    
    Returns:
        double: mean of the Gaussian distribution
        double: relative weight of the specific feature in the mixture
    """
    bounds          = draws[0].bounds
    x               = np.linspace(bounds[0,0], bounds[0,1], 10000)
    prob            = np.median([d.pdf(x) for d in draws], axis = 0)
    prob            = prob/np.sum(prob*(x[1]-x[0]))
    interp_med      = interp1d(x, prob)
    interp_med_full = interp1d(x, prob)
    mu = []
    w  = []
    infer_sigma = False
    if vel_disp is None:
        vel_disp    = []
        sigma       = np.linspace(0.5,10,200)
        infer_sigma = True
    elif hasattr(vel_disp, '__iter__'):
        if not len(vel_disp) == n_populations:
            raise Exception('Please provide eiter one or n_population values for vel_disp')
    else:
        vel_disp = np.ones(n_populations)*vel_disp
    for j in range(int(n_populations)):
        mu.append(x[np.argmax(prob)])
        max_p  = interp_med_full(mu[-1])
        if infer_sigma:
            dist       = np.zeros(len(sigma))
            for i, s in enumerate(sigma):
                xi = np.linspace(mu[j]-s, mu[j]+s,100)
                dist[i] = js(norm(mu[j],s).pdf(xi), interp_med(xi))
            vel_disp.append(sigma[np.argmin(dist)])
        w.append(max_p/norm(mu[j], vel_disp[j]).pdf(mu[j]))
        prob      -= w[j]*norm(mu[j], vel_disp[j]).pdf(x)
        prob       = prob/np.sum(prob*(x[1]-x[0]))
        interp_med = interp1d(x, prob)
    return mu, w, vel_disp

def probability_population(rv_cm, pop_pars, bounds, idx):
    if idx == -1:
        return 1./np.diff(bounds)[0]
    else:
        mu_pop = pop_pars[idx][0]
        s_pop  = pop_pars[idx][1]
        return np.mean([np.sum(d.w*norm_pdf(d.means.flatten(), mu_pop, np.sqrt(d.covs.flatten() + s_pop**2))) for d in rv_cm])
    

def probability_single_star(epochs, rv_star, rv_dist, population, bounds, n_pts = 1e4):
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
        logP_single = np.log(MC_integral(gaussian_product(epochs), population, n_draws = int(n_pts), error = False))
        logP_binary = np.sum([np.log(MC_integral(rv_star, e, n_draws = int(n_pts), error = False)) for e in epochs])
    # Numerical stability
    norm_val = np.max([logP_binary, logP_single])
    return np.exp(logP_single - norm_val), np.exp(logP_binary - norm_val)

def variability_test_single(means, covs, threshold = 4, min_variation = 20):
    """
    Variability test introduced in Eq. 1 of Sana et al. (2013) (https://arxiv.org/pdf/1209.4638.pdf)
    """
    sigmas   = np.array([np.abs(mc[0] - mc[1])/np.sqrt(sc[0] + sc[1]) for mc, sc in zip(combinations(means, 2), combinations(covs, 2))])
    deltas   = np.array([np.abs(mc[0] - mc[1]) for mc in combinations(means, 2)])
    if len(sigmas) == 0:
        return None
    else:
        s_detect = np.max(sigmas)
    if s_detect > threshold:
        # Binary star
        if np.max(deltas) > min_variation:
            return (False, True)
        else:
            return (False, False)
    else:
        # Single star
        return (True, None)
