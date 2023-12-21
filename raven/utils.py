import numpy as np
from pathlib import Path
from itertools import combinations
from scipy.stats import norm
from scipy.interpolate import interp1d
from figaro.utils import make_gaussian_mixture

class median_reconstruction:
    """
    Wrapper to add the .pdf method to SciPy's interp1d
    """
    def __init__(self, x, p_x):
        self.interpolant = interp1d(x, p_x, bounds_error = False, fill_value = 'extrapolate', kind = 'cubic')
    
    def pdf(self, x):
        return self.interpolant(x)

def make_mixtures(folder, rel_error = 0.03, error = False):
    """
    Create gaussian mixture for stars observed in globular clusters.
    If no associated uncertainty is available, a relative error is used.
    
    Arguments:
        str or Path Folder: path to stars
        double rel_error:   relative uncertainty, default 3%
        bool error:         whether to return a no-errors flag or not
    
    Returns:
        list:       mixtures
        np.ndarray: bounds
        list:       names
        bool:       error flag (optional)
    """
    parent_folder = Path(folder).parent
    means = []
    covs  = []
    names = []
    err   = []
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
    mixtures = make_gaussian_mixture(means, covs, bounds = bounds, names = names, out_folder = parent_folder, save = True, save_samples = True, probit = False)
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

def probability_single_star(star, median, mu, weight, bounds, vel_disp = 5., n_pts = 1e5):
    """
    Probability for an object to be a single star.
    
    Arguments:
        mixture star:      mixture representing the (potentially multi-epoch) observations
        callable median:   (H)DPGMM median reconstructions of the radial velocity distribution
        double mu:         mean of the single star distribution. If None, is inferred from draws directly.
        double weight:     relative weight of the single star component in the overall radial velocity distribution. If None, is inferred from draws directly.
        np.ndarray bounds: bounds for integration
        double vel_disp:   velocity dispersion for single stars
        int n_pts:         number of points for the integral calculation
    
    Returns:
        double: probability for the object to be a single star
    """
    # Integrals
    bounds      = np.atleast_2d(bounds)
    single_dist = norm(mu, vel_disp)
    x           = np.linspace(bounds[0,0], bounds[0,1], int(n_pts))
    dx          = x[1]-x[0]
    non_par     = np.sum(median.pdf(x)*star.pdf(x)*dx)
    p_single    = np.sum(single_dist.pdf(x)*star.pdf(x)*dx)
    p_binary    = (non_par - weight*p_single)/(1. - weight)
    return p_single, np.max([p_binary, 0]) # accounts for potential numerical error when p_single ~ p_nonpar (w ~ 1)

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
