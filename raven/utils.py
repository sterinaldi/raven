import numpy as np
from pathlib import Path
from scipy.stats import norm
from scipy.interpolate import interp1d
from figaro.utils import make_gaussian_mixture
from figaro.montecarlo import MC_integral

class median_reconstruction:
    """
    Wrapper to add the .pdf method to SciPy's interp1d
    """
    def __init__(self, x, p_x):
        self.interpolant = interp1d(x, p_x, bounds_error = False, fill_value = 0.)
    
    def pdf(self, x):
        return self.interpolant(x)

def make_mixtures(folder, rel_error = 0.03):
    """
    Create gaussian mixture for stars observed in globular clusters.
    If no associated uncertainty is available, a relative error is used.
    
    Arguments:
        str or Path Folder: path to stars
        double rel_error:   relative uncertainty, default 3%
    
    Returns:
        list:       mixtures
        np.ndarray: bounds
        list:       names
    """
    parent_folder = Path(folder).parent
    means = []
    covs  = []
    names = []
    for star in Path(folder).glob('*.[tc][xs][tv]'):
        star_name = star.parts[-1].split('.')[0]
        data = np.atleast_2d(np.genfromtxt(star))
        # Check for measurements with no uncertainty
        flag_no_errors = False
        if (data[:,1] == 0.).any():
            flag_no_errors = True
            idx = np.where(data[:,1] == 0)[0]
            data[:,1][idx] = data[:,0][idx]*rel_error
        if flag_no_errors == True:
            with open(Path(parent_folder, 'no_errors.txt'), 'a') as f:
                f.write('{0}: star with no uncertainty. {1:.0f}% by hand.\n'.format(star_name, rel_error*100.))
        # Store means and covariances
        means.append(data[:,0])
        covs.append(np.atleast_2d(data[:,1])**2)
        names.append(star_name)
    # Build mixtures
    bounds   = np.array([np.min([m.min() for m in means])-3*np.sqrt(np.max([c.max() for c in covs])), np.max([m.max() for m in means])+3*np.sqrt(np.max([c.max() for c in covs]))])
    mixtures = make_gaussian_mixture(means, covs, bounds = bounds, names = names, out_folder = parent_folder, save = True, save_samples = True, probit = False)
    return mixtures, bounds, names

def find_mean_weight(draws, vel_disp = 5.):
    """
    Find the mean (maximum of the median) and relative weight of the 5 km/s std Gaussian distribution of single stars
    
    Arguments:
        list draws: list of figaro.mixture instances
        double vel_disp: velocity dispersion for single stars
    
    Returns:
        double: mean of the Gaussian distribution
        double: relative weight of the specific feature in the mixture
    """
    bounds = draws[0].bounds
    x      = np.linspace(bounds[0,0], bounds[0,1], 1000)
    prob   = np.median([d.pdf(x) for d in draws], axis = 0)
    prob   = prob/np.sum(prob*(x[1]-x[0]))
    mu     = x[np.argmax(prob)]
    max_p  = np.max(prob)
    return mu, max_p/norm(mu, vel_disp).pdf(mu)

def probability_single_star(star, median, mu = None, weight = None, vel_disp = 5., n_draws = 1e6):
    """
    Probability for an object to be a single star.
    
    Arguments:
        mixture star: mixture representing the (potentially multi-epoch) observations
        callable median: (H)DPGMM median reconstructions of the radial velocity distribution
        double mu: mean of the single star distribution. If None, is inferred from draws directly.
        double weight: relative weight of the single star component in the overall radial velocity distribution. If None, is inferred from draws directly.
        double vel_disp: velocity dispersion for single stars
        double n_draws: number of MC draws
    
    Returns:
        double: probability for the object to be a single star
    """
    # TODO: rewrite with analytical integral (now p slightly > 1 from time to time)
    if mu is None or weight is None:
        mu, weight = find_mean_weight(draws, vel_disp)
    # Integrals
    single_dist = norm(mu, vel_disp)
    non_par     = MC_integral(median, star, n_draws = n_draws, error = False)
    p_single    = MC_integral(single_dist, star, n_draws = n_draws, error = False)
    p_binary    = (non_par - weight*p_single)/(1. + weight)
    return p_single/(p_single + p_binary)
