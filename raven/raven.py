import numpy as np
import optparse as op
import matplotlib.pyplot as plt

from pathlib import Path
from tqdm import tqdm
from scipy.stats import norm

from figaro.mixture import HDPGMM
from figaro.utils import get_priors, rvs_median
from figaro.load import save_density, load_density, load_data
from figaro.plot import plot_median_cr as figaro_plot_median_cr

from raven.utils import make_mixtures, find_mean_weight, probability_single_star, variability_test_single
from raven.gibbs import Gibbs
from raven.plot import plot_single_fraction, plot_median_cr, plot_p_single

def main():

    parser = op.OptionParser()
    parser.add_option("-i", "--input", type = "string", dest = "stars_folder", help = "Folder with star files", default = None)
    parser.add_option("-b", "--bounds", type = "string", dest = "plot_bounds", help = "Plot bounds. Must be a string formatted as '[xmin, xmax]'. Quotation marks are required and scientific notation is accepted", default = None)
    parser.add_option("-o", "--output", type = "string", dest = "output", help = "Output folder. Default: same directory as samples folder", default = None)
    parser.add_option("--name", type = "string", dest = "h_name", help = "Name to be given to hierarchical inference files. Default: same name as samples folder parent directory", default = None)
    parser.add_option("-p", "--postprocess", dest = "postprocess", action = 'store_true', help = "Postprocessing", default = False)
    parser.add_option("--skip_stars", dest = "skip_stars", action = 'store_true', help = "Skip reconstruction of individual stars", default = False)
    parser.add_option("--symbol", type = "string", dest = "symbol", help = "LaTeX-style quantity symbol, for plotting purposes", default = 'v_{cm}')
    parser.add_option("--unit", type = "string", dest = "unit", help = "LaTeX-style quantity unit, for plotting purposes", default = '\\mathrm{km/s}')
    parser.add_option("-d", "--draws", type = "int", dest = "n_draws", help = "Number of draws for hierarchical distribution", default = 10000)
    parser.add_option("--star_draws", type = "int", dest = "n_star_draws", help = "Number of draws for individual stars", default = 1000)
    parser.add_option("-s", "--sigma_prior", dest = "sigma_prior", type = "float", help = "Expected standard deviation (prior) for hierarchical inference.", default = 2.)
    parser.add_option("--sigma_prior_star", dest = "sigma_prior_star", type = "float", help = "Expected standard deviation (prior) for object RV inference.", default = 1.)
    parser.add_option("-r", "--rel_error", dest = "rel_error", type = "float", help = "Relative error for measurements without uncertainties. Default 10%", default = 0.1)
    parser.add_option("-v", "--vel_disp", dest = "vel_disp", type = "float", help = "Velocity dispersion for single stars. Default 5 km/s", default = 5.)
    parser.add_option("--n_pts", dest = "n_pts", type = "float", help = "Number of points for probability calculation", default = 10000)
    parser.add_option("--n_samples", dest = "n_samples", type = "float", help = "Number of samples for single fraction calculation", default = 1000)
    parser.add_option("--sana_variability", dest = "sana_variability", action = 'store_true', help = "Use variability test from Sana et al (2012)", default = False)
    
    (options, args) = parser.parse_args()

    # Paths
    if options.stars_folder is not None:
        options.stars_folder = Path(options.stars_folder).resolve()
    else:
        raise Exception("Please provide path to samples.")
    options.output = options.stars_folder.parent
    output_draws = Path(options.output, 'draws')
    if not output_draws.exists():
        output_draws.mkdir()
    output_individual = Path(options.output, 'individual_objects')
    if not output_individual.exists():
        output_individual.mkdir()
    # Read hierarchical name
    if options.h_name is None:
        options.h_name = options.stars_folder.parent.parts[-1]
    # Prepare mixtures and load samples
    star_mixtures, bounds, names, provided_errors = make_mixtures(options.stars_folder, output_individual, options.rel_error, error = True)
    if options.sana_variability and provided_errors:
        options.single = {name: variability_test_single([ss[0].means for ss in mix], [ss[0].covs for ss in mix]) for mix, name in zip(star_mixtures, names)}
    else:
        options.single = None
    # Plot bounds
    if options.plot_bounds is not None:
        options.plot_bounds = np.array(np.atleast_2d(eval(options.plot_bounds)), dtype = np.float64)[0]
    else:
        options.plot_bounds = bounds
    # Reconstruction
    if not options.postprocess:
        if not options.skip_stars:
            individual_star_draws_fine   = []
            individual_star_draws_coarse = []
            for star in tqdm(names, desc = 'NP RV objects'):
                # Load data
                star_folder = Path(output_individual, star)
                fine_star_folder = Path(star_folder, 'fine')
                if not fine_star_folder.exists():
                    fine_star_folder.mkdir()
                coarse_star_folder = Path(star_folder, 'broad')
                if not coarse_star_folder.exists():
                    coarse_star_folder.mkdir()
                epochs, _   = load_data(Path(star_folder, 'events'), verbose = False)
                mixtures    = load_density(Path(star_folder, 'draws', 'posteriors_single_event.json'), make_comp = False)
                # Issue with single-epoch stars
                if len(np.shape(mixtures)) == 1:
                    mixtures = [mixtures]
                # Individual star analysis
                # Fine structure
                prior_pars  = get_priors(np.atleast_2d(bounds), std = options.sigma_prior_star, probit = False, hierarchical = True)
                mix         = HDPGMM(np.atleast_2d(bounds), prior_pars = prior_pars, probit = False)
                draws_fine  = [mix.density_from_samples(mixtures) for _ in range(options.n_star_draws)]
                # Save density
                save_density(draws_fine, folder = fine_star_folder, name = 'draws_'+star, ext = 'json')
                individual_star_draws_fine.append(draws_fine)
                # Coarse structure
                prior_pars   = get_priors(np.atleast_2d(bounds), probit = False, hierarchical = True)
                mix          = HDPGMM(np.atleast_2d(bounds), prior_pars = prior_pars, probit = False)
                draws_coarse = [mix.density_from_samples(mixtures) for _ in range(options.n_star_draws)]
                # Save density
                save_density(draws_coarse, folder = coarse_star_folder, name = 'draws_'+star, ext = 'json')
                individual_star_draws_coarse.append(draws_coarse)
                # Plot
                figaro_plot_median_cr(draws_fine, bounds = bounds, out_folder = star_folder, name = star, label = options.symbol, unit = options.unit)
            individual_star_draws_fine   = np.array(individual_star_draws_fine)
            individual_star_draws_coarse = np.array(individual_star_draws_coarse)
            save_density(individual_star_draws_fine, folder = output_individual, name = 'draws_individual_objects_fine', ext = 'json')
            save_density(individual_star_draws_coarse, folder = output_individual, name = 'draws_individual_objects_coarse', ext = 'json')
        else:
            individual_star_draws_fine   = load_density(Path(output_individual, 'draws_individual_objects_fine.json'), make_comp = False)
            individual_star_draws_coarse = load_density(Path(output_individual, 'draws_individual_objects_coarse.json'), make_comp = False)
        # Draws from median RV posteriors
        events = np.array([rvs_median(d, 1000) for d in individual_star_draws_fine])
        # Run hierarchical analysis
        individual_star_draws_copy = np.copy(individual_star_draws_fine) # Avoid shuffling issues
        prior_pars = get_priors(np.atleast_2d(bounds), samples = events, std = options.sigma_prior, probit = False, hierarchical = True)
        mix        = HDPGMM(np.atleast_2d(bounds), prior_pars = prior_pars, probit = False)
        draws      = np.array([mix.density_from_samples(individual_star_draws_copy) for _ in tqdm(range(options.n_draws), desc = 'NP RV cluster')])
        # Save draws
        save_density(draws, folder = output_draws, name = 'draws_'+options.h_name, ext = 'json')
    else:
        individual_star_draws_coarse = load_density(Path(output_individual, 'draws_individual_objects_coarse.json'), make_comp = False)
        draws                        = load_density(Path(output_draws, 'draws_'+options.h_name+'.json'), make_comp = False)
    mu, weight = find_mean_weight(draws, options.vel_disp)
    # Plot
    plot_median_cr(draws, mu = mu, vel_disp = options.vel_disp, weight = weight, bounds = options.plot_bounds, out_folder = options.output, name = options.h_name, label = options.symbol, unit = options.unit)
    
    # Compute probability for each object (Gibbs sampling)
    single_fraction = []
    prob_single     = []
    for _ in tqdm(range(int(options.n_samples)), desc = 'p(single)'):
        probability = np.array([probability_single_star(epochs, rv_star, draws, mu = mu, weight = weight, bounds = options.plot_bounds, vel_disp = options.vel_disp, n_pts = int(options.n_pts)) for epochs, rv_star in zip(star_mixtures, individual_star_draws_coarse)])
        # Gibbs sampling
        sampler = Gibbs(probability[:,0], probability[:,1])
        sf, p_s = sampler.run()
        single_fraction.append(sf)
        prob_single.append(p_s)
    prob_single = np.array(prob_single)
    single_fraction = np.random.choice(np.concatenate([f for f in single_fraction]), size = int(10*options.n_samples))
    np.savetxt(Path(options.output, 'samples_fraction_{}.txt'.format(options.h_name)), single_fraction)
    # Save probabilities
    ll, mm, uu = np.percentile(prob_single, [16,50,84], axis = 0)
    idx = np.argsort(mm)
    with open(Path(options.output, 'p_single_{}.txt'.format(options.h_name)), 'w') as f:
        f.write('# star p_single\n')
        [f.write('{0}\t{1:.2f}+{2:.2f}-{3:.2f}\n'.format(name, p, u-p, p-l)) for name, p, u, l in zip(np.array(names)[idx[::-1]], mm[idx[::-1]], uu[idx[::-1]], ll[idx[::-1]])]
    # Plot single star fractions samples
    plot_single_fraction(single_fraction, out_folder = options.output, name = options.h_name)
    plot_p_single(probs = prob_single.T[idx], stars = np.array(names)[idx], out_folder = options.output, name = options.h_name, single = options.single)

if __name__ == '__main__':
    main()
