import numpy as np
import optparse as op
import matplotlib.pyplot as plt

from pathlib import Path
from tqdm import tqdm
from scipy.stats import norm

from figaro.mixture import HDPGMM
from figaro.utils import get_priors
from figaro.plot import plot_median_cr
from figaro.load import save_density, load_density, load_data
from figaro import plot_settings

from raven.utils import make_mixtures, find_mean_weight, probability_single_star, median_reconstruction
from raven.gibbs import Gibbs

def main():

    parser = op.OptionParser()
    parser.add_option("-i", "--input", type = "string", dest = "stars_folder", help = "Folder with star files", default = None)
    parser.add_option("-b", "--bounds", type = "string", dest = "plot_bounds", help = "Plot bounds. Must be a string formatted as '[xmin, xmax]'. Quotation marks are required and scientific notation is accepted", default = None)
    parser.add_option("-o", "--output", type = "string", dest = "output", help = "Output folder. Default: same directory as samples folder", default = None)
    parser.add_option("--name", type = "string", dest = "h_name", help = "Name to be given to hierarchical inference files. Default: same name as samples folder parent directory", default = None)
    parser.add_option("-p", "--postprocess", dest = "postprocess", action = 'store_true', help = "Postprocessing", default = False)
    parser.add_option("--symbol", type = "string", dest = "symbol", help = "LaTeX-style quantity symbol, for plotting purposes", default = 'v_r')
    parser.add_option("--unit", type = "string", dest = "unit", help = "LaTeX-style quantity unit, for plotting purposes", default = '\\mathrm{km/s}')
    parser.add_option("-d", "--draws", type = "int", dest = "n_draws", help = "Number of draws for hierarchical distribution", default = 1000)
    parser.add_option("-s", "--sigma_prior", dest = "sigma_prior", type = "float", help = "Expected standard deviation (prior) for hierarchical inference.", default = 3.)
    parser.add_option("-r", "--rel_error", dest = "rel_error", type = "float", help = "Relative error for measurements without uncertainties. Default 3%", default = 0.03)
    parser.add_option("-v", "--vel_disp", dest = "vel_disp", type = "float", help = "Velocity dispersion for single stars. Default 5 km/s", default = 5.)
    parser.add_option("--n_pts", dest = "n_pts", type = "float", help = "Number of points for probability calculation", default = 10000)
    parser.add_option("--n_samples", dest = "n_samples", type = "float", help = "Number of samples for single fraction calculation", default = 10000)
    
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
    # Read hierarchical name
    if options.h_name is None:
        options.h_name = options.stars_folder.parent.parts[-1]
    # Prepare mixtures and load samples
    mixtures, bounds, names = make_mixtures(options.stars_folder, options.rel_error)
    events, _ = load_data(Path(options.output, 'events'))
    # Plot bounds
    if options.plot_bounds is not None:
        options.plot_bounds = np.array(np.atleast_2d(eval(options.plot_bounds)), dtype = np.float64)
    else:
        options.plot_bounds = bounds
    # Reconstruction
    if not options.postprocess:
        # Run hierarchical analysis
        prior_pars = get_priors(np.atleast_2d(bounds), samples = events, std = options.sigma_prior, probit = False, hierarchical = True)
        mix        = HDPGMM(np.atleast_2d(bounds), prior_pars = prior_pars, probit = False)
        draws      = np.array([mix.density_from_samples(mixtures) for _ in tqdm(range(options.n_draws), desc = options.h_name)])
        # Save draws
        save_density(draws, folder = output_draws, name = 'draws_'+options.h_name, ext = 'pkl')
    else:
        draws = load_density(Path(output_draws, 'draws_'+options.h_name+'.pkl'))
    mu, weight  = find_mean_weight(draws, options.vel_disp)
    # Plot
    single_pdf = lambda x: weight*norm(mu, options.vel_disp).pdf(x)
    plot_median_cr(draws, injected = single_pdf, bounds = options.plot_bounds[0], out_folder = options.output, name = options.h_name, label = options.symbol, unit = options.unit, hierarchical = True, injected_label = '\\sigma = {0}'.format(options.vel_disp)+'\ \\mathrm{km/s}', n_pts = int(options.n_pts))
    # Compute probability for each object
    prob_hdpgmm = np.genfromtxt(Path(options.output, 'prob_'+options.h_name+'.txt'), names = True)
    median      = median_reconstruction(prob_hdpgmm['x'], prob_hdpgmm['50'])
    probability = np.array([probability_single_star(star[0], median, mu = mu, weight = weight, bounds = options.plot_bounds, vel_disp = options.vel_disp, n_pts = int(options.n_pts)) for star in tqdm(mixtures, desc = 'p(single)')])
    # Gibbs sampling
    sampler = Gibbs(probability[:,0], probability[:,1])
    single_fraction, prob_single = sampler.run(int(options.n_samples))
    l_f, m_f, u_f = np.percentile(single_fraction, [16,50,84])
    # Save probabilities
    idx = np.argsort(prob_single)[::-1]
    with open(Path(options.output, 'p_single_{}.txt'.format(options.h_name)), 'w') as f:
        f.write('# star p_single\n')
        [f.write('{0}\t{1}\n'.format(name, p)) for name, p in zip(np.array(names)[idx], prob_single[idx])]
    # Plot single star fractions samples
    label = '$w_\\mathrm{single} ='+'{0:.2f}'.format(m_f)+'^{+'+'{0:.2f}'.format(u_f-m_f)+'}_{-'+'{0:.2f}'.format(m_f-l_f)+'}$'
    fig, ax = plt.subplots()
    ax.hist(single_fraction, histtype = 'step', density = True, label = label)
    ax.set_xlabel('$w_\\mathrm{single}$')
    ax.set_ylabel('$p(w_\\mathrm{single})$')
    leg = ax.legend(handlelength=0, handletextpad=0)
    for item in leg.legendHandles:
        item.set_visible(False)
    fig.savefig(Path(options.output, 'single_fraction_{}.pdf'.format(options.h_name)), bbox_inches = 'tight')

if __name__ == '__main__':
    main()
