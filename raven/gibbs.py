import numpy as np
from scipy.stats import beta

class Gibbs:
    """
    Gibbs sampler to infer the binary fraction and the probability for each star to be single or part of a binary system.
    
    Arguments:
        np.ndarray p_single: probability for each star to be drawn from the single distribution
        np.ndarray p_binary: probability for each star to be drawn from the binary distribution
        double alpha:        symmetric Dirichlet distribution concentration parameter
    """
    def __init__(self, p_single,
                       p_binary,
                       alpha = 2.
                       ):
        self.p_single        = p_single
        self.p_binary        = p_binary
        self.alpha           = alpha
        self.samples         = []
        self.single_fraction = []
    
    def initialise(self):
        """
        Initialise the sampler
        """
        self.assignments = np.zeros(len(self.p_single))
        self.n_single    = 0
        self.n_binary    = 0
    
    def update_assignments(self, idx_star, single):
        """
        Update clusters and assignments
        
        Arguments:
            int idx_star: star to be updated
            bool single:  whether the star is assigned to the single stars cluster or not
        """
        if single:
            self.assignments[idx_star] = 1
            self.n_single             += 1
        else:
            self.assignments[idx_star] = 0
            self.n_binary             += 1
    
    def draw_assignment(self, idx_star):
        """
        Draw random assignment
        
        Arguments:
            int idx_star: star to be updated
        """
        p_s = self.p_single[idx_star]*(self.n_single + self.alpha/2.)/(self.n_single + self.n_binary + self.alpha)
        p_b = self.p_binary[idx_star]*(self.n_binary + self.alpha/2.)/(self.n_single + self.n_binary + self.alpha)
        N   = p_s + p_b
        return np.random.choice([False, True], p = [p_b/N, p_s/N])
    
    def draw_sample(self):
        """
        Draw a new single_fraction and assignments sample
        """
        self.initialise()
        order = np.arange(len(self.p_single))
        np.random.shuffle(order)
        for idx_star in order:
            self.update_assignments(idx_star, self.draw_assignment(idx_star))
        self.single_fraction.append(beta(self.n_single + self.alpha/2., self.n_binary + self.alpha/2.).rvs())
        self.samples.append(np.copy(self.assignments))
    
    def run(self, n_samples = 1e2):
        """
        Sample from the assignments space.
        
        Arguments:
            int n_samples: number of desired samples
        
        Returns:
            np.ndarray: single_fraction samples
            np.ndarray: p_single for each star
        """
        for _ in range(int(n_samples)):
            self.draw_sample()
        return self.single_fraction, np.mean(self.samples, axis = 0)
