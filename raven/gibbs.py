import numpy as np
from collections import Counter
from scipy.stats import beta

class Gibbs:
    """
    Gibbs sampler to infer the binary fraction and the probability for each star to be single or part of a binary system.
    
    Arguments:
        np.ndarray p_single: probability for each star to be drawn from the single distribution
        np.ndarray p_binary: probability for each star to be drawn from the binary distribution
        double alpha:        symmetric Dirichlet distribution concentration parameter
    """
    def __init__(self, p_pop,
                       p_single,
                       alpha = 2.
                       ):
        self.p_pop           = p_pop
        self.p_single        = p_single
        self.alpha           = alpha
        self.n_pop           = len(self.p_single[0])
        self.samples_pop     = []
        self.samples_single  = []
        self.single_fraction = []
    
    def initialise(self):
        """
        Initialise the sampler
        """
        self.assignments_pop    = np.zeros(len(self.p_single))
        self.assignments_single = np.zeros(len(self.p_single))
        self.n_single           = np.zeros(self.n_pop)
        self.n_binary           = np.zeros(self.n_pop)
        self.n_stars            = np.zeros(self.n_pop)
    
    def update_assignments_pop(self, idx_star, idx_pop):
        """
        Update clusters and assignments
        
        Arguments:
            int idx_star: star to be updated
            bool single:  whether the star is assigned to the single stars cluster or not
        """
        self.assignments_pop[idx_star] = idx_pop
        self.n_stars[idx_pop]         += 1
    
    def update_assignments_single(self, idx_star, single, idx_pop):
        """
        Update clusters and assignments
        
        Arguments:
            int idx_star: star to be updated
            bool single:  whether the star is assigned to the single stars cluster or not
        """
        if single:
            self.assignments_single[idx_star] = 1
            self.n_single[idx_pop]           += 1
        else:
            self.assignments_single[idx_star] = 0
            self.n_binary[idx_pop]           += 1
    
    def draw_assignment_pop(self, idx_star):
        """
        Draw random assignment for population
        
        Arguments:
            int idx_star: star to be updated
        
        Returns:
            int: population to assign the star to
        """
        p = self.p_pop[idx_star]*(self.n_stars + 1)/(np.sum(self.n_stars) + self.n_pop)
        return np.random.choice(self.n_pop, p = p/np.sum(p))
    
    def draw_assignment_single(self, idx_star, idx_pop):
        """
        Draw random assignment
        
        Arguments:
            int idx_star: star to be updated
        """
        p_s = self.p_single[idx_star,idx_pop,0]*(self.n_single[idx_pop] + self.alpha/2.)/(self.n_single[idx_pop] + self.n_binary[idx_pop] + self.alpha)
        p_b = self.p_single[idx_star,idx_pop,1]*(self.n_binary[idx_pop] + self.alpha/2.)/(self.n_single[idx_pop] + self.n_binary[idx_pop] + self.alpha)
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
            idx_pop = self.draw_assignment_pop(idx_star)
            self.update_assignments_pop(idx_star, idx_pop)
            ss = self.draw_assignment_single(idx_star, idx_pop)
            self.update_assignments_single(idx_star, ss, idx_pop)
        self.single_fraction.append(beta(self.n_single + self.alpha/2., self.n_binary + self.alpha/2.).rvs())
        self.samples_single.append(np.copy(self.assignments_single))
        self.samples_pop.append(np.copy(self.assignments_pop))
    
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
        count = [Counter(s) for s in np.array(self.samples_pop).T]
        draw_p_pop = np.array([np.array([c[idx_pop] for idx_pop in range(self.n_pop)])/n_samples for c in count])
        self.single_fraction = np.atleast_2d(self.single_fraction)
        if self.n_pop == 1:
            self.single_fraction = self.single_fraction.T
        return np.mean(self.single_fraction, axis = 0), np.mean(self.samples_single, axis = 0), draw_p_pop
