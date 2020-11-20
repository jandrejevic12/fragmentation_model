import numpy as np
import scipy as sp
import scipy.stats as sps
import math
from mpmath import hyp1f2
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.special import gamma, gammainc, gammaincc
from scipy.special import kv as bessel
import operator as op
from functools import reduce

# define scaling function distribution
class scaling_dist(sp.stats.rv_continuous):

    def _pdf(self, x, g):
        # must fit with fscale = 1, floc = 0
        z = g*(g+1)*x
        return g*(g+1)/gamma(g)/2.*pow(z,g/2.-1)*np.exp(-pow(z,0.5))

    def _cdf(self, x, g):
        # must fit with fscale = 1, floc = 0
        z = g*(g+1)*x
        return gammainc(g,pow(z,0.5))

    def _argcheck(self, g):
        shape = g > 0
        return shape

class len_dist(sp.stats.rv_continuous):

    def _pdf(self, x, g, s):
        # must fit with fscale = 1, floc = 0
        # s is the mean of the area data
        th = np.sqrt(s/(g*(g+1)))
        return 1./th/gamma(g)*pow(x/th,g-1)*np.exp(-x/th)

    def _cdf(self, x, g, s):
        # must fit with fscale = 1, floc = 0
        th = np.sqrt(s/(g*(g+1)))
        return gammainc(g, x/th)

    def _argcheck(self, g, s):
        shape = (g > 0) & (s > 0)
        return shape

def label_fun(i):
    return 10**(-i)

class rw_gamma(sp.stats.rv_continuous):
    # use with loc=0, scale=1
    def _pdf(self, x, k, g):
        th = 0.5/k/g
        return 2./np.sqrt(np.pi)/th/gamma(k*g)*pow(np.abs(x)/2/th,k*g-0.5)*bessel(0.5-k*g,np.abs(x)/th)

    def _sf(self, x, k, g):
        try:
            len(k)
            k = k[0]
            g = g[0]
        except:
            k = k
            g = g
        th = 0.5/k/g
        nu = 0.5-k*g
        #nu = k*g - 0.5
        y = np.zeros(len(x))
        for i,xx in enumerate(x):
            pow1 = pow(4,nu)
            pow2 = pow(xx/th,-nu)
            pows = pow2*pow1*pow2
            hyp1 = hyp1f2(0.5-nu,1-nu,1.5-nu,0.25*np.square(xx/th))
            hyp2 = hyp1f2(0.5,1.5,nu+1,0.25*np.square(xx/th))
            y[i] = 1+np.sqrt(np.pi)*(xx/th)/np.sin(np.pi*nu)/gamma(k*g)*(pows/(2*nu-1)/gamma(1-nu)*hyp1+1./gamma(nu+1)*hyp2)
            #y[i] = 1-1./np.sin(np.pi*nu)*(xx*np.pi/(2*th*gamma(nu+0.5)*gamma(1-nu)*gamma(1.5))*hyp1f2(0.5,1-nu,1.5,0.25*np.square(xx/th))-pow(xx/2/th,2*nu+1)*np.sqrt(np.pi)/(gamma(nu+1)*gamma(nu+1.5))*hyp1f2(nu+0.5,nu+1,nu+1.5,0.25*np.square(xx/th)))
        y[y>1] = 0
        return y

    def _argcheck(self, k, g):
        shape = k*g > 0
        return shape

# generate a random walk with gamma distributed step lengths and return the extremal positions (positive and negative)
def gamma_walk(a, th, n, seed):
    # a is the shape, th the scale, and n the number of steps
    steps = sp.stats.gamma.rvs(a, scale=th, size=n, random_state=seed)
    start = np.random.choice([-1,1])
    dirs = start*(-1)**np.arange(len(steps)) # alternating directions
    folds = steps*dirs
    pos = np.cumsum(folds)
    return pos

# generate an ensemble of gamma distributed random walks
def gamma_walk_samples(a, th, n, samples, w, seed):
    pos = np.zeros((samples,n))
    # generate all random walks
    for i in range(samples):
        pos[i,:] = gamma_walk(a, th, n, seed+i)
    # array to count how many paths escape each |w|
    counts = np.zeros(len(w))
    for j,ww in enumerate(w):
        # boolean array to indicate which samples escape
        c = np.zeros(samples).astype(np.bool)
        for k in range(n):
            # at each step, check which walks escaped
            c = c|((pos[:,k]>=ww)|(pos[:,k]<=-ww)) # if already escaped, stay True
        counts[j] += np.count_nonzero((c).astype(np.int))
    return counts


