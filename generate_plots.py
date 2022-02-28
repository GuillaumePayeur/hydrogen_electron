import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm
from emcee import EnsembleSampler

# Function for getting the spherical harmonic Y_lm
def spherical_harmonic(l,m):
    def Y_lm(theta,phi):
        return sph_harm(m,l,phi,theta)
    return Y_lm

# Function for getting Laguerre polynomial coefficients for quantum numbers n,l
def laguerre(n,l):
    n_r = n-l-1
    coeffs = np.zeros((n_r+1))
    coeffs[0] = 1
    for k in range(1,coeffs.shape[0]):
        coeffs[k] = coeffs[k-1]*(k-1-n_r)/(k*(k+2*l+1))
    def F(r):
        total = 0
        for order in range(coeffs.shape[0]):
            total += coeffs[order]*r**order
        total = total*(r**(l)*np.e**(-r/2))
        return total
    return F

# Function for getting the density function in real space for quantum numbers n,l,m
def log_density_function(n,l,m):
    def log_density(x):
        r,theta,psi = np.abs(x[0]),x[1],x[2]
        R = laguerre(n,l)
        Y = spherical_harmonic(l,m)
        return np.log(r**2*np.abs(np.sin(theta))*np.abs(R(r)*Y(theta,psi))**2+1e-10)
    return log_density

# Function for sampling from a wave function
def sample(log_density,num_samples):
    ndim, nwalkers = 3, 100
    p0 = np.random.randn(nwalkers, ndim)

    sampler = EnsembleSampler(nwalkers, ndim, log_density)
    sampler.run_mcmc(p0, num_samples//nwalkers)
    samples = sampler.get_chain()
    samples = np.transpose(samples,[2,0,1]).reshape((3,-1))
    return samples

# Function for plotting electron density given samples
def make_plot(samples):
    R = np.abs(samples[0])
    theta = samples[1]
    phi = samples[2]
    x = R*np.cos(phi)*np.sin(theta)
    y = R*np.sin(phi)*np.sin(theta)
    z = R*np.cos(theta)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x,y,z,s=1,c=R,alpha=0.5,cmap='hot')
    plt.show()

# Function for plotting electron density given quantum numbers
def plot(n,l,m,num_samples):
    log_density = log_density_function(n,l,m)
    r = np.random.rand(1)
    theta = 2*np.pi*np.random.rand(1)
    phi = np.pi*np.random.rand(1)
    samples = sample(log_density,num_samples)
    make_plot(samples)

if __name__ == '__main__':
    n = 3
    l = 2
    m = 0
    num_samples = 10000
    plot(n,l,m,num_samples)
