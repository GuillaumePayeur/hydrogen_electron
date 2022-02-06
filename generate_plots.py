import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm

# Function for getting the spherical harmonic Y_lm
def spherical_harmonic(l,m):
    def Y_lm(theta,phi):
        return sph_harm(m,l,theta,phi)
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
            total += r**order
        total = total*(r**(l+1)*np.e**(-r/2))
        return total
    return F

# Function to normalize a wave function
def normalize(psi):
    pass
    # to do

# Function for getting the wave function in real space for quantum numbers n,l,m
def wave_function(n,l,m):
    def psi(r,theta,psi):
        R = laguerre(n,l)
        Y = spherical_harmonic(l,m)
        return R(r)*Y(theta,psi)
    return normalize(psi)

# Function for sampling from a wave function
def sample(psi,num_samples):
    samples = np.zeros((3,num_samples))
    # to do

# Function for plotting electron density given samples
def make_plot(samples):
    pass
    # to do

# Function for plotting electron density given quantum numbers
def plot(n,l,m,num_samples):
    psi = wave_function(n,l,m)
    samples = sample(psi,num_samples))
    make_plot(samples)

if __name__ == '__main__':
    n = 1
    l = 0
    m = 0
    num_samples = 1000
    plot(n,l,m,num_samples)
