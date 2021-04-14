# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 18:40:00 2021

@author: Quirijn B. van Woerkom
Code for Phase 2 of AE4889 Special Topics in Astrodynamics
"""
###### Imports ######
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import time


###### Constants ######
m_S = 1.9891e30 # kg, mass of the Sun
m_E = 6.0477e24 # kg, mass of the Earth + Moon system
AU = 149597870.7e3 # m, astronomical unit (average distance Earth-Sun)
sigma_star = 1.53e-3 # kg/m^2, critical solar-sail loading parameter


###### Dimensionless units ######
dim_mass = m_E + m_S # kg
dim_length = AU # m
dim_time = 5022415 # s




###### Important (dimensionless) constants ######
mu = m_E/(m_E+m_S)

#%% Task 1.1a

# f_mu: the function whose root is to be found
def f_mu(x):
    return x - (1-mu)/(x + mu)**2 - mu/(x- 1 + mu)**2

# f_mu_prime: the derivative of f_mu
def f_mu_prime(x):
    return 1 + 2*(1-mu)/(x + mu)**3 + 3*mu/(x - 1 + mu)**3

def NewtonsMethod(func, funcderiv, initial_guess, maxit=1000, convdiff = 1e-16):
    """
    A function calculating a zero of func (which has known derivative funcderiv)
    with initial guess initial_guess (a float). Continues until the maximum 
    number of iterations maxit is reached or until the difference between two 
    successive approximations is less than convdiff. Returns the final estimate.
    """
    # Initialise variables
    it = 0
    x = initial_guess
    # Initialise convergence variable
    not_conv = True
    while not_conv:
        xold = x
        x = xold - func(xold)/funcderiv(xold)
        it += 1
        if it >= maxit:
            not_conv = False
        if np.abs(xold - x) < convdiff:
            not_conv = False
    return x
x_L2 = NewtonsMethod(f_mu,f_mu_prime, 1.1)

# Verification
print(f_mu(x_L2))

def BisectionRoot(func, interval, n_it):
    """
    Define a function that uses the bisection method to find the interval
    on which a root of a function can be found
    """
    # Set initial values
    xleft = interval[1]
    xright = interval[0]
    for it in range(n_it):
        midpoint = (xleft + xright)/2
        if np.sign(func(midpoint)) == np.sign(func(xleft)):
            xleft = midpoint
        else:
            xright = midpoint
        if np.abs(xleft - xright) < 1e-16:
            break
    return [xleft, xright]


x_L2_ver = BisectionRoot(f_mu, [0, 100], 10000)
