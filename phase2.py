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