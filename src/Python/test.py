import numpy as np
from numba import jit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import sys
#sys.path.append("../") 

import LISA as l
import Data as d
import Galactic_Binary as gb

f     = 0.005
fdot  = 1.1112290000e-15 

# These point the source along the -z axis
theta = np.pi
phi   = 0.0

# these are the test values I have been using for a while
# theta = np.pi/2. - (-0.048536)
# phi   = 4.576966 

A     = 5.494741e-23 
iota  = 0.345696 
psi   = 2.997043 
phi0  = 1.584432
fddot = 11./3.*fdot**2/f

T = l.YEAR

params = np.array([f*T, np.cos(theta), phi, np.log(A), np.cos(iota), psi, phi0, fdot*T**2, fddot*T**3])

orb = l.Orbit(l.YEAR)
for i in range(0,1000):
	binary = gb.GB(params, orb)