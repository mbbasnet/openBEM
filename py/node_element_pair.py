import numpy as np
import isoparametric
from numpy.polynomial.legendre import leggauss

from openBEM.py.isoparametric import jac_normal_1D, serendip_1D

ncoord = np.array([10, 20])
elcoord = np.array(np.array([[-1, 0], [1,0], [0,3]]))
inci = 2 # quadratic element
gauss_pt, gauss_wt = leggauss(10)

def radius_param(col_pt, int_pt):

    dxr = int_pt - col_pt # differences dx
    r = np.linalg.norm(dxr) # radius as length of the vector
    dxr = dxr/r # normalize cx/dr
    return r, dxr

if(True):
    '''
    gauss integration for node element pair
    '''

    for xsi in gauss_pt:

        # Calculate element parameters
        Ni = serendip_1D(xsi, inci) # Shape function
        jac, Un = jac_normal_1D(xsi, inci, elcoord) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord) # Cartesian coordinate at the integration point

        # Calculate radius r and derivative dx/dr = [dx1/dr, dx2/dr] between collocation and integration points
        r, dxr = radius_param(ncoord, int_pt)

        # Calling fundamental solution for frequency domain elastodynamics
        
    