import numpy as np
import isoparametric
from numpy.polynomial.legendre import leggauss

from openBEM.py.isoparametric import jac_normal_1D, serendip_1D
from openBEM.py.fundamental_solution import FundamentalSolution

ncoord = np.array([10, 20])
elcoord = np.array(np.array([[-1, 0], [1,0], [0,3]]))
inci = 2 # quadratic element
gauss_pt, gauss_wt = leggauss(10)
funda = FundamentalSolution(c1, c2, mu, nu, omega) # fundamental solution for elastic isotropic case

def radius_param(col_pt, int_pt):

    dxr = int_pt - col_pt # differences dx
    r = np.linalg.norm(dxr) # radius as length of the vector
    dxr = dxr/r # normalize cx/dr
    return r, dxr
        
        
def bem_elm_eldyn():
    '''
    Computes element BEM matrices for each element node pair
    for normal case i.e. collocation node not in the element
    '''
    # Reset/initialize the local G and H matrices for new element
    G = 0.0
    Hs = 0.0
    H = 0.0
    
    
    # For normal element node pair (node outside element)
    for xsi in gauss_pt:

        # Calculate element parameters
        Ni = serendip_1D(xsi, inci) # Shape function
        jac, Un = jac_normal_1D(xsi, inci, elcoord) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord) # Cartesian coordinate at the integration point

        # Calculate radius r and derivative dx/dr = [dx1/dr, dx2/dr] between collocation and integration points
        r, dxr = radius_param(ncoord, int_pt)
        drn = np.dot(dxr, Un) # dr/dn
        
        # Calculate reciprocal displacement and tractions for normal case (no singularity)
        funda.calculateDisplacement(dxr)
        funda.calculateTraction_dyn(r, dxr, Un, drn)
        funda.calculateTraction_st(r, dxr, Un, drn)
        
        # create local variables of reciprocal displacement and tractions
        Uf = funda.displacement
        Tdf = funda.Traction_dyn
        Ts = funda.Traction_st
        
        # Calculate G and H matrices
        elm_par += jac*gauss_wt*Ni
        G += np.dot(Uf,elm_par)
        Hs += np.dot(Tdf,elm_par)
        H += np.dot(Ts,elm_par)
    
    return G, Hs, H
        
def bem_elm_eldyn_st():
    '''
    Computes element BEM matrices for each element node pair
    for collocation node start node of the quadratic element
    '''
    # Reset/initialize the local G and H matrices for new element
    G = 0.0
    Hs = 0.0
    H = 0.0
    
    # singular part
    for eta in gauss_pt_log: 
        xsi = 2*eta-1 #eta = (1+xsi)/2 => xsi = 2*eta-1 # for start node

        # Calculate element parameters
        Ni = serendip_1D(xsi, inci) # Shape function
        jac, Un = jac_normal_1D(xsi, inci, elcoord) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord) # Cartesian coordinate at the integration point

        # Calculate radius r and derivative dx/dr = [dx1/dr, dx2/dr] between collocation and integration points
        r, dxr = radius_param(ncoord, int_pt)
        drn = np.dot(dxr, Un) # dr/dn
        
        # Calculate reciprocal displacement and tractions for normal case (no singularity)
        funda.calculateDisplacementSS()
        
        # create local variables of reciprocal displacement and tractions
        Uss = funda.DisplacementSS
        
        # Calculate singular part of G matrix
        elm_par += jac*gauss_wt*Ni
        G += np.dot(Uss,elm_par)
        
    # non-singular part
    for xsi in gauss_pt:
        eta = (1+xsi)/2 #eta = (1+xsi)/2 => xsi = 2*eta-1 # fro start node 
        # Calculate element parameters
        Ni = serendip_1D(xsi, inci) # Shape function
        jac, Un = jac_normal_1D(xsi, inci, elcoord) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord) # Cartesian coordinate at the integration point

        # Calculate radius r and derivative dx/dr = [dx1/dr, dx2/dr] between collocation and integration points
        r, dxr = radius_param(ncoord, int_pt)
        drn = np.dot(dxr, Un) # dr/dn
        
        # Calculate reciprocal displacement and tractions for normal case (no singularity)
        funda.calculateDisplacementSN(r, dxr, eta)
        funda.calculateTraction_dyn(r, dxr, Un, drn)
        funda.calculateTraction_st(r, dxr, Un, drn)
        
        # create local variables of reciprocal displacement and tractions
        Uf = funda.displacement
        Tdf = funda.Traction_dyn
        Ts = funda.Traction_st
        
        # Calculate G and H matrices
        elm_par += jac*gauss_wt*Ni
        G += np.dot(Uf,elm_par)
        Hs += np.dot(Tdf,elm_par)
        H += np.dot(Ts,elm_par)
        
    return G, Hs, H
        
        
def bem_elm_eldyn_end():
    '''
    Computes element BEM matrices for each element node pair
    for collocation node end node of the quadratic element
    '''
    # Reset/initialize the local G and H matrices for new element
    G = 0.0
    Hs = 0.0
    H = 0.0

    # singular part
    for eta in gauss_pt_log: 
        xsi = 2*eta+1 #eta = (xsi-1)/2 => xsi = 2*eta+1 # for end node 

        # Calculate element parameters
        Ni = serendip_1D(xsi, inci) # Shape function
        jac, Un = jac_normal_1D(xsi, inci, elcoord) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord) # Cartesian coordinate at the integration point

        # Calculate radius r and derivative dx/dr = [dx1/dr, dx2/dr] between collocation and integration points
        r, dxr = radius_param(ncoord, int_pt)
        drn = np.dot(dxr, Un) # dr/dn
        
        # Calculate reciprocal displacement and tractions for normal case (no singularity)
        funda.calculateDisplacementSS()
        
        # create local variables of reciprocal displacement and tractions
        Uss = funda.DisplacementSS
        
        # Calculate singular part of G matrix
        elm_par += jac*gauss_wt*Ni
        G += np.dot(Uss,elm_par)
        
    # non-singular part
    for xsi in gauss_pt:
        eta = (xsi-1)/2 #eta = (xsi-1)/2 => xsi = 2*eta+1 # fro start node 
        # Calculate element parameters
        Ni = serendip_1D(xsi, inci) # Shape function
        jac, Un = jac_normal_1D(xsi, inci, elcoord) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord) # Cartesian coordinate at the integration point

        # Calculate radius r and derivative dx/dr = [dx1/dr, dx2/dr] between collocation and integration points
        r, dxr = radius_param(ncoord, int_pt)
        drn = np.dot(dxr, Un) # dr/dn
        # Calculate element parameters
        Ni = serendip_1D(xsi, inci) # Shape function
        jac, Un = jac_normal_1D(xsi, inci, elcoord) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord) # Cartesian coordinate at the integration point
        # Calculate reciprocal displacement and tractions for normal case (no singularity)
        funda.calculateDisplacementSN(r, dxr, eta)
        funda.calculateTraction_dyn(r, dxr, Un, drn)
        funda.calculateTraction_st(r, dxr, Un, drn)
        
        # create local variables of reciprocal displacement and tractions
        Uf = funda.displacement
        Tdf = funda.Traction_dyn
        Ts = funda.Traction_st
        
        # Calculate G and H matrices
        elm_par += jac*gauss_wt*Ni
        G += np.dot(Uf,elm_par)
        Hs += np.dot(Tdf,elm_par)
        H += np.dot(Ts,elm_par)
        
    return G, Hs, H

def bem_elm_eldyn_int():
    '''
    Computes element BEM matrices for each element node pair
    for collocation node intermediate node of the quadratic element
    '''
            
    # The element should be divided into two parts
    # 0----A------2-------B------1
    # lets divide the element into 
        
        