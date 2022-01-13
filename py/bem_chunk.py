'''
Compute BEM matrices for a chunk or a layer or a part
'''
import numpy as np
from openBEM.py.fundamental_solution import FundamentalSolution
from openBEM.py.node_element_pair import bem_elm_eldyn, bem_elm_eldyn_end, bem_elm_eldyn_int, bem_elm_eldyn_st

coord = np.zeros(100,2)
cnodes = np.arange(0,100)
celm = np.arange(0,100)
elcoord = np.zeros((3,2))
elm = np.arange(100,3)

if True:
    funda = FundamentalSolution(c1=2000.0, c2=700.0, mu=700*700*1800, nu=0.33, omega=6.28)
    for nn in cnodes: # list of current node indices
        ncoord = coord[nn]
        for ee in celm: # list of current element indices
            # for quadratic elements
            inci = 2
            # Getting element coordinates from the node indices
            elcoord[0][:] = coord[elm[ee][0]]
            elcoord[1][:] = coord[elm[ee][1]]
            elcoord[2][:] = coord[elm[ee][2]]

            # Check singularity and perform bem computation
            if nn == elm[ee][0]: # Start node singularity
                G_e, Hs_e, H_e = bem_elm_eldyn_st(ncoord,elcoord,inci,funda,ngauss=10)

            elif nn == elm[ee][1]: # End node singularity
                G_e, Hs_e, H_e = bem_elm_eldyn_end(ncoord,elcoord,inci,funda,ngauss=10)

            elif nn == elm[ee][2]: # Intermediate node singularity
                G_e, Hs_e, H_e = bem_elm_eldyn_int(ncoord,elcoord,inci,funda,ngauss=10)

            else: # No singularity
                G_e, Hs_e, H_e = bem_elm_eldyn(ncoord,elcoord,inci,funda,ngauss=10)

            # Map the eleG_e, Hs_e, H_e = ment matrices to BEM chunk level matrices G and H
                

