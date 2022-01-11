#%%
import numpy as np


def serendip_1D(xsi, inci):
    '''
    Computes Shape functions for serindipity elements
    1D isoparametric representation
    '''
    if inci == 0: 
        # Constant element 
        # 0-----------1 # isoparametric coordinates
        # -----0------- # node location
        Ni = np.array([1])
    elif inci == 1:
        # Linear element 
        # 0-----------1 # isoparametric coordinates
        # 0-----------1 # node location
        Ni = np.array([0.5*(1-xsi), 0.5*(1+xsi)])
    elif inci == 2:
        # Quadratic element 
        # 0-----------1 # isoparametric coordinates
        # 0-----2-----1 # node location
        Ni = np.array([-0.5*xsi*(1-xsi), 0.5*xsi*(1+xsi), 1.0 - xsi*xsi])
    else:
        ValueError('The element type is not supported')
    return Ni # The shape function


def serendip_deriv_1D(xsi, inci):
    '''
    Computes derivatives of the Shape functions for serindipity elements
    1D isoparametric representation
    '''
    if inci == 0: 
        # Constant element 
        # 0-----------1 # isoparametric coordinates
        # -----0------- # node location
        dNi = np.array([0])
    elif inci == 1:
        # Linear element 
        # 0-----------1 # isoparametric coordinates
        # 0-----------1 # node location
        dNi = np.array([-0.5, 0.5])
    elif inci == 2:
        # Quadratic element 
        # 0-----------1 # isoparametric coordinates
        # 0-----2-----1 # node location
        dNi = np.array([-0.5+xsi, 0.5+xsi,  -2.0*xsi])
    else:
        ValueError('The element type is not supported')
    return dNi # The shape function

def jac_normal_1D(xsi, inci, coords):
    '''
    calculates jacobian and unit normal at the given isoparametric position
    '''
    # call derivative of shape function
    dNi = serendip_deriv_1D(xsi,inci)
    dNx = np.matmul(dNi, coords)
    
    # calculate unit normal for for 2d cartesian coordinates
    # for 1D isoparametric element with 2D coords plane 
    Un = np.array([-dNx[1], dNx[0]]) # This case corresponds to clockwise outward normal
    jac = np.linalg.norm(dNx) # Jacobian of the vector
    Un = Un/jac # Unit normal of the vector
    
    return jac, Un


'''
for xsi in np.arange(-1, 1, 0.1):
    #Ni = serendip_1D(ii, 1)
    coords = np.array([[-1, 0], [1,0], [0,3]])
    Ni = serendip_1D(xsi, 2)
    jac, Un = jac_normal_1D(xsi, 2, coords)
    ipt = np.matmul(Ni, coords)
    print('data:', xsi, jac, Un, ipt)
'''

# %%
