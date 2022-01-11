'''
Created on Jan 11, 2022

@author: mbb
'''


from scipy.special._ufuncs import kn
import math, cmath

EUGAM = math.pi/(2.0*math.exp(1.0)) # Euler Gamma Number

def fundaEl():
    '''
    Calculates fundamental solution for elastodynamics
    '''
    pass


def psi_kappa_plane(omega, r, small_arg=False):
    '''
    psi and kappa functions and their derivatives for 2D inplane case
    Dominguez equations (2.263), (2.264), (2.266) and (2.267)
    '''
    
    z1 = 1.0j*omega*r/c1
    z2 = 1.0j*omega*r/c2
    
    if small_arg:
        # Modified bessel functios of second kind for small arguments (near singular cases)
        K0_z1 = kn(0, z1) - (-cmath.log(0.5*z1)-EUGAM)
        K1_z1 = kn(1, z1) - ((1.0/z1)+(0.5*z1)*(cmath.log(0.5*z1)+EUGAM-0.5))
        K2_z1 = kn(2, z1) - ((2.0/(z1*z1))-0.5)
        
        K0_z2 = kn(0, z2) - (-cmath.log(0.5*z2)-EUGAM)
        K1_z2 = kn(1, z2) - ((1.0/z2)+(0.5*z2)*(cmath.log(0.5*z2)+EUGAM-0.5))
        K2_z2 = kn(2, z2) - ((2.0/(z2*z2))-0.5)
        
        AA = (z2*z2/(2.0*r))*(cmath.log(0.5*z2)+EUGAM-0.5)
        BB = (c2*c2/(c1*c1))*(z1*z1/(2.0*r))*(cmath.log(0.5*z1)+EUGAM-0.5)
        
        dK = 0.5 -0.5*(c2*c2/(c1*c1))
        
    else:
        # Modified bessel functios of second kind
    
        K0_z1 = kn(0, z1)
        K1_z1 = kn(1, z1)
        K2_z1 = kn(2, z1)
        
        K0_z2 = kn(0, z2)
        K0_z2 = kn(0, z2)
        K2_z2 = kn(2, z2)
        
        AA = 0.0
        BB = 0.0
        dK = 0.0 
        
    
    # psi and kappa functions
    psi = K0_z2 + (c2/(1.0j*omega*r))*(K1_z2 - (c2/c1)*K1_z1)
   
    
    
    # This is used for non singular part of singular case of G matrix  in singular case                 
    self.PSIss = complex128(0.5*((K2_z2-0.5)-(c2*c2/(c1*c1)) *((K2_z1-0.5)) +(K0_z2 -EUGAM) \
                                    +(c2*c2/(c1*c1))*(K0_z1-EUGAM)))\
                                -0.5*(cmath.log(0.5*z2)+(c2*c2/(c1*c1))*cmath.log(self.k1/2))
    
    kappa = K2_z2 - (c2*c2/(c1*c1))*K2_z1 + dK
    
    
    # derivatives of psi and kappa functions wrt r dpsi/dr and dkappa/dr
    dpsi = 1.0j*omega*K1_z2*(2.0*c2/(omega*omega*r*r) - (1.0/c2)) \
            + K1_z1*(2.0*c2*c2/(c1*1.0j*omega*r*r)) \
                - K0_z2/r + K0_z1*(c2*c2/(c1*c1*r)) \
                    - AA
                
    dkappa = - K1_z2*1.0j*omega/c2 -(2.0/r)*K2_z2 \
        + (c2*c2*1.0j*omega/(c1*c1*c1))* K1_z1 + (2.0*c2*c2/(r*c1*c1))*K2_z1 \
            -AA + BB

    