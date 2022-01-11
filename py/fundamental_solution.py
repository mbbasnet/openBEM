'''
Created on Sep 28, 2015

@author: MBB
'''
import numpy as np
from numpy import complex128, float64
from scipy.special._ufuncs import kv
import math, cmath

from BEM_system.bessel import BesselFunction
bessel = BesselFunction()

KroneckerDelta = lambda a,b:1 if int(round(a)) == int(round(b)) else 0 # define kroneckerDelta for all classes

class FundamentalSolution(object):
    '''
    Contains methods to calculate fundamental solutions for Traction and Displacement
    calculates PSI, KAPPA and U* for singular and non singular case and P*_static
    '''

    def __init__(self, c1, c2, mu, nu, angularFrequency):
        '''
        Constructor
        '''
        
        self.pWaveVelocity = c1
        self.sWaveVelocity = c2
        self.shearModulus = mu
        self.poissonRatio = nu
        self.omega = angularFrequency
        self.alpha = 2. # for 2D case alpha is 2
        
        self.eulerGamma = float64(math.pi/(2*math.exp(1)))        
        #self.deltaKAPPA = complex128(0.5*(1-self.sWaveVelocity**2/self.pWaveVelocity**2))
        #self.AA = lambda z,r: complex128((z**2/(2*r))*(cmath.log(z/2)+self.eulerGamma-1/2))
        #self.BB = lambda z,r: complex128((self.sWaveVelocity/self.pWaveVelocity)**2*(z**2/(2*r))\
        #                                 *(cmath.log(z/2)+self.eulerGamma-1/2))
        
    def calculateBasicParameters(self,r):
        '''
        '''
        self.k1 = complex128(-self.omega*1j/self.pWaveVelocity) # Dominguez page 130
        self.k2 = complex128(-self.omega*1j/self.sWaveVelocity) # Dominguez page 130
        self.z1 = self.k1*r 
        self.z2 = self.k2*r 
        self.VelRatio = self.pWaveVelocity/self.sWaveVelocity
        
        self.AA = complex128((self.z2**2/(2*r))*(cmath.log(self.z2/2)+self.eulerGamma-1/2))
        self.BB = complex128((self.sWaveVelocity/self.pWaveVelocity)**2*(self.z1**2/(2*r))\
                                         *(cmath.log(self.z1/2)+self.eulerGamma-1/2))
        self.deltaKAPPA = 0.5*(1-(1/self.VelRatio)**2)
        
        '''
        self.FK0_Z1,self.FK1_Z1,self.FK2_Z1,\
        self.FK0M_Z1,self.FK1M_Z1,self.FK2M_Z1 = bessel.besselK123(self.k1*r)
        self.FK0_Z2,self.FK1_Z2,self.FK2_Z2,\
        self.FK0M_Z2,self.FK1M_Z2,self.FK2M_Z2 = bessel.besselK123(self.k2*r)
        '''
        self.FK0_Z1 = kv(0,self.k1*r)
        self.FK1_Z1 = kv(1,self.k1*r)
        self.FK2_Z1 = kv(2,self.k1*r)
        self.FK0_Z2 = kv(0,self.k2*r)
        self.FK1_Z2 = kv(1,self.k2*r)
        self.FK2_Z2 = kv(2,self.k2*r)
        
        self.FK0M_Z1 = self.FK0_Z1-(-cmath.log(self.z1/2)-self.eulerGamma)
        self.FK1M_Z1 = self.FK1_Z1-((1/self.z1)+(self.z1/2)*(cmath.log(self.z1/2)+self.eulerGamma-1/2))
        self.FK2M_Z1 = self.FK2_Z1-((2/self.z1**2)-(1/2))
        self.FK0M_Z2 = self.FK0_Z2-(-cmath.log(self.z2/2)-self.eulerGamma)
        self.FK1M_Z2 = self.FK1_Z2-((1/self.z2)+(self.z2/2)*(cmath.log(self.z2/2)+self.eulerGamma-1/2)) 
        self.FK2M_Z2 = self.FK2_Z2-((2/self.z2**2)-(1/2))
           
    def calculatePSI(self, r): # checked OK
        '''
        calculate PSI coefficient for 2D inplane wave propagation. Uses modified bessel functions of second kind
        singular and non-singular cases
        '''
        
        # bessel.besselK is modivied bessel function of second type
        # this PSI is used for calculation of G for non singular case that is r does not tend to zero
        self.PSI = complex128(self.FK0_Z2 \
                              + (self.FK1_Z2-self.FK1_Z1/self.VelRatio)/(self.k2*r)) # Domingueq eqn 2.209
        
        #print('PSI = ', self.PSI)
               
        ## the PSI derivative is used for computation of total traction only 
        self.dPSI_dr = complex128(-self.k2*self.FK1_Z2 + (self.k2*(-self.FK0_Z2/2 \
                           - self.FK2_Z2/2)- self.k1*(-self.FK0_Z1/2\
                            - self.FK2_Z1/2)/self.VelRatio)/(self.k2*r)- (self.FK1_Z2 \
                            - self.FK1_Z1/self.VelRatio)/(self.k2*r**2))
        #print('dPSI_dr = ', self.dPSI_dr)
                           
        # This is used for non singular part of singular case of G matrix  in singular case                 
        self.PSIss = complex128(0.5*((self.FK2M_Z2-0.5)-(1/self.VelRatio)**2\
                                *((self.FK2M_Z1-0.5)) +(self.FK0M_Z2\
                                -self.eulerGamma)+(1/self.VelRatio)**2*(self.FK0M_Z1-self.eulerGamma)))\
                                -0.5*(cmath.log(self.k2/2)+(1/self.VelRatio)**2*cmath.log(self.k1/2))

        # this singular term is treated separately  #-0.5*(1+(c2/c1)**2)*cmath.log(r) 
        #self.PSI = self.PSIss-0.5*(1+(1/self.VelRatio)**2)*cmath.log(r)
        #print('PSI nonsingular part = ', self.PSIss)
                                             
        # This dPSIss_dr is used for Dynamic part of H matrix only        
        self.dPSIss_dr = complex128(-self.k2*self.FK1M_Z2 + (self.k2*(-self.FK0M_Z2/2\
                             - self.FK2M_Z2/2) - self.k1*(-self.FK0M_Z1/2\
                             - self.FK2M_Z1/2)/self.VelRatio)/(self.k2*r)- (self.FK1M_Z2 \
                             - self.FK1M_Z1/self.VelRatio)/((self.k2)*r**2) -self.AA)
                            # Singular part, omitted from calculation self.PSIss = self.dPSIss_dr - (1+1/self.VelRatio**2)/(2*r)
                            # this singular part is acocounted in static part of H matrix
                            
                            
        #print('dPSI_dr nonsingular part = ',self.dPSIss_dr)
                                                
    def calculateKAPPA(self, r): # checked OK
        '''
        calculates kappa for real r. singular and non-singular cases
        '''
        
        # This KAPPAss is used for dynamic part of H matrix only
        self.KAPPAss = complex128(self.FK2M_Z2-self.FK2M_Z1/self.VelRatio**2) # Dominquez eqn 2.209 
        self.KAPPA = complex128(self.KAPPAss-self.deltaKAPPA)
        
        #print('KAPPA nonsingular part = ', self.KAPPAss)

        # This dKAPPA_dr is used for computation of Dynamic part of H matrix
        self.dKAPPA_dr = complex128((-self.k2*self.FK1M_Z2-2*self.FK2M_Z2/r\
                                  +(self.k1*self.FK1M_Z1+2*self.FK2M_Z1/r)/self.VelRatio**2)\
                                    -(self.AA-self.BB))                         
        #print('dKAPPA_dr nonsingular part = ',self.dKAPPAss_dr)
        
    def calculateDisplacement(self, rDiff): # checked OK
        '''
        calculates reicprocal displacement U* using already calculated PSI and KAPPA
        r - distance from the node to collocation point
        rDiff - vector of derivatives of r 
        eta - unit normal
        l,k - displacement direction of cartesian coordinate
        
        '''     
        #omega = self.omega
        # Dominguez equation 2.199
        TERM1 = complex128(1/(self.alpha*math.pi*self.shearModulus))
        
        U = np.zeros((2,2),dtype= complex128)
        '''
        for l in range(0,2):
            for k in range(0,2):
                U[l][k] = complex128(TERM1*(self.PSI*KroneckerDelta(l,k)\
                                            -self.KAPPA*rDiff[l]*rDiff[k]))
        '''
        U[0][0]= complex128(TERM1*(self.PSI-self.KAPPA*rDiff[0]*rDiff[0])) 
        U[0][1]= complex128(-TERM1*self.KAPPA*rDiff[0]*rDiff[1]) 
        U[1][0]= U[0][1]
        U[1][1]= complex128(TERM1*(self.PSI-self.KAPPA*rDiff[1]*rDiff[1]))  
        
        self.Displacement = np.asmatrix(U, dtype= complex128)
        
    def calculateDisplacementSN(self,r, rDiff, NDZ): # checked OK
        '''
        calculates reicprocal displacement U* using already calculated PSI and KAPPA
        r - distance from the node to collocation point
        rDiff - vector of derivatives of r 
        eta - unit normal
        l,k - displacement direction of cartesian coordinate
        NDZ is dependent upon collocation point and local coordinates
        ''' 
    
        R = r/NDZ # NDZ is dependent upon collocation node 
        # for collocation point at 1, NDZ = abs(1+zeta)/2,   for 2: NDZ = abs(zeta )  for 3: NDZ = abs(1-zeta)/2
        # Adding non singular part as a part of partition of singular part log(r) = log(R)-log(1/NDZ)
        #-0.5*(1+(c2/c1)**2)*cmath.log(r) this term divided into R and NDZ
        PSI_sn = complex128(self.PSIss -0.5*(1+(1/self.VelRatio)**2)*cmath.log(R))
        KAPPA_sn = self.KAPPAss - self.deltaKAPPA # to avoid numerical errors
        TERM1 = complex128(1/(self.alpha*math.pi*self.shearModulus))
        U = np.zeros((2,2),dtype= complex128)
        #print('Kappa', self.KAPPA, r)
        '''
        for l in range(0,2):
            for k in range(0,2):
                # PSIss has only non singular part . singular part is dealt separately
                U[l][k] = complex128(TERM1*(PSI_sn*KroneckerDelta(l,k)-KAPPA_sn*rDiff[l]*rDiff[k]))
        '''
        U[0][0]= complex128(TERM1*(PSI_sn-KAPPA_sn*rDiff[0]*rDiff[0])) 
        U[0][1]= complex128(-TERM1*KAPPA_sn*rDiff[0]*rDiff[1]) 
        U[1][0]= complex128(-TERM1*KAPPA_sn*rDiff[1]*rDiff[0]) 
        U[1][1]= complex128(TERM1*(PSI_sn-KAPPA_sn*rDiff[1]*rDiff[1]))                
        self.Displacement = np.asmatrix(U, dtype= complex128)
        
    def calculateDisplacementSS(self): # 
        '''
        calculates reicprocal displacement U* using already calculated PSI and KAPPA
        R = r/NDZ where ndz is local coordinates for singular cases of transform
        point
        
        '''   
        #Singular term is -0.5*(1+(c2/c1)**2)*cmath.log(1/NDZ)
        # then the log r term is again divided into two nonsingular and non singular terms R and (1+zeta)/2
        # then logarithmic gauss quadrature scheme can be applied on it.
        
        PSI_ss = 0.5*(1+(self.sWaveVelocity/self.pWaveVelocity)**2) # the term log(1/NDZ) is contained in the formula of integration
        # the sign is +ve as -ve sign *-log(1/NDZ) make is positive
        TERM1 = complex128(1/(self.alpha*math.pi*self.shearModulus))
        '''
        U_ss = np.zeros((2,2),dtype= complex128)
        for l in range(0,2):
            for k in range(0,2):
                U_ss[l][k] = complex128(TERM1*PSI_ss*KroneckerDelta(l,k))

        self.DisplacementSS = np.asmatrix(U_ss, dtype= complex128)
        '''
        U_ss = complex128(TERM1*PSI_ss)
        self.DisplacementSS = np.asmatrix([[U_ss,0],[0,U_ss]], dtype= complex128)
                    
    def calculateTraction_st(self,r,rDiff,eta,dr_dn): # Checked OK
        '''
        calculates P* for static case (used for calculation of H_static )
        '''

        S2 = -1/(4*math.pi*(1-self.poissonRatio)*r)
        S1 = 1-2*self.poissonRatio
        P_st = np.zeros((2,2),dtype= complex128)
        
        '''
        for l in range(0,2):
            for k in range(0,2):
                P_st[l][k]= complex128(S2*(dr_dn*(S1*KroneckerDelta(k,l)+2*rDiff[k]*rDiff[l])\
                                                          + S1*(eta[l]*rDiff[k]-eta[k]*rDiff[l])))
   
        '''       
        
        P_st[0][0]= complex128(S2*dr_dn*(S1+2*rDiff[0]*rDiff[0]))
        P_st[0][1]= complex128(S2*(dr_dn*(2*rDiff[1]*rDiff[0])+ S1*(eta[0]*rDiff[1]-eta[1]*rDiff[0])))
        P_st[1][0]= complex128(S2*(dr_dn*(2*rDiff[0]*rDiff[1])+ S1*(eta[1]*rDiff[0]-eta[0]*rDiff[1])))
        P_st[1][1]= complex128(S2*dr_dn*(S1+2*rDiff[1]*rDiff[1]))
        
        
        self.Traction_st = np.asmatrix(P_st, dtype = complex128)
               
    def calculateTraction_dyn(self,r,rDiff,eta,dr_dn): # checked Ok.
        '''
        P* - for small arguments (r-->0)
        '''
         
        T0 = 1/(self.alpha*math.pi)
        T1 = (self.dPSIss_dr-self.KAPPAss/r)
        T2 = (2/r)*self.KAPPAss
        T3 = 2*self.dKAPPA_dr*dr_dn
        T4 = (self.dPSIss_dr-self.dKAPPA_dr-self.alpha*self.KAPPAss/(2*r))
        PR = np.zeros((2,2),dtype=complex128)
        '''
        for l in range(0,2):
            for k in range(0,2):
               
                PR[l][k] = complex128((1/(self.alpha*math.pi))\
                                     *((self.dPSIss_dr-self.KAPPAss/r)*(KroneckerDelta(l,k)*dr_dn+rDiff[k]*eta[l])\
                                        -(2/r)*self.KAPPAss*(eta[k]*rDiff[l]-2*rDiff[l]*rDiff[k]*dr_dn)\
                                        -2*self.dKAPPA_dr*rDiff[l]*rDiff[k]*dr_dn+((self.VelRatio**2)-2)\
                                        *(self.dPSIss_dr-self.dKAPPA_dr-self.alpha*self.KAPPAss/(2*r))*rDiff[l]*eta[k]))
        '''        
               
        PR[0][0] = complex128(T0\
                            *(T1*(dr_dn+rDiff[0]*eta[0])\
                            -T2*(eta[0]*rDiff[0]-2*rDiff[0]*rDiff[0]*dr_dn)\
                            -T3*rDiff[0]*rDiff[0]+((self.VelRatio**2)-2)\
                            *T4*rDiff[0]*eta[0]))
        
        PR[0][1] = complex128(T0\
                            *(T1*(rDiff[1]*eta[0])\
                            -T2*(eta[1]*rDiff[0]-2*rDiff[0]*rDiff[1]*dr_dn)\
                            -T3*rDiff[0]*rDiff[1]+((self.VelRatio**2)-2)\
                            *T4*rDiff[0]*eta[1]))
        
        PR[1][0] = complex128(T0\
                            *(T1*(rDiff[0]*eta[1])\
                            -T2*(eta[0]*rDiff[1]-2*rDiff[1]*rDiff[0]*dr_dn)\
                            -T3*rDiff[1]*rDiff[0]+((self.VelRatio**2)-2)\
                            *T4*rDiff[1]*eta[0]))
        
        PR[1][1] = complex128(T0\
                            *(T1*(dr_dn+rDiff[1]*eta[1])\
                            -T2*(eta[1]*rDiff[1]-2*rDiff[1]*rDiff[1]*dr_dn)\
                            -T3*rDiff[1]*rDiff[1]+((self.VelRatio**2)-2)\
                            *T4*rDiff[1]*eta[1]))
        
        
        self.Traction_dyn = np.asmatrix(PR, dtype = complex128)
        
       
    def calculateTraction(self,r,rDiff,eta,dr_dn): # checked Ok. not used as we calculate dynamic part separately
        '''
        P* - for small arguments (r-->0)
        NOT USED SO FAR FOR THIS PURPOSE HERE
        '''
        
        T0 = 1/(self.alpha*math.pi)
        T1 = (self.dPSI_dr-self.KAPPA/r)
        T2 = (2/r)*self.KAPPA
        T3 = 2*self.dKAPPA_dr*dr_dn
        T4 = (self.dPSI_dr-self.dKAPPA_dr-self.alpha*self.KAPPA/(2*r))
        PR = np.zeros((2,2),dtype=complex128)
        
        #for l in range(0,2):
        #    for k in range(0,2):

        #       PR[l][k] = complex128((1/(self.alpha*math.pi))\
        #                            *((self.dPSI_dr-self.KAPPA/r)*(KroneckerDelta(l,k)*dr_dn+rDiff[k]*eta[l])\
        #                            -(2/r)*self.KAPPA*(eta[k]*rDiff[l]-2*rDiff[l]*rDiff[k]*dr_dn)\
        #                            -2*self.dKAPPA_dr*rDiff[l]*rDiff[k]*dr_dn\
        #                            +((self.VelRatio**2)-2)*(self.dPSI_dr-self.dKAPPA_dr\
        #                            -self.alpha*self.KAPPA/(2*r))*rDiff[l]*eta[k]))
        
                
        PR[0][0] = complex128(T0\
                            *(T1*(dr_dn+rDiff[0]*eta[0])\
                            -T2*(eta[0]*rDiff[0]-2*rDiff[0]*rDiff[0]*dr_dn)\
                            -T3*rDiff[0]*rDiff[0]+((self.VelRatio**2)-2)\
                            *T4*rDiff[0]*eta[0]))
        
        PR[0][1] = complex128(T0\
                            *(T1*(rDiff[1]*eta[0])\
                            -T2*(eta[1]*rDiff[0]-2*rDiff[0]*rDiff[1]*dr_dn)\
                            -T3*rDiff[0]*rDiff[1]+((self.VelRatio**2)-2)\
                            *T4*rDiff[0]*eta[1]))
        
        PR[1][0] = complex128(T0\
                            *(T1*(rDiff[0]*eta[1])\
                            -T2*(eta[0]*rDiff[1]-2*rDiff[1]*rDiff[0]*dr_dn)\
                            -T3*rDiff[1]*rDiff[0]+((self.VelRatio**2)-2)\
                            *T4*rDiff[1]*eta[0]))
        
        PR[1][1] = complex128(T0\
                            *(T1*(dr_dn+rDiff[1]*eta[1])\
                            -T2*(eta[1]*rDiff[1]-2*rDiff[1]*rDiff[1]*dr_dn)\
                            -T3*rDiff[1]*rDiff[1]+((self.VelRatio**2)-2)\
                            *T4*rDiff[1]*eta[1]))

        self.Traction = np.asmatrix(PR, dtype = complex128)
            
