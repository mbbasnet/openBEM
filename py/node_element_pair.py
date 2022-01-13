import numpy as np
import isoparametric
from numpy.polynomial.legendre import leggauss

from openBEM.py.isoparametric import jac_normal_1D, serendip_1D

from openBEM.py.gauss_quadrature import GaussQuadrature

ncoord = np.array([10, 20])
elcoord = np.array(np.array([[-1, 0], [1,0], [0,3]]))
inci = 2 # quadratic element

#funda = FundamentalSolution(c1, c2, mu, nu, omega) # fundamental solution for elastic isotropic case

def radius_param(col_pt, int_pt):

    dxr = int_pt - col_pt # differences dx
    r = np.linalg.norm(dxr) # radius as length of the vector
    dxr = dxr/r # normalize cx/dr
    return r, dxr

def logarithmicQuadrature(deg):
    '''
    '''
    
    if deg == 2:
        logGaussPoint = np.array([0.112008806166976, 0.602276908118738], dtype = float64)
        logGaussWeight = np.array([0.718539319030385,0.281460680969616], dtype = float64)
        
    elif deg == 3:
        logGaussPoint = np.array([0.638907930873254e-1, 0.368997063715619, \
                                    0.766880303938942], dtype = float64)
        logGaussWeight = np.array([0.513404552232363, 0.391980041201488, \
                                    0.946154065661491e-1], dtype = float64)
        
    elif deg == 4:
        logGaussPoint = np.array([0.414484801993832,0.245274914320602,0.556165453560276,\
                                    0.848982394532985], dtype = float64)
        logGaussWeight = np.array([0.383464068145135,0.386875317774763,0.190435126950142,\
                                    0.392254871299598e-1], dtype = float64)
        
    elif deg == 5:
        logGaussPoint = np.array([0.291344721519721e-1,0.173977213320898,0.411702520284902,\
                                    0.677314174582820,0.894771361031008], dtype = float64)
        logGaussWeight = np.array([0.297893471782894,0.349776226513224,0.234488290044052,\
                                    0.989304595166331e-1,0.189115521431958e-1], dtype = float64)
        
    elif deg == 6:
        logGaussPoint = np.array([0.216340058441169e-1,0.129583391154951,0.314020449914766\
                                    ,0.538657217351802,0.756915337377403,0.922668851372120], dtype = float64)
        logGaussWeight = np.array([0.238763662578548,0.308286573273947,0.245317426563210,\
                                    0.142008756566477,0.554546223248863e-1,0.101689586929323e-1], dtype = float64)
        
    elif deg == 7:
        logGaussPoint = np.array([0.167193554082585e-1,0.100185677915675,0.246294246207931,0.433463493257033,\
                                    0.632350988047766,0.811118626740106,0.940848166743348], dtype = float64)
        logGaussWeight = np.array([0.196169389425248,0.270302644247273,0.239681873007691,0.165775774810433,\
                                    0.889432271376580e-1,0.331943043565711e-1,0.593278701512593e-2], dtype = float64)
        
    elif deg == 8:
        logGaussPoint = np.array([0.133202441608925e-1,0.797504290138949e-1,0.197871029326188,0.354153994351909,\
                                    0.529458575234917,0.701814529939100,0.849379320441107,0.953326450056360], dtype = float64)
        logGaussWeight = np.array([0.164416604728003,0.237525610023306,0.226841984431919,0.175754079006070,\
                                    0.112924030246759,0.578722107177821e-1,0.209790737421330e-1,0.368640710402762e-2], dtype = float64)
        
    elif deg == 9:
        logGaussPoint = np.array([0.108693360841755e-1,0.649836663380079e-1,0.162229398023883,0.293749903971675,0.446631881905468,\
                                    0.605481662776129,0.754110137157164,0.877265828835838,0.962250559410282], dtype = float64)
        logGaussWeight = np.array([0.140068438748135,0.209772205201030,0.211427149896603,0.177156233938080,0.127799228033206,\
                                    0.784789026115622e-1,0.390225049853991e-1,0.138672955495930e-1,0.240804103639231e-2], dtype = float64)
        
    elif deg == 10:
        logGaussPoint = np.array([0.904263096219965e-2,0.539712662225006e-1,0.135311824639251,0.247052416287160,\
                                    0.380212539609332,0.523792317971843,0.665775205516425,0.794190416011966,0.898161091219004,\
                                    0.968847988718634], dtype = float64)
        '''
        logGaussWeight = np.array([0.120955131954571,0.186363542564072,0.195660873277760,0.173577142182907,0.135695672995484,\
                                    0.936467585381105e-1,0.557877273514159e-1,0.271598108992333e-1,0.951518260284852e-2,\
                                    0.163815763359826e-2], dtype = float64)
        '''
        logGaussWeight = np.array([0.1209551319,0.1863635425,0.1956608732,0.1735771421,0.1356956729,0.0936467585,\
                                    0.0557877273,0.0271598109,0.0095151826,0.0016381576], dtype = float64)
        
    elif deg == 12:
        logGaussPoint = np.array([0.654872227908006e-2,0.389468095604500e-1,0.981502631060066e-1,0.181138581590632,\
                                    0.283220067667373,0.398434435163437,0.519952626792353,0.640510916716107,\
                                    0.752865012051831,0.850240024162302,0.926749683223914,0.977756129689998], dtype = float64)
        logGaussWeight = np.array([0.931926914439313e-1,0.149751827576322,0.166557454364593,0.159633559436988,\
                                    0.138424831864836,0.110016570635721,0.799618217708290e-1,0.524069548246418e-1,\
                                    0.300710888737612e-1,0.142492455879983e-1,0.489992458232176e-2,0.834029038056903e-3], dtype = float64)
        
    elif deg == 14:
        logGaussPoint = np.array([0.496600357386854e-2,0.294325401188852e-1,0.743762922245358e-1,0.138138491989186,0.218055648498959,\
                                    0.310662083918102,0.411872475177750,0.517179307398654,0.621864859728511,0.721220745208109,\
                                    0.810765988071590,0.886454038034435,0.944859139461819,0.983331026485678], dtype = float64)
        logGaussWeight = np.array([0.742912250675104e-1,0.122988772469323,0.142199306562523,0.143229297641264,0.132345083772085,\
                                    0.114135875736676,0.922830380790736e-1,0.697536732939376e-1,0.488303236005136e-1,0.311017960644161e-1,\
                                    0.174628119501961e-1,0.814242342987594e-2,0.276843641856394e-2,0.467935914040560e-3], dtype = float64)
        
    elif deg == 16:
        logGaussPoint = np.array([0.389783448711592e-2,0.230289456168732e-1,0.582803983062404e-1,0.108678365091054,\
                                    0.172609454909844,0.247937054470578,0.332094549129917,0.422183910581949,\
                                    0.515082473381463,0.607556120447729,0.696375653228214,0.778432565873265,\
                                    0.850850269715391,0.911086857222272,0.957025571703542,0.987047800247984], dtype = float64)
        logGaussWeight = np.array([0.607917100435912e-1,0.102915677517582,0.122355662046009,0.127569246937016,\
                                    0.123013574600071,0.111847244855486,0.965963851521243e-1,0.793566643514731e-1,\
                                    0.618504945819652e-1,0.454352465077267e-1,0.310989747515818e-1,0.194597659273608e-1,\
                                    0.107762549632055e-1,0.497254289008764e-2,0.167820111005119e-2,0.282353764668436e-3], dtype = float64)
        
    return logGaussPoint,logGaussWeight

        
def bem_elm_eldyn(ncoord, elcoord, inci, funda, ngauss):
    '''
    Computes element BEM matrices for each element node pair
    for normal case i.e. collocation node not in the element
    '''
    # Reset/initialize the local G and H matrices for new element
    G = 0.0
    Hs = 0.0
    H = 0.0
    
    gauss_pt, gauss_wt = leggauss(ngauss)

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
        elm_par = jac*gauss_wt*Ni
        G += np.dot(Uf,elm_par)
        Hs += np.dot(Tdf,elm_par)
        H += np.dot(Ts,elm_par)
    
    return G, Hs, H
        
def bem_elm_eldyn_st(ncoord, elcoord, inci, funda, ngauss):
    '''
    Computes element BEM matrices for each element node pair
    for collocation node start node of the quadratic element
    '''
    # Reset/initialize the local G and H matrices for new element
    G = 0.0
    Hs = 0.0
    H = 0.0
    
    # singular part
    gauss_pt_log, gauss_wt_log = logarithmicQuadrature(ngauss)
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
        elm_par = jac*gauss_wt_log*Ni
        G += np.dot(Uss,elm_par)
        
    # non-singular part
    gauss_pt, gauss_wt = leggauss(ngauss)
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
        elm_par = jac*gauss_wt*Ni
        G += np.dot(Uf,elm_par)
        Hs += np.dot(Tdf,elm_par)
        H += np.dot(Ts,elm_par)
        
    return G, Hs, H
        
        
def bem_elm_eldyn_end(ncoord, elcoord, inci, funda, ngauss):
    '''
    Computes element BEM matrices for each element node pair
    for collocation node end node of the quadratic element
    '''
    # Reset/initialize the local G and H matrices for new element
    G = 0.0
    Hs = 0.0
    H = 0.0

    # singular part
    gauss_pt_log, gauss_wt_log = logarithmicQuadrature(ngauss)
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
        elm_par = jac*gauss_wt_log*Ni
        G += np.dot(Uss,elm_par)
        
    # non-singular part
    gauss_pt, gauss_wt = leggauss(ngauss)
    for xsi in gauss_pt:
        eta = (xsi-1)/2 #eta = (xsi-1)/2 => xsi = 2*eta+1 # fro start node 
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
        elm_par = jac*gauss_wt*Ni
        G += np.dot(Uf,elm_par)
        Hs += np.dot(Tdf,elm_par)
        H += np.dot(Ts,elm_par)
        
    return G, Hs, H


def bem_elm_eldyn_int(ncoord, elcoord, inci, funda, ngauss):
    '''
    Computes element BEM matrices for each element node pair
    for collocation node intermediate node of the quadratic element
    '''

    # Reset/initialize the local G and H matrices for new element
    G_a = 0.0
    Hs_a = 0.0
    H_a = 0.0

    G_b = 0.0
    Hs_b = 0.0
    H_b = 0.0

            
    # The element should be divided into two parts
    # 0----A------2-------B------1
    # lets divide the element into Two 

    # Dividing the quadratic element into two linear elements
    # Coordinates for element 0----(A)---2
    elcoord_A = np.zeros((2,2))
    elcoord_A[0][:] = elcoord[0][:]
    elcoord_A[1][:] = elcoord[2][:]

    # Coordinates for element 0----(B)---2
    elcoord_B = np.zeros((2,2))
    elcoord_B[0][:] = elcoord[2][:]
    elcoord_B[1][:] = elcoord[1][:]

    # singular part
    gauss_pt_log, gauss_wt_log = logarithmicQuadrature(ngauss)
    for eta in gauss_pt_log: # element A
        xsi = 2*eta+1 #eta = (xsi-1)/2 => xsi = 2*eta+1 # for end node 

        # Calculate element parameters
        Ni = serendip_1D(xsi, inci-1) # Shape function for linear element
        jac, Un = jac_normal_1D(xsi, inci-1, elcoord_A) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord_A) # Cartesian coordinate at the integration point

        # Calculate radius r and derivative dx/dr = [dx1/dr, dx2/dr] between collocation and integration points
        r, dxr = radius_param(ncoord, int_pt)
        drn = np.dot(dxr, Un) # dr/dn
        
        # Calculate reciprocal displacement and tractions for normal case (no singularity)
        funda.calculateDisplacementSS()
        
        # create local variables of reciprocal displacement and tractions
        Uss = funda.DisplacementSS
        
        # Calculate singular part of G matrix
        elm_par = jac*gauss_wt_log*Ni
        G_a += np.dot(Uss,elm_par)
        
    # non-singular part
    gauss_pt, gauss_wt = leggauss(ngauss)
    for xsi in gauss_pt: # Element A
        eta = (xsi-1)/2 #eta = (xsi-1)/2 => xsi = 2*eta+1 # fro start node 
        # Calculate element parameters
        Ni = serendip_1D(xsi, inci-1) # Shape function
        jac, Un = jac_normal_1D(xsi, inci-1, elcoord_A) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord_A) # Cartesian coordinate at the integration point

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
        elm_par = jac*gauss_wt*Ni
        G_a += np.dot(Uf,elm_par)
        Hs_a += np.dot(Tdf,elm_par)
        H_a += np.dot(Ts,elm_par)

    # singular part (Element B)
    for eta in gauss_pt_log: # Element B
        xsi = 2*eta-1 #eta = (1+xsi)/2 => xsi = 2*eta-1 # for start node

        # Calculate element parameters
        Ni = serendip_1D(xsi, inci-1) # Shape function
        jac, Un = jac_normal_1D(xsi, inci-1, elcoord_B) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord_B) # Cartesian coordinate at the integration point

        # Calculate radius r and derivative dx/dr = [dx1/dr, dx2/dr] between collocation and integration points
        r, dxr = radius_param(ncoord, int_pt)
        drn = np.dot(dxr, Un) # dr/dn
        
        # Calculate reciprocal displacement and tractions for normal case (no singularity)
        funda.calculateDisplacementSS()
        
        # create local variables of reciprocal displacement and tractions
        Uss = funda.DisplacementSS
        
        # Calculate singular part of G matrix
        elm_par = jac*gauss_wt_log*Ni
        G_b += np.dot(Uss,elm_par)
        
    # non-singular part (Element B)
    for xsi in gauss_pt:
        eta = (1+xsi)/2 #eta = (1+xsi)/2 => xsi = 2*eta-1 # fro start node 
        # Calculate element parameters
        Ni = serendip_1D(xsi, inci-1) # Shape function
        jac, Un = jac_normal_1D(xsi, inci-1, elcoord_B) # Jacobian and unit normal
        int_pt = np.matmul(Ni, elcoord_B) # Cartesian coordinate at the integration point

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
        elm_par = jac*gauss_wt*Ni
        G_b += np.dot(Uf,elm_par)
        Hs_b += np.dot(Tdf,elm_par)
        H_b += np.dot(Ts,elm_par)
    
    # Combining the Variables from two matrices
    G = np.zeros((2,3))
    Hs = np.zeros((2,3))
    H = np.zeros((2,3))

    G[:][0] = G_a[:][0]
    G[:][1] = G_a[:][1] + G_b[:][0]
    G[:][2] = G_b[:][1]

    Hs[:][0] = Hs_a[:][0]
    Hs[:][1] = Hs_a[:][1] + Hs_b[:][0]
    Hs[:][2] = Hs_b[:][1]

    H[:][0] = H_a[:][0]
    H[:][1] = H_a[:][1] + H_b[:][0]
    H[:][2] = H_b[:][1]
        
    return G, Hs, H
     
        