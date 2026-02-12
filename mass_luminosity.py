"""
Calculate Mass-Luminosity relations for different models of line
emission.

All functions take a vector of masses in M_sun and return luminosities
in L_sun.

Model parameter values are given by a dictionary called MLpar.  Each
function also takes a value of the redshift z even if the L(M) model is not
redshift dependent.  This allows the functions to be called consistently by
LineModel()

TODO:
Add in models from Matlab code
"""

import numpy as np
import astropy.units as u
import astropy.constants as cu
from scipy.interpolate import interp2d,interp1d

from limlam_mocker.limlam_mocker import add_log_normal_scatter

def MassPow(Mvec, MLpar, z):
    """
    Power law L(M)/L_sun = A*(M/M_sun)^b (See Breysse et al. 2015)

    Parameters:
    A	      Overall amplitude, dimensionless
    b         Power law slope, dimensionless
    
    Assumed to be redshift independent

    >>> Mvec = np.array([1e10,1e11,1e12]) * u.Msun
    >>> MLpar = {'A':2e-6, 'b':1.}
    >>> z = 3.0
    >>> print MassPow(Mvec,MLpar,z)
    [   20000.   200000.  2000000.] solLum
    """

    A = MLpar['A']
    b = MLpar['b']
    L = A * np.array(Mvec)**b*u.Lsun
    return L
    
def DblPwr(Mvec, MLpar, z):
    """
    Double power law with redshift dependence 
    L(M)/Lsun = A * 10^(b1*z) * (M/1e8 Msun)^b2 * (1+M/M_*)^b3
    
    Parameters:
    A	      Overall amplitude, dimensionless
    b1        Redshift slope, dimensionless
    b2        Low mass power law, dimensionless
    b3        High mass power law, dimensionless
    Mstar     Power law turnover mass, in M_sun
    
    >>> Mvec = np.array([1e10,1e11,1e12]) * u.Msun
    >>> MLpar = {'A':5.8e-3, 'b1':0.35, 'b2':1.97, 'b3':-2.92, \
        'Mstar':8.e11*u.Msun}
    >>> z = 3.0
    >>> print DblPwr(Mvec,MLpar,z)
    [    546.6...   37502...  462439...] solLum
    """
    
    A = MLpar['A']
    b1 = MLpar['b1']
    b2 = MLpar['b2']
    b3 = MLpar['b3']
    Mstar = MLpar['Mstar']
    
    L = A * 10.**(b1*z) * (Mvec/(1.e8*u.Msun))**b2 * (1.+Mvec/Mstar)**b3
    L = L*u.Lsun
    
    return L
   


def COMAP_Fid(Mvec,MLpar,z):
    '''
    New COMAP fiducial double-power law model
    L(M) = C/((M/Ms)^A+(M/Ms)^B)
    
    Parameters:
    A   Low-mass slope
    B   High-mass slope
    C   Overall normalization
    Ms  Turnover mass
    
    Predicted values for (A, B, log C, log (Ms/Msol), sigma):
    * pessimistic: (-3.7, 7.0, 11.1, 12.5, 0.36)
    * realistic: (-2.75, 0.05, 10.61, 12.3, 0.42)
    * realistic-plus: (-2.85, -0.42, 10.63, 12.3, 0.42)
    * optimistic: (-2.4, -0.5, 10.45, 12.21, 0.36))
    '''
   
    # Set LIM model parameters for the COMAP fiducial cosmology, the values are
    # taken as the realistic+ prediction, shown as data-driven prior “UM+COLDz+COPSS”
    # in Table 1 of 2111.05931. Model parameters are gone over in more detail in
    # appendix A, the actual L_{CO} function is eq. 8

    A = MLpar['A']
    B = MLpar['B']
    C = MLpar['C']
    Ms = MLpar['Ms']*u.Msun
    
    L = C/((Mvec/Ms)**A+(Mvec/Ms)**B)
    return L*4.9e-5*u.Lsun
   

    
def FIRE_Forged_CO(Mvec,MLpar,z):
    '''
    Model from Tolgay et al 2026 using FIRE galaxies to estimate how L_CO 
    changes with halo mass and redshift.
    '''
    # default interpollation table for line luminosity model parameters
    interp_table_file = ('/home/njcarlson/Repositories/CITA_LIM/limlam_mo'+
            'cker/tables/FIRE-forged_CO.txt')
    # interp_table = np.array([
    #     # z,   A   ,   B   , log C, log(M*/Msun), sigma
    #     [ 0, -3.731, -1.249, 7.178,    11.037   , 0.74 ],
    #     [ 1, -3.357, -1.125, 7.957,    11.643   , 0.80 ],
    #     [ 2, -3.442, -2.119, 7.475,    11.558   , 0.67 ],
    #     [ 3, -3.175, -3.175, 7.205,    11.727   , 0.54 ] ])

    # Load interpollation table file from model parameters if present
    if 'interp_table_file' in MLpar:
        if not MLpar['interp_table_file'] is None:
            interp_table_file = MLpar['interp_table_file']

    # Load the interpollation file
    interp_table = np.genfromtxt( interp_table_file, delimiter=',' )

    # For now, we do a linear interpollation between model parameters 
    A     = np.interp( z, interp_table[:,0], interp_table[:,1] )
    B     = np.interp( z, interp_table[:,0], interp_table[:,2] )
    logC  = np.interp( z, interp_table[:,0], interp_table[:,3] )
    logM  = np.interp( z, interp_table[:,0], interp_table[:,4] )
    Ms    = 10**logM * u.Msun
    sigma = np.interp( z, interp_table[:,0], interp_table[:,5] )

    # Need to implement scatter

    # Luminosity in K km/s pc^2
    L_CO_prime = C / ( ( Mvec/Ms )**A + ( Mvec/Ms )**B )

    # Compute and return Luminosity in solar luminosities
    L_CO = L_CO_prime * 4.9e-5 * u.Lsun
    return L_CO



def TonyLi(Mvec, MLpar, z):
    '''
    CO emission model from Li et al. (2016).  Uses Behroozi et al. SFR(M)
    results.
    
    NOTE ON THIS MODEL: The Li et al. model has two types of scatter: one on
    SFR(M) and one on LCO(SFR), denoted as sigma_SFR and sigma_LCO.  The LCO
    scatter should be entered into LineModel() as the usual sigma_scatter
    input.  However, the SFR scatter behaves differently in that it does not
    preserve mean(LCO), but preserves mean(SFR) instead.  Thus it should be
    given as part of MLpar, there are specific hacks added to LineModel() to
    account for this.
    
    Parameters:
    alpha         Slope of logLIR/logLCO relation, dimensionless
    beta          Intercept of logLIR/logLCO relation, dimensionless
    dMF           10^10 times SFR/LIR normalization (See Li et al. Eq 1), 
                    dimensionless
    BehrooziFile  Filename where Behroozi et al. data is stored, default
                    'sfr_release.dat'. File can be downloaded from
                    peterbehroozi.com/data, (string)
    Mcut_min  Minimum mass below which L=0 (in M_sun)
    Mcut_max  Maximum mass above which L=0 (in M_sun)
    
    >>> Mvec = np.array([1e10,1e11,1e12]) * u.Msun
    >>> MLpar = {'alpha':1.17, 'beta':0.21, 'dMF':1.0,\
        'BehrooziFile':'sfr_release.dat'}
    >>> z = 3.0
    >>> print TonyLi(Mvec,MLpar,z)
    [  2.05...e+02   7.86...e+03   4.56...e+05] solLum
    '''
    
    alpha = MLpar['alpha']
    beta = MLpar['beta']
    dMF = MLpar['dMF']
    BehrooziFile = MLpar['BehrooziFile']
    
    # Read and interpolate Behroozi SFR(M) data
    x = np.loadtxt(BehrooziFile)
    zb = np.unique(x[:,0])-1.
    logMb = np.unique(x[:,1])
    logSFRb = x[:,2].reshape(137,122,order='F')
    
    logSFR_interp = interp2d(logMb,zb,logSFRb,bounds_error=False,fill_value=0.)
    
    # Compute SFR(M) in Msun/yr
    logM = np.log10((Mvec.to(u.Msun)).value)
    if np.array(z).size>1:
        SFR = np.zeros(logM.size)
        for ii in range(0,logM.size):
            SFR[ii] = 10.**logSFR_interp(logM[ii],z[ii])
    else:
        SFR = 10.**logSFR_interp(logM,z)
    
    # Compute IR luminosity in Lsun
    LIR = SFR/(dMF*1e-10)
    
    # Compute L'_CO in K km/s pc^2
    Lprime = (10.**-beta * LIR)**(1./alpha)
    
    # Compute LCO
    L = (4.9e-5*u.Lsun)*Lprime

    return L
    
def SilvaCII(Mvec, MLpar, z):
    '''
    Silva et al. (2015) CII model, relates CII luminosity and SFR by
    log10(L_CII/Lsun) = a_LCII*log10(SFR/(Msun/yr)) + b_LCII
    
    SFR(M) computed from the double power law fit in their Eq. (8), with
    parameters interpolated from their Table 2.
    
    Note that the L(M) relations derived from this model are a variant on the
    DblPwr model above, but with the input parameters changed to match the
    Silva et al. numbers
    
    Parameters:
    a   a_LCII parameter in L(SFR), dimensionless
    b   b_LCII parameter in L(SFR)
    
    >>> Mvec = np.array([1e10,1e11,1e12]) * u.Msun
    >>> MLpar = {'a':0.8475, 'b':7.2203}
    >>> z = 7.5
    >>> print SilvaCII(Mvec,MLpar,z)
    [  4.58...e+06   1.61...e+08   6.89...e+08] solLum
    '''
    
    aLCII = MLpar['a']
    bLCII = MLpar['b']
    
    # Interpolate SFR from Table 2 of Silva et al. 2015
    #SFR = Silva_SFR(Mvec,z)
    SFR = Behroozi_SFR('sfr_release.dat',Mvec,z)*u.Msun/u.yr
    
    # LCII relation
    L = 10**(aLCII*np.log10(SFR/(1*u.Msun/u.yr))+bLCII)*u.Lsun
    
    return L

def LichenCII(Mvec, MLpar, z):

    M0 = MLpar['M0'] * u.Msun
    Mmin = MLpar['Mmin'] * u.Msun
    alpha_MH1 = MLpar['alpha_MH1']
    alpha_LCII = MLpar['alpha_LCII']

    def MH1_fit(M, M_0, M_min, alphaMH1):
        x = M/M_min
        return M_0 * ((x)**alphaMH1) * np.exp(-1/((x)**0.35))

    MH1 = MH1_fit(Mvec, M0, Mmin, alpha_MH1)

    L_CII = alpha_LCII*(MH1.value) * u.Lsun

    '''
    from limlam_mocker import limlam_mocker as llm
    import scipy as sp
    import scipy.interpolate
    sigma_sfr = 0.3


    M_halos = Mvec
    z_z = np.linspace(0, 5, 6)

    tablepath = '/cita/h/home-2/horlaville/clara_limlam/limCode2020-master_clara_2/limlam_mocker/tables/sfr_behroozi_release.dat'

    dat_zp1, dat_logm, dat_logsfr, _ = np.loadtxt(tablepath, unpack=True)

    dat_logzp1 = np.log10(dat_zp1)
    dat_sfr    = 10.**dat_logsfr

    dat_logzp1  = np.unique(dat_logzp1)    # log(z), 1D 
    dat_logm    = np.unique(dat_logm)    # log(Mhalo), 1D        
    dat_sfr     = np.reshape(dat_sfr, (dat_logm.size, dat_logzp1.size))
    
    sfr_interp_tab = sp.interpolate.RectBivariateSpline(dat_logm, dat_logzp1, dat_sfr, kx=1, ky=1)

    z_s = [0 for i in range(len(z_z))]
    for i in range(len(z_z)):
        z_s[i] = np.linspace(z_z[i], z_z[i], len(M_halos))

    
    fit_ps = [0 for i in range(len(z_z))]

    fit_ps[0] = ((4.3*(10**10)), (2*(10**12)), (0.24))
    fit_ps[1] = ((1.5*(10**10)), (6*(10**11)), (0.53)) 
    fit_ps[2] = ((1.3*(10**10)), (3.6*(10**11)), (0.6))
    fit_ps[3] = ((2.9*(10**9)), (6.7*(10**10)), (0.76))
    fit_ps[4] = ((1.4*(10**9)), (2.1*(10**10)), (0.79))
    fit_ps[5] = ((1.9*(10**9)), (2*(10**10)), (0.74))


    MH1_z0 = MH1_fit(M_halos, fit_ps[0]) 
    MH1_z1 = MH1_fit(M_halos, fit_ps[1])
    MH1_z2 = MH1_fit(M_halos, fit_ps[2])
    MH1_z3 = MH1_fit(M_halos, fit_ps[3]) 
    MH1_z4 = MH1_fit(M_halos, fit_ps[4])
    MH1_z5 = MH1_fit(M_halos, fit_ps[5])

    '''

    return L_CII

def LichenCII_v2(Mvec, MLpar, z):

    M0 = MLpar['M0'] * u.Msun
    Mmin = MLpar['Mmin'] * u.Msun
    alpha_MH1 = MLpar['alpha_MH1']
    alpha_LCII = MLpar['alpha_LCII']
   
    # First, load the SFR for each halo following the prescription from the TonyLi model and using the data from the parameter 'BehrooziFile'.

    BehrooziFile = MLpar['BehrooziFile']

    # Read and interpolate Behroozi SFR(M) data
    x = np.loadtxt(BehrooziFile)
    zb = np.unique(x[:,0])-1.
    logMb = np.unique(x[:,1])
    logSFRb = x[:,2].reshape(137,122,order='F')

    logSFR_interp = interp2d(logMb,zb,logSFRb,bounds_error=False,fill_value=0.)

    # Compute SFR(M) in Msun/yr
    logM = np.log10((Mvec.to(u.Msun)).value)
    if np.array(z).size>1:
        SFR = np.zeros(logM.size)
        for ii in range(0,logM.size):
            SFR[ii] = 10.**logSFR_interp(logM[ii],z[ii])
    else:
        SFR = 10.**logSFR_interp(logM,z)

    # Then, we compute the stellar mass for each halo using the prescription by Behroozi+2013 [code snippet from DTC]:

    a = lambda z:1/(1.+z)
    # Behroozi et al 2013a parameters
    nu = lambda z:np.exp(-4*a(z)**2)
    log10_eps = lambda z:-1.777-0.006*(a(z)-1)*nu(z)-0.119*(a(z)-1)
    log10_M1 = lambda z:11.514+(-1.793*(a(z)-1)-0.251*z)*nu(z)
    alpha = lambda z:-1.412+0.731*(a(z)-1)*nu(z)
    delta = lambda z:3.508+(2.608*(a(z)-1)-0.043*z)*nu(z)
    gamma = lambda z:0.316+(1.319*(a(z)-1)+0.279*z)*nu(z)
    f = lambda x,z:-np.log10(10**(alpha(z)*x)+1)+delta(z)*(np.log10(1+np.exp(x)))**gamma(z)/(1+np.exp(10**(-x)))
    xi = lambda z:0.218-0.023*(a(z)-1)

    def stellar_m(halo_m,z,scatter=True):
        sm = 10**(log10_eps(z)+log10_M1(z)+f(np.log10((halo_m/(10**log10_M1(z))).value),z)-f(0,z))
        if scatter:
            rand = np.random.lognormal(-0.5*(xi(z)*np.log(10))**2,xi(z)*np.log(10))
            return sm*rand
        else:
            return sm

    stellar_mass = stellar_m(Mvec, z)

    # Then, retrieve metallicity with FMR from Heintz+2021 (from Curti+2020)

    gamma = 0.31
    beta = 2.1
    m_0 = 10.11
    m_1 = 0.56
    Z_0 = 8.779

    def M_0(sfr):
        return (10**(m_0))*(sfr**(m_1))

    halo_M0 = M_0(SFR)

    def ps_metal(stellar_m, M_0):
        return Z_0 - (gamma/beta)*np.log10(1 + (stellar_m/M_0)**(-beta))

    # The ps_metal is the FMR but for the 'pseudo'-metallicity quantity described in Heintz+21; it's not the direct metallicity. It's defined as 12 + log(O/H).
    # In solar metallicity, log(Z/Z_sol) = 0 when 12 + log(O/H) = 8.69, so 'real'Z = 10**('pseudo'Z - 8.69):

    def metal(ps_m):
        return 10**(ps_m - 8.69)

    halo_psZ = ps_metal(stellar_mass, halo_M0)
    halo_Z = metal(halo_psZ)

    def MH1_fit(M, M_0, M_min, alphaMH1):
        x = M/M_min
        return M_0 * ((x)**alphaMH1) * np.exp(-1/((x)**0.35))

    MH1 = MH1_fit(Mvec, M0, Mmin, alpha_MH1)

    L_CII = alpha_LCII*(MH1.value)*halo_Z * u.Lsun
    
    return L_CII


def LichenCII_v3(Mvec, MLpar, z):

    M0 = MLpar['M0'] * u.Msun
    Mmin = MLpar['Mmin'] * u.Msun
    alpha_MH1 = MLpar['alpha_MH1'] # try different values and look for effects on relative entropy plot
    alpha_LCII = MLpar['alpha_LCII']
    zdex = MLpar['zdex']
    alpha_0 = MLpar['alpha0']
    gamma_0 = MLpar['gamma0']
   
    # First, load the SFR for each halo following the prescription from the TonyLi model and using the data from the parameter 'BehrooziFile'.

    BehrooziFile = MLpar['BehrooziFile']

    # Read and interpolate Behroozi SFR(M) data
    x = np.loadtxt(BehrooziFile)
    zb = np.unique(x[:,0])-1.
    logMb = np.unique(x[:,1])
    logSFRb = x[:,2].reshape(137,122,order='F')

    logSFR_interp = interp2d(logMb,zb,logSFRb,bounds_error=False,fill_value=0.)

    # Compute SFR(M) in Msun/yr
    logM = np.log10((Mvec.to(u.Msun)).value)
    if np.array(z).size>1:
        SFR = np.zeros(logM.size)
        for ii in range(0,logM.size):
            SFR[ii] = 10.**logSFR_interp(logM[ii],z[ii])
    else:
        SFR = 10.**logSFR_interp(logM,z)
   
    # Checking for the effect of SFR variation on the CII luminosities:
    # SFR = MLpar[sfr_frac]*SFR


    # Then, we compute the stellar mass for each halo using the prescription by Behroozi+2013 [code snippet from DTC]:

    a = lambda z:1/(1.+z)
    # Behroozi et al 2013a parameters
    nu = lambda z:np.exp(-4*a(z)**2)
    log10_eps = lambda z:-1.777-0.006*(a(z)-1)*nu(z)-0.119*(a(z)-1)
    log10_M1 = lambda z:11.514+(-1.793*(a(z)-1)-0.251*z)*nu(z)
    alpha = lambda z:alpha_0+0.731*(a(z)-1)*nu(z)
    delta = lambda z:3.508+(2.608*(a(z)-1)-0.043*z)*nu(z)
    gamma = lambda z:0.316+(1.319*(a(z)-1)+0.279*z)*nu(z)
    f = lambda x,z:-np.log10(10**(alpha(z)*x)+1)+delta(z)*(np.log10(1+np.exp(x)))**gamma(z)/(1+np.exp(10**(-x)))
    xi = lambda z:0.218-0.023*(a(z)-1)

    def stellar_m(halo_m,z,scatter=False):
        sm = 10**(log10_eps(z)+log10_M1(z)+f(np.log10((halo_m/(10**log10_M1(z))).value),z)-f(0,z))
        if scatter:
            rand = np.random.lognormal(-0.5*(xi(z)*np.log(10))**2,xi(z)*np.log(10))
            return sm*rand
        else:
            return sm

    stellar_mass = stellar_m(Mvec, z)

    # Then, retrieve metallicity with FMR from Heintz+2021 (from Curti+2020)

    gamma = gamma_0
    beta = 2.1
    m_0 = 10.11
    m_1 = 0.56
    Z_0 = 8.779

    def M_0(sfr):
        return (10**(m_0))*(sfr**(m_1))

    halo_M0 = M_0(SFR)

    def ps_metal(stellar_m, M_0):
        return Z_0 - (gamma/beta)*np.log10(1 + (stellar_m/M_0)**(-beta))

    # The ps_metal is the FMR but for the 'pseudo'-metallicity quantity described in Heintz+21; it's not the direct metallicity. It's defined as 12 + log(O/H).
    # In solar metallicity, log(Z/Z_sol) = 0 when 12 + log(O/H) = 8.779, so 'real'Z = 10**('pseudo'Z - 8.779):

    def metal(ps_m):
        return 10**(ps_m - 8.779)

    halo_psZ = ps_metal(stellar_mass, halo_M0)
    halo_Z = metal(halo_psZ)
    
    Z_scatter = add_log_normal_scatter(halo_Z, zdex, seed = 23)

    def MH1_fit(M, M_0, M_min, alphaMH1):
        x = M/M_min
        return M_0 * ((x)**alphaMH1) * np.exp(-1/((x)**0.35))

    MH1 = MH1_fit(Mvec, M0, Mmin, alpha_MH1)

    L_CII = (alpha_LCII*(MH1.value)*Z_scatter)*u.Lsun
    
    return L_CII



def SilvaCO(Mvec, MLpar, z):
    '''
    Silva et al. (2015) CO model, relates CO luminosity to halo mass via a
    quadruple power law.  Can be called with either the values of the power
    law parameters or with the name of the transition (e.g. 'trans':'3-2'),
    in which case values from Table 4 of Silva et al. will be used.
    
    If MLpar contains both 'trans' argument and 
    
    Doctest uses parameters for CO(3-2).
    
    Parameters:
    L0                Overall normalization of L(M), dimensionless
    M1,M2,M3,M4       Locations of power law turnovers
    d0,d1,d2,d3,d4    Power law slopes
    
    >>> Mvec = np.array([1e10,1e11,1e12]) * u.Msun
    >>> MLpar = {'L0':3.0e-24,'M1':6.0e11*u.Msun,'M2':5.0e12*u.Msun,'M3':\
                 4e13*u.Msun,'M4':1.0*u.Msun,'d0':2.6,'d1':-3.5,'d2':0.2,\
                 'd3':2.2,'d4':0.0}
    >>> z = 0.5
    >>> print SilvaCO(Mvec,MLpar,z)
    [  2.83...e+02   7.02...e+04   1.68...e+06] solLum
    >>> MLpar_1 = {'trans':'3-2'}
    >>> print SilvaCO(Mvec,MLpar_1,z)
    [  2.83...e+02   7.02...e+04   1.68...e+06] solLum
    '''
    
    if 'trans' in MLpar:
        if MLpar['trans']=='2-1':
            par = dict(L0=4.7e-29,M1=1.0e11*u.Msun,M2=6.0e11*u.Msun,
                       M3=5.0e12*u.Msun,M4=5.0e14*u.Msun,d0=3.05,d1=-2.0,
                       d3=1.9,d4=5.0)
        elif MLpar['trans']=='3-2':
            par = dict(L0=3.0e-24,M1=6.0e11*u.Msun,M2=5.0e12*u.Msun,
                       M3=4.0e13*u.Msun,M4=1.0*u.Msun,d0=2.6,d1=-3.5,d2=0.2,
                       d3=2.2,d4=0.0)
        elif MLpar['trans']=='4-3':
            par = dict(L0=8.0e-18,M1=9.0e11*u.Msun,M2=5.0e12*u.Msun,
                       M3=3.0e13*u.Msun,M4=1.0*u.Msun,d0=2.05,d1=-1.7,d2=-1.8,
                       d3=2.3,d4=0.0)
        elif MLpar['trans']=='5-4':
            par = dict(L0=3.5e-18,M1=2.0e12*u.Msun,M2=4.0e12*u.Msun,
                       M3=1.0e13*u.Msun,M4=1.0*u.Msun,d0=2.05,d1=-2.0,d2=-3.9,
                       d3=4.5,d4=0.0)
        elif MLpar['trans']=='6-5':
            #par = dict(L0=4.0e-18,M1=9.0e10*u.Msun,M2=6.0e11*u.Msun,
            par = dict(L0=8.0e-18,M1=9.0e11*u.Msun,M2=5.0e12*u.Msun,
                       M3=3.0e13*u.Msun,M4=1.0*u.Msun,d0=2.05,d1=-1.7,d2=-1.8,
                       d3=2.3,d4=0.0)
        elif MLpar['trans']=='5-4':
            par = dict(L0=3.5e-18,M1=2.0e12*u.Msun,M2=4.0e12*u.Msun,
                       M3=1.0e13*u.Msun,M4=1.0*u.Msun,d0=2.05,d1=-2.0,d2=-3.9,
                       d3=4.5,d4=0.0)
        elif MLpar['trans']=='6-5':
            par = dict(L0=4.0e-18,M1=9.0e10*u.Msun,M2=6.0e11*u.Msun,
                       M3=4.0e13*u.Msun,M4=1.0*u.Msun,d0=2.6,d1=-3.5,d2=0.2,
                       d3=2.2,d4=0.0)
        elif MLpar['trans']=='4-3':
            par = dict(L0=8.0e-18,M1=9.0e11*u.Msun,M2=5.0e12*u.Msun,
                       M3=3.0e13*u.Msun,M4=1.0*u.Msun,d0=2.05,d1=-1.7,d2=-1.8,
                       d3=2.3,d4=0.0)
        elif MLpar['trans']=='5-4':
            par = dict(L0=3.5e-18,M1=2.0e12*u.Msun,M2=4.0e12*u.Msun,
                       M3=1.0e13*u.Msun,M4=1.0*u.Msun,d0=2.05,d1=-2.0,d2=-3.9,
                       d3=4.5,d4=0.0)
        elif MLpar['trans']=='6-5':
            par = dict(L0=4.0e-18,M1=9.0e10*u.Msun,M2=6.0e11*u.Msun,
                       M3=1.0*u.Msun,M4=1.0*u.Msun,d0=2.0,d1=1.5,d2=-3.75,
                       d3=0.0,d4=0.0)
        else:
            raise ValueError("Invalid transition")
    else:
        par = MLpar
    
    
    L0 = par['L0']
    M1 = par['M1']
    M2 = par['M2']
    M3 = par['M3']
    M4 = par['M4']
    d0 = par['d0']
    d1 = par['d1']
    d2 = par['d2']
    d3 = par['d3']
    d4 = par['d4']
    
    L = (L0 * Mvec**d0 * (1+Mvec/M1)**d1 * (1+Mvec/M2)**d2 * (1+Mvec/M3)**d3 *
         (1+Mvec/M4)**d4).value*u.Lsun
         
    return L
    
def MHI_21cm(Mvec, MLpar, z):
    '''
    Obuljen et al. (2018) 21cm MHI(M) model, relates MHI to halo mass by
    MHI = M0 * (M/Mmin)^alpha * exp(-Mmin/M)
    
    NOTE that the best fit values given by Obuljen et al. for M0 and Mmin are
    in Msun/h units
    
    Parameters
    M0      Overall normalization of MHI(M) (in Msun)
    Mmin    Location of low-mass exponential cutoff (in Msun)
    alpha   Slope at high-mass (dimensionless)
    
    >>> Mvec = np.array([1e10,1e11,1e12]) * u.Msun
    >>> MLpar = {'M0':4.73e8*u.Msun,'Mmin':2.66e11*u.Msun,'alpha':0.44}
    >>> z = 0.03
    >>> print MHI_21cm(Mvec,MLpar,z)
    [  1.94...e-12   1.33...e-01   4.03...e+00] solLum
    '''
    
    M0 = MLpar['M0']
    Mmin = MLpar['Mmin']
    alpha = MLpar['alpha']
    
    CLM = 6.215e-9*u.Lsun/u.Msun # Conversion factor btw MHI and LHI
    
    MHI = M0*(Mvec/Mmin)**alpha*np.exp(-Mmin/Mvec)
    L = CLM*MHI
    return L
    
def Constant_L(Mvec, MLpar, z):
    '''
    Model where every halo has a constant luminosity independent of mass.
    Still has cutoffs at Mcut_min and Mcut_max.
    
    Intended primarily for sanity checks and debugging.
    
    Parameters:
    L0  Luminosity of every halo
    
    >>> Mvec = np.array([1e10,1e11,1e12]) * u.Msun
    >>> MLpar = {'L0':1*u.Lsun}
    >>> z = 1
    >>> print Constant_L(Mvec,MLpar,z)
    [ 1.  1.  1.] solLum
    '''
    
    L0 = MLpar['L0']
    
    return L0*np.ones(Mvec.size)

def Padmanabhan_CII(Mvec, MLpar, z):
    '''
    CII Model from Padmanabhan (2018)

    L_CII = F(z) * (M/M1)^beta * exp(-N1/M)
    F(z) = ((1+z)^2.7/(1+[(1+z)/2.9]^5.6))^alpha

    Parameters:
    M1     Overall normalization (Msun)
    beta   Slope of low-mass power law
    N1     Location of low-mass cutoff (Msun)
    alpha  Slope of redshift scaling   

    >>> Mvec = np.array([1e10,1e11,1e12]) * u.Msun
    >>> MLpar = dict(M1=2.39e-5*u.Msun,N1=4.19e11*u.Msun,alpha=1.79,beta=0.49)
    >>> z = 3.
    >>> print Padmanabhan_CII(Mvec,MLpar,z)
    42
    '''

    M1 = MLpar['M1']
    beta = MLpar['beta']
    N1 = MLpar['N1']
    alpha = MLpar['alpha']

    F = ((1+z)**2.7/(1+((1+z)/2.9)**5.6))**alpha

    return F*(Mvec/M1)**beta*np.exp(-N1/Mvec)*u.Lsun

def Be13_Pwr(Mvec, MLpar, z):
    '''
    Assumes a power law relation between L and SFR, and SFR(M) from Behroozi
    et al. 2013
    
    Parameters:
    K       Normalization between L in Lsun and SFR^gamma
    gamma   Power law slope between L and SFR
    '''
    
    K = MLpar['K']
    gamma = MLpar['gamma']
    
    # Read and interpolate Behroozi SFR(M) data
    BehrooziFile = 'sfr_release.dat'
    x = np.loadtxt(BehrooziFile)
    zb = np.unique(x[:,0])-1.
    logMb = np.unique(x[:,1])
    logSFRb = x[:,2].reshape(137,122,order='F')
    
    logSFR_interp = interp2d(logMb,zb,logSFRb,bounds_error=False,fill_value=0.)
    
    # Compute SFR(M) in Msun/yr
    logM = np.log10((Mvec.to(u.Msun)).value)
    if np.array(z).size>1:
        SFR = np.zeros(logM.size)
        for ii in range(0,logM.size):
            SFR[ii] = 10.**logSFR_interp(logM[ii],z[ii])
    else:
        SFR = 10.**logSFR_interp(logM,z)
        
    return (K*SFR**gamma)*u.Lsun

    
def CII_Metallicity(Mvec,MLpar,z):
    '''
    Metallicity-based collisional excitation model from Gong 2012 and
    Silva 2015
    
    Parameters:
    Te  Electron kinetic temperature, in K
    ne  Electron number density, in cm**-3
    '''
    
    Te = MLpar['Te']
    ne = MLpar['ne']
    
    lCII = 158*u.um
    nuCII = (cu.c/lCII).to(u.GHz)
    
    TCMB = 2.7*u.K
    xul = 1./(np.exp(cu.h*nuCII/(cu.k_B*TCMB))-1)
    
    Tei = np.array([100,1000,10000])
    Cul_i = np.array([3.41e-7,1.09e-7,4.55e-8])

    Cul = interp1d(Tei,Cul_i)((Te/u.K).decompose().value)*u.cm**3/u.s
    
    Tstar = 91*u.K
    Aul = 2.36e-6/u.s
    
    Tstar_TS = np.log((Aul*(1+xul)+ne*Cul)/(Aul*xul+ne*Cul*np.exp(-Tstar/Te)))
    
    #M0 = 10**9*u.Msun
    #Mc = 3.5e11*u.Msun
    #MZ = M0*(Mvec/Mc)**3.6*(1+Mvec/Mc)**-3.25
    
    M0 = (z-1)*u.Msun
    Ma = 10**8*u.Msun
    Mb = 9e9*u.Msun
    Mc = 2e12*u.Msun
    Md = 2e13*u.Msun
    a = 1.7
    b = 1.0
    c = -5
    d = 2.5
    MZ = M0*(Mvec/Ma)**a*(1+Mvec/Mb)**b*(1+Mvec/Mc)**c*(1+Mvec/Md)**d
    
    mc = 12*cu.m_p
    
    return (2*cu.h*nuCII*Aul*np.exp(-Tstar_TS)*MZ/(3*(1+z)**3*mc)).to(u.Lsun)
    #return (2*cu.h*nuCII*Aul*MZ/(3*(1+z)**3*mc)).to(u.Lsun)
    
def Chung_Lya(Mvec,MLpar,z):
    '''
    Lya model for CO-HETDEX cross-correlation from Chung+2018
    '''
    A = MLpar['A']
    BehrooziFile = MLpar['BehrooziFile']
#    z0 = MLpar['z0'] 
    
    SFR = Behroozi_SFR(BehrooziFile,Mvec,z)
    
    fesc = (0.18+0.82/(1+0.8*SFR**0.875))**2/np.sqrt(1+np.exp(-1.6*z+5))
    
    return (A*SFR*fesc*u.erg/u.s).to(u.Lsun)
    
def Chung_Lya_Counts(Mvec,MLpar,z):
    '''
    Returns counts of Lya emitters above detection threshold using Chung_Lya
    model
    '''
    p = dict(A=MLpar['A'],BehrooziFile=MLpar['BehrooziFile'])
    Lmin = MLpar['Lmin']
    XLT = MLpar['XLT']
    z0 = MLpar['z0']
    
    L0 = Chung_Lya(Mvec,MLpar,z)
    ct = np.zeros(Mvec.shape)*u.Lsun
    
    ct[L0>=Lmin] = 1*u.Lsun
    
    
    return ct
    
    
def HamsaTest(Mvec,MLpar,z):
    '''
    Test CII model from Hamsa
    '''
    mhalo=(Mvec.to(u.Msun)).value
    M1 = 2.39e-5
    N1 = 4.19e11
    beta = 0.49
    alpha = 1.79
    term1 = (1+z)**2.7/(1 + ((1+z)/2.9)**5.6) # redshift evolution of SFR
    mhiab = ((term1)**(alpha))*((mhalo/M1)**(beta))*np.exp(-N1/mhalo)
    return mhiab*u.Lsun
    
def Simon_file(Mvec,MLpar,z):
    '''
    Reads mass-luminosity model from a file in the format Simon sent me
    '''
    fname = MLpar['filename']
    h = MLpar['hubble']
    
    x = np.loadtxt(fname)
    
    zf = np.unique(x[:,0])
    Mf = np.unique(x[:,1]) # Table halo masses in Msun/h
    Lf = x[:,2].reshape(zf.size,Mf.size) # Table lum's in erg/s

    L_interp = interp2d(Mf,zf,Lf,bounds_error=False,fill_value=0)
    
    Msunh = u.Msun/h
    
    return (L_interp((Mvec.to(Msunh)).value,z)*u.erg/u.s).to(u.Lsun)
    L0 = interp1d(z0,np.array([4.3,6.2,4.0])*1e6)(z)*u.Lsun
    b = interp1d(z0,np.array([2.4,2.6,2.8]))(z)
    Mc = interp1d(z0,np.array([3.5,3.0,2.0])*1e11)(z)*u.Msun
    d = interp1d(z0,np.array([2.8,3.4,3.3]))(z)
    
    return R*L0*(Mvec/Mc)**b*(1+Mvec/Mc)**-d
    
    
    
###################
# Other functions #
###################
def Behroozi_SFR(BehrooziFile, M, z):
    '''
    Returns SFR(M,z) interpolated from Behroozi et al.
    '''
    x = np.loadtxt(BehrooziFile)
    zb = np.unique(x[:,0])-1.
    logMb = np.unique(x[:,1])
    logSFRb = x[:,2].reshape(137,122,order='F')
    
    logSFR_interp = interp2d(logMb,zb,logSFRb,bounds_error=False,fill_value=0.)
    
    logM = np.log10((M.to(u.Msun)).value)
    if np.array(z).size>1:
        SFR = np.zeros(logM.size)
        for ii in range(0,logM.size):
            SFR[ii] = 10.**logSFR_interp(logM[ii],z[ii])
    else:
        SFR = 10.**logSFR_interp(logM,z)
    
    return SFR
    

def Silva_SFR(M,z):
    '''
    Returns SFR(M,z) interpolated from values in Table 2 of Silva et al. 2015
    '''
    x = np.loadtxt('Silva15_SFR_params.dat')
    z0 = x[0,:]
    M0 = x[1,:]*u.Msun/u.yr
    Ma = x[2,:]*u.Msun
    Mb = x[3,:]*u.Msun
    a = x[4,:]
    b = x[5,:]
    SFR0 = np.zeros((z0.size,M.size))*u.Msun/u.yr
    for ii in range(0,z0.size):
        SFR0[ii,:] = M0[ii]*(M/Ma[ii])**a[ii]*(1+M/Mb[ii])**b[ii]
    return interp1d(z0,SFR0.value,axis=0)(z)*u.Msun/u.yr

    """
    z0 = x[0,:]
    M0 = interp1d(z0,x[1,:])(z)*u.Msun/u.yr
    Ma = interp1d(z0,x[2,:])(z)*u.Msun
    Mb = interp1d(z0,x[3,:])(z)*u.Msun
    a = interp1d(z0,x[4,:])(z)
    b = interp1d(z0,x[5,:])(z)
    
    return M0*(M/Ma)**a*(1+M/Mb)**b
    """

if __name__ == "__main__":
    import doctest

    doctest.testmod(optionflags=doctest.ELLIPSIS |
                    doctest.NORMALIZE_WHITESPACE)
