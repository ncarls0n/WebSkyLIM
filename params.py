import astropy.units as u

default_par = dict(cosmo_model = 'Planck15',   # Input cosmological model
              model_type = 'LF',    # Type of line model, either 'LF' or 'ML'
              model_par = {'phistar':8.7e-11*u.Lsun**-1*u.Mpc**-3,
                           'Lstar':2.1e6*u.Lsun,
                           'alpha':-1.87,
                           'Lmin':5000*u.Lsun},    # Model parameters
              nu = 115*u.GHz,    # Rest frequency of target line
              nuObs = 30*u.GHz,    # Central frequency of observation
              Mmin = 1e9*u.Msun,    # Minimum mass for line emission
              Mmax = 1e15*u.Msun,    # Maximum mass for line emission
              nM = 5000,    # Number of points to compute mass functions
              hmf_model = 'Tinker08',    # Mass function model. Choose from hmf
              Lmin = 10*u.Lsun,    # Minimum luminosity for dn/dL
              Lmax = 1e8*u.Lsun,    # Maximum luminosity for dn/dL
              nL = 5000,    # Number of points to compute luminosity functions
              kmin = 1e-2/u.Mpc,    # Minimum wavenumber for power spectra
              kmax = 7/u.Mpc,    # Maximum wavenumber for power spectra
              nk = 100,    # Number of wavenumber points for power spectra
              sigma_scatter = 0.,    # Scatter in mass-luminosity relation
              fduty = 1.,    # Fraction of halos emitting at any given time
              do_onehalo = False,    # Bool to include one-halo contributions
              do_Jysr = False,    # Compute quantities in brightness temp if
                                 # False, Jy/sr if True
                                 
              # PARAMETERS USED FOR Line_Obs MODELS
              Tsys_NEFD = 40*u.K,    # System temperature
              Nfeeds = 19,    # Number of detectors
              beam_FWHM = 4.1*u.arcmin,     # Beam full width at half max
              Delta_nu = 8*u.GHz,    # Total frequency bandwidth
              dnu = 15.6*u.MHz,    # Width of single frequency channel
              tobs = 6000*u.hr,    # Observing time per field
              Omega_field = 2.25*u.deg**2,     # Solid angle of single field
              Nfield = 1,    # Number of fields to be observed
              Tmin_VID = 1e-3*u.uK,    # Minimum intensity for VID
              Tmax_VID = 1000*u.uK,    # Maximum intensity for VID
              nT = 10**5,    # Number of intensity points for VID
              do_fast_VID = True, # Bool, do VID with FFT's if True, brute-
                                 # force integration if False
              sigma_G = 1.6,    # Gaussian variance parameter for VID
              Ngal_max = 100,    # Max sources/voxel for VID integrals
              Nbin_hist = 101,   # Number of bins to predict voxel counts
              subtract_VID_mean = False,    # Bool to use absolute or relative
                                           # intensity for VID
              linear_VID_bin = False,    # Output predicted histogram in linear
                                        # or logarithmically spaced bins
              
              # PARAMETERS USED FOR LIMLAM SIMULATIONS
              catalogue_file = 
                  'limlam_mocker/catalogues/default_catalogue.npz',
                  # Location of peak-patch catalog file
              map_output_file = 'limlam_output.npz' # Output file location
              )

# Tony Li model and COMAP1 CURRENTLY WITHOUT SCATTER
TonyLi_PhI = dict(cosmo_model = 'Planck15',
                  model_type = 'ML',
                  model_name = 'TonyLi',
                  model_par = {'alpha':1.17,'beta':0.21,'dMF':1.0,
                               'BehrooziFile':'sfr_release.dat','sig_SFR':0.0},
                  nu = 115.27*u.GHz,
                  nuObs = 15*u.GHz,
                  Mmin = 1e9*u.Msun, # Set high to compare to incomplete sims
                  Mmax = 1e15*u.Msun,                    
                  nM = 5000,
                  nL = 5000,
                  kmin = 1e-2/u.Mpc,
                  kmax = 7./u.Mpc,
                  nk = 100,
                  sigma_scatter = 0.0,
                  Tsys_NEFD = 44*u.K,
                  Nfeeds = 19,
                  beam_FWHM = 0.8*u.arcmin,
                  Delta_nu = 4*u.GHz,
                  dnu = 20*u.MHz,
                  tobs = 1500*u.hr,
                  Omega_field = 4*u.deg**2,
                  do_Jysr = False,
                  Nfield = 4,
                  catalogue_file = '/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz',
                  )

TonyLi_PhII = dict(cosmo_model = 'Planck15',
                  model_type = 'ML',
                  model_name = 'TonyLi',
                  model_par = {'alpha':1.17,'beta':0.21,'dMF':1.0,
                               'BehrooziFile':'sfr_release.dat','sig_SFR':0.0},
                  nu = 115.27*u.GHz,
                  nuObs = 26*u.GHz,
                  Mmin = 1e9*u.Msun, # Set high to compare to incomplete sims
                  Mmax = 1e15*u.Msun,
                  nM = 5000,
                  nL = 5000,
                  kmin = 1e-2/u.Mpc,
                  kmax = 7./u.Mpc,
                  nk = 100,
                  sigma_scatter = 0.0,
                  Tsys_NEFD = 44*u.K,
                  Nfeeds = 19,
                  beam_FWHM = 0.8*u.arcmin,
                  Delta_nu = 13.7*u.GHz,
                  dnu = 300*u.MHz,
                  tobs = 1500*u.hr,
                  Omega_field = 538*u.deg**2,
                  do_Jysr = False,
                  Nfield = 4,
                  catalogue_file = '/mnt/scratch-lustre/horlaville/nate_sims/may_5/ngaussian/ngaussian_lc_mlambda10.npz',
                  )

# COMAP Fiducial
COMAP_Fid = dict(cosmo_model = 'Planck15',
                  model_type = 'ML',
                  model_name = 'COMAP_Fid',

                  # COMAP_Fid model parameters using 
                  model_par = {'A':-2.85,'B':-0.42,'C':10**10.63,'Ms':10**12.3*u.Msun},
                  # Set LIM model parameters for the COMAP fiducial cosmology, the values are
                  # taken as the realistic+ prediction, shown as data-driven prior “UM+COLDz+COPSS”
                  # in Table 1 of 2111.05931. Model parameters are gone over in more detail in
                  # appendix A, the actual L_{CO} function is eq. 8

                  # Observatory parameters
                  # CO(1-0) rest frame frequency
                  nu    = 115.27*u.GHz,

                  # CO(2-1) rest frame freqeuncy
                  # nu = 230.5380*u.GHz,

                  # CO(3-2) rest frame frequency
                  # nu = 345.7960*u.GHz,

                  # Ku band frequency
                  nuObs = 15*u.GHz,

                  # Ka band frequency
                  # nuObs = 30*u.GHz,

                  # Q band frequency
                  # nuObs = 45*u.GHz,

                  # Dark matter halo masses
                  Mmin = 1e9*u.Msun, # Set high to compare to incomplete sims
                  Mmax = 1e15*u.Msun,
                  nM = 5000,
                  nL = 5000,
                  kmin = 1e-2/u.Mpc,
                  kmax = 7./u.Mpc,
                  nk = 100,
                  sigma_scatter = 0.0,
                  Tsys_NEFD = 44*u.K,
                  Nfeeds = 19,
                  beam_FWHM = 0.8*u.arcmin,
                  Delta_nu = 13.7*u.GHz,
                  dnu = 300*u.MHz,
                  tobs = 1500*u.hr,
                  Omega_field = 538*u.deg**2,
                  do_Jysr = False,
                  Nfield = 4,
                  catalogue_file = '/mnt/scratch-lustre/horlaville/nate_sims/may_5/ngaussian/ngaussian_lc_mlambda10.npz',
                  )




# /home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz
# Gaussian lightcone: nate_sims/gaussian/gauss_lightcone.npz
# Non-Gaussian Instability lightcone: nate_sims/ngauss_instab/ngauss_instab_lightcone.npz
# Non-Gaussian Power Spectrum lightncone: nate_sims/ngauss_pspec/ngauss_pspec_lightcone.npz

# Silva et al. (2015) CII model m2 with CCATp parameters at z~6
Silva_m1_z6_CCATp = dict(model_type = 'ML',
                         model_name = 'SilvaCII',
                         model_par = {'a':0.8475,'b':7.2203},
                         nu = 1897*u.GHz,
                         nuObs = 250*u.GHz,
                         Mmin = 1e10*u.Msun,
                         Mmax = 1e14*u.Msun,
                         Tsys_NEFD=48.5*u.mJy*u.s**(1./2),
                         Nfeeds=1.7,
                         beam_FWHM=46*u.arcsec,
                         Delta_nu=60*u.GHz, # was 20 before, why?
                         dnu=2.7*u.GHz,
                         tobs=4000*u.hr,
                         Omega_field=2*u.deg**2,
                         do_Jysr=True,
                         Nfield = 4,
                         catalogue_file = '/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz'
                         )
Silva_m2_z6_CCATp = dict(model_type = 'ML',
                         model_name = 'SilvaCII',
                         model_par = {'a':1.0,'b':6.9647},
                         nu = 1897*u.GHz,
                         nuObs = 250*u.GHz,
                         Mmin = 1e10*u.Msun,
                         Mmax = 1e14*u.Msun,
                         Tsys_NEFD=48.5*u.mJy*u.s**(1./2),
                         Nfeeds=1.7,
                         beam_FWHM=46*u.arcsec,
                         Delta_nu=60*u.GHz, # was 20 before, why?
                         dnu=2.7*u.GHz,
                         tobs=4000*u.hr,
                         Omega_field=2*u.deg**2,
                         do_Jysr=True,
                         Nfield = 4,
                         catalogue_file = '/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz')

Lichen = dict(model_type = 'ML',
                         model_name = 'LichenCII',
                         model_par = {'M0':(1.9*(10**9)),'Mmin':(2*(10**10)), 'alpha_MH1':0.74, 'alpha_LCII':0.005},
                         nu = 1897*u.GHz,
                         nuObs = 250*u.GHz,
                         Mmin = 1e10*u.Msun,
                         Mmax = 1e14*u.Msun,
                         Tsys_NEFD = 72.5*u.mJy*u.s**(1./2),
                         Nfeeds = 120,
                         beam_FWHM = 48*u.arcsec,
                         Delta_nu = 60*u.GHz, # was 20 before, why?
                         dnu = 2.8*u.GHz,
                         tobs = 2000*u.hr,
                         Omega_field = 4*u.deg**2,
                         do_Jysr = True,
                         Nfield = 4,
                         catalogue_file = '/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz'
                         )

Lichen_v2 = dict(model_type = 'ML',
                         model_name = 'LichenCII_v2',
                         model_par = {'M0':(1.9*(10**9)),'Mmin':(2*(10**10)), 'alpha_MH1':0.74, 'alpha_LCII':0.017, 'BehrooziFile': 'sfr_reinterp.dat'},
                         nu = 1897*u.GHz,
                         nuObs = 250*u.GHz,
                         Mmin = 1e10*u.Msun,
                         Mmax = 1e14*u.Msun,
                         Tsys_NEFD = 72.5*u.mJy*u.s**(1./2),
                         Nfeeds = 120,
                         beam_FWHM = 48*u.arcsec,
                         Delta_nu = 60*u.GHz, # was 20 before, why?
                         dnu = 2.8*u.GHz,
                         tobs = 2000*u.hr,
                         Omega_field = 4*u.deg**2,
                         do_Jysr = True,
                         Nfield = 4,
                         catalogue_file = '/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz'
                         )

Lichen_v3 = dict(model_type = 'ML',
                         model_name = 'LichenCII_v3',
                         model_par = {'zdex': 0.05, 'M0':(1.9*(10**9)),'Mmin':(2*(10**10)), 'alpha_MH1':0.74, 'alpha_LCII':0.017, 'BehrooziFile': 'sfr_reinterp.dat'},
                         nu = 1897*u.GHz,
                         nuObs = 270*u.GHz,
                         Mmin = 1e10*u.Msun,
                         Mmax = 1e14*u.Msun,
                         Tsys_NEFD = 72.5*u.mJy*u.s**(1./2),
                         Nfeeds = 120,
                         beam_FWHM = 48*u.arcsec,
                         Delta_nu = 40*u.GHz, # was 20 before, why?
                         dnu = 1.11*u.GHz,
                         tobs = 2000*u.hr,
                         Omega_field = 4*u.deg**2,
                         do_Jysr = True,
                         Nfield = 4,
                         catalogue_file = '/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz'
                         )

Lichen_v4 = dict(model_type = 'ML',
                         model_name = 'LichenCII_v3',
                         model_par = {'zdex': 0.05, 'M0':(1.9*(10**9)),'Mmin':(2*(10**10)), 'alpha_MH1':0.74, 'alpha_LCII':0.017, 'alpha0':-1.412, 'gamma0':0.31, 'BehrooziFile': 'sfr_reinterp.dat'},
                         nu = 1897*u.GHz,
                         nuObs = 400*u.GHz,
                         Mmin = 1e10*u.Msun,
                         Mmax = 1e14*u.Msun,
                         Tsys_NEFD = 72.5*u.mJy*u.s**(1./2),
                         Nfeeds = 120,
                         beam_FWHM = 48*u.arcsec,
                         Delta_nu = 40*u.GHz, # was 20 before, why?
                         dnu = 1.11*u.GHz,
                         tobs = 2000*u.hr,
                         Omega_field = 16*u.deg**2,
                         do_Jysr = True,
                         Nfield = 4,
                         catalogue_file = '/mnt/scratch-lustre/horlaville/nate_sims/may_5/gaussian/gaussian_lc_M_gt_1.129e+11Msol.npz'
                         )




Silva_m3_z6_CCATp = dict(model_type = 'ML',
                         model_name = 'SilvaCII',
                         model_par = {'a':0.8727,'b':6.7250},
                         nu = 1897*u.GHz,
                         nuObs = 250*u.GHz,
                         Mmin = 1e10*u.Msun,
                         Mmax = 1e14*u.Msun,
                         Tsys_NEFD=48.5*u.mJy*u.s**(1./2),
                         Nfeeds=1.7,
                         beam_FWHM=46*u.arcsec,
                         Delta_nu=60*u.GHz, # was 20 before, why?
                         dnu=2.7*u.GHz,
                         tobs=4000*u.hr,
                         Omega_field=2*u.deg**2,
                         do_Jysr=True,
                         Nfield = 4,
                         catalogue_file = '/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz'
                         )

Silva_m4_z6_CCATp = dict(model_type = 'ML',
                         model_name = 'SilvaCII',
                         model_par = {'a':0.9231,'b':6.5234},
                         nu = 1897*u.GHz,
                         nuObs = 250*u.GHz,
                         Mmin = 1e10*u.Msun,
                         Mmax = 1e14*u.Msun,
                         Tsys_NEFD=48.5*u.mJy*u.s**(1./2),
                         Nfeeds=1.7,
                         beam_FWHM=46*u.arcsec,
                         Delta_nu=60*u.GHz, # was 20 before, why?
                         dnu=2.7*u.GHz,
                         tobs=4000*u.hr,
                         Omega_field=2*u.deg**2,
                         do_Jysr=True,
                         Nfield = 4,
                         catalogue_file = '/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz'
                         )




                         
# Silva et al. (2015) CO model with CII Stage II parameters
Silva_CO32_StageII = dict(model_type = 'ML',
                          model_name = 'SilvaCO',
                          
                          # CO(2-1)
                          model_par = {'L0':4.7e-29,'M1':1.0e11*u.Msun,
                                       'M2':6.0e11*u.Msun,'M3':5e12*u.Msun,
                                       'M4':5.0e14*u.Msun,'d0':3.05,'d1':-2.0,
                                       'd2':-2.3,'d3':1.9,'d4':5.0},
                          nu = 230.5*u.GHz,
                          nuObs = 215.25*u.GHz,
                          Delta_nu = 69.5*u.GHz,
                          
                          # CO(3-2)
                          # model_par = {'L0':3.0e-24,'M1':6.0e11*u.Msun,
                          # 'M2':5.0e12*u.Msun,'M3':4e13*u.Msun,
                          # 'M4':1.0*u.Msun,'d0':2.6,'d1':-3.5,'d2':0.2,
                          # 'd3':2.2,'d4':0.0},
                          # nu = 345.8*u.GHz,
                          # nuObs = 250*u.GHz,
                          # Delta_nu = 100*u.GHz,
                          
                          Mmin = 1e8*u.Msun,
                          Mmax = 1e14*u.Msun,
                          Tsys_NEFD = 5*u.mJy*u.s**(1./2),
                          Nfeeds = 16000,
                          beam_FWHM = 38*u.arcsec,
                          dnu = 0.4*u.GHz,
                          tobs = 2000*u.hr,
                          Omega_field = 10*u.deg**2,
                          do_Jysr=True)

# Hamsa's CII parameters from her paper with parameter for George's sim data
# most of these parameters are just copied directly from TonyLi 2016 parameters
Hamsa_CII = dict(cosmo_model = 'Planck15',
                  model_type = 'ML',
                  model_name = 'Padmanabhan_CII',
                  model_par = dict(M1=2.39e-5*u.Msun,N1=4.19e11*u.Msun,alpha=1.79,beta=0.49), # from Hamsa's paper
                  nu = 1897*u.GHz,       # dummy (from George's CO sim) --> USE 1897 GHz for CII line always!
                  nuObs = 296*u.GHz,     # dummy
                  Mmin = 4e9*u.Msun,    # George's sim used Mmin = 4e9; Hamsa's paper used Mmin = 1e10 h^-1 Msun
                  Mmax = 1e15*u.Msun,        
                  nM = 5000,
                  nL = 5000,
                  kmin = 1e-2/u.Mpc,
                  kmax = 7./u.Mpc,
                  nk = 100,
                  sigma_scatter = 0.3,
                  Tsys_NEFD = 40*u.K,
                  Nfeeds = 19,
                  beam_FWHM = 4.1*u.arcmin,
                  Delta_nu = 60*u.GHz,
                  dnu = 10.*u.GHz,
                  tobs = 1500*u.hr,
                  Omega_field = 4.*u.deg**2, # from Patrick
                  Nfield = 4,
                  catalogue_file = '/home/dongwooc/scratchspace/gstein_COMAP_backup/z4.77-6.21_700Mpc/npz/COMAP_z4.77-6.21_700Mpc_seed_13579.npz',
                  do_Jysr = True
                  )
              
                                 
              
