from peakpatchtools import PeakPatch

# Add parameters to the lim object, each parameter is explained in the following 
# markdown box 
def update_map_COMAP_Fid( 
    
        # Peak Patch run
        r : PeakPatch | None ,
        nu_binning  = None,
    
        # Telescope specs
        Omega_field = None,#'COMAP', 
        tobs        = 29000*u.h,     # observing time (2111.05933, end of page 3)        
        Nfeeds      = 19,            # number of feeds (2111.05927, Table 2)
        Nfield      = 3,             # number of fields (2111.05929, Table 1)
        Mmin        = 1.0e8 *u.Msun, # minimum DM halo mass
        Mmax        = 1.0e15*u.Msun, # maximum DM halo mass
        
        # CO model
        model_name = 'COMAP_Fid', # 'LichenCII_v3',
        band       = 'Ku',
        line       = 'CO(1-0)',
    
        # model_type = 'ML' to use halo mass function or 'LF' to use luminosity function
        model_type = 'ML',
    
        # Set LIM model parameters for the COMAP fiducial cosmology, the values are
        # taken as the realistic+ prediction, shown as data-driven prior “UM+COLDz+COPSS”
        # in Table 1 of 2111.05931. Model parameters are gone over in more detail in
        # appendix A, the actual L_{CO} function is eq. 8
        model_par = { 
            'A'  : -2.85,
            'B'  : -0.42,
            'C'  :  10**10.63,
            'Ms' :  10**12.3 # in solar masses
        },

        catalogue_file = None
        ):
    # r is a peak patch run object
    # Omega_field is the solid angle of RA x Dec x nu cubes
    # nu_binning determines how nu bins are made. This allows us to ensure all bins lie
    #     fully within the redshift range of the peak patch runs nu'_max <= nu/(1+z_min)
    #     and nu'_min >= nu/(1+z_max). Here are options
    #     None = just does its thing
    #     'centre' = bins not exceeding PP bounds such that 
    #         (nu'_max + nu'_min)/2 = ( nu/(1+z_min) + nu/(1+z_min) )/2
    #     'high' = bins not exceeding PP bounds such that nu'_max = nu/(1+z_min)
    #     'low' = bins not exceeding PP bounds such that nu'_min = nu/(1+z_max)
    
    # Automatically set the field of view for the cube
    if Omega_field is None:
        Omega_field = pkp.solid_angle_subtended_by_square(
                r.boxsize, r.cenz+r.boxsize/2, u.deg**2 )
    elif Omega_field == 'COMAP':
        Omega_field = 4.0 * u.deg**2

    # make sure COMAP is sensitive to the emission line in this band
    if ( ( line=='CO(2-1)' and band=='Ku' ) 
            or ( line=='CO(3-2)' and ( band=='Ku' or band=='Ka') ) ):
        raise ValueError('COMAP is not sensitive to the {0} line in the {1} band.'
                         .format(line,band) )
        
    # rest frame CO(1-0) emission line
    if   line == 'CO(1-0)': nu = 115.2712*u.GHz
    elif line == 'CO(2-1)': nu = 230.5380*u.GHz
    elif line == 'CO(3-2)': nu = 345.7960*u.GHz
    else: 
        raise ValueError(('emission line {0} not recognised. Supported values are:'+
                         '\n\'CO(1-0)\', \'CO(2-1)\' and \'CO(3-2)\'.').format(line))
        
    # redshift of an emission line is given by
    def nu_obs(nu,z):           return nu/(1+z)
    def z_of_nu_obs(nu,nu_obs): return nu/nu_obs-1

    # frequency range of the intensity cube, the angular resolution (see 2309.03184,
    # first paragraph of page 11), and Noise-equivalent flux density estimate (based
    # on a wild guess)
    if band == 'Ku': 
        nu_min, nu_max = 14.0*u.GHz, 16.0*u.GHz
        beam_FWHM      = 3.7*u.arcmin
        T_sys          = 22.*u.K
    elif band == 'Ka': 
        nu_min, nu_max = 28.0*u.GHz, 32.0*u.GHz
        beam_FWHM      = 3.9*u.arcmin
        T_sys          = 44.*u.K
    elif band == 'Q' : 
        nu_min, nu_max = 42.0*u.GHz, 48.0*u.GHz
        beam_FWHM      = 4.0*u.arcmin
        T_sys          = 60.*u.K
    else:
        raise ValueError(('frequency band {0} not recognised. Supported values are:'+
                          '\nKu, Ka and Q.').format(band))
    
    # redshift range of the intensity cube
    zmin, zmax = z_of_nu_obs(nu,nu_max) , z_of_nu_obs(nu,nu_min)
    
    # For COMAP, frequency channelisation is (taken from Table 2 of 2111.05927)
    dnu = 31.25*u.MHz
    
    # Minimum and maximum frequencies at which we observe the redshifted CO
    # emission line for this band
    nu_min , nu_max = nu_obs(nu,zmax) , nu_obs(nu,zmin)
    
    # check frequency range by comparing against Peak Patch box z range
    print(r.boxsize,r.cenz)
    chi = np.array([-.5,.5])*r.boxsize+r.cenz
    z_chi_tab = pkp.z_of_r_comoving_table( Omega_m0=r.OmB+r.Omx,
                                           Omega_Lambda=r.Omvac, H_0=r.h*100 )
    zmin_pp,zmax_pp = pkp.z_of_r_comoving( chi, z_chi_table=z_chi_tab )
    print(chi,zmin_pp,zmax_pp)
    if (zmin > zmax_pp) or (zmax < zmin_pp):
        raise ValueError( ('Peak Patch z range [{0},{1}] does not coincide with '+
                           'line/band redshfits [{2},{3}].'
                          ).format(zmin_pp,zmax_pp,zmin,zmax) ) 
    elif ( zmin_pp > zmin ) or ( zmax_pp < zmax ):
        print( ('Warning: Peak Patch range [{0},{1}] only partially coincides with'+
                ' line/band redshifts [{2},{3}].'
               ).format(zmin_pp,zmax_pp,zmin,zmax) )
    
    if nu_binning==None:
        
        # Observed frequency range for redshifted CII emission
        Delta_nu = nu_max-nu_min

        # Mean frequency for CII emission from halos in this catalogue 
        nuObs = (nu_max+nu_min)/2
    
    else:
        
        N_bins = (nu_max-nu_min) // dnu
            
        # Mean frequency for CII emission from halos in this catalogue 
        nuObs = (nu_max+nu_min)/2
    
        if nu_binning == 'centre':
            nu_min,nu_max = nuObs-N_bins/2*dnu , nuObs+N_bins/2*dnu

        elif nu_binning == 'low' or nu_binning == 'high-z':
            nu_max = nu_min + N_bins * dnu
            nuObs  = (nu_max+nu_min)/2
            
        elif nu_binning == 'high' or nu_binning == 'low-z':
            nu_min = nu_max - N_bins * dnu
            nuObs  = (nu_max+nu_min)/2

        elif nu_binning == 'centre/overlap':
            temp_max = nu_min + N_bins * dnu
            temp_min = nu_max - N_bins * dnu
            if temp_max < nu_max:
                nu_max = temp_max + dnu
            if temp_min > nu_min:
                nu_min = temp_min - dnu
            nuObs = (nu_max+nu_min)/2
            
        elif nu_binning == 'low/overlap' or nu_binning == 'high-z/overlap':
            temp_max = nu_min + N_bins * dnu
            if temp_max < nu_max:
                nu_max = temp_max + dnu
            nuObs = (nu_max+nu_min)/2
            
        elif nu_binning == 'high/overlap' or nu_binning == 'low-z/overlap':
            temp_min = nu_max - N_bins * dnu
            if temp_min > nu_min:
                nu_min = temp_min - dnu
            nuObs = (nu_max+nu_min)/2
            
        else:
            raise ValueError('nu_binning={0} not recognised.'.format(nu_binning))
        
        # Observed frequency range for redshifted CII emission
        Delta_nu = nu_max-nu_min
    
    # Dark matter halo catalogue file
    if catalogue_file is None:
        if not hasattr(r,'catalogue_file'):
            r.add_halos()
        catalogue_file = f'{r.catalogue_file}.npz'
    
    # Create objects of class lim 
    m_co = lim( model_params='COMAP_Fid', doSim=True )
    
    # Add parameters to the lim object, each parameter is explained in the following 
    # markdown box 
    m_co.update(

        # CO model
        model_name = model_name,
    
        # model_type = 'ML' to use halo mass function or 'LF' to use luminosity function
        model_type = model_type,
    
        # Set LIM model parameters
        model_par = model_par,
    
        # Dark matter halo catalogue file
        catalogue_file = catalogue_file,
    
        # Observatory parameters
        nu          = nu,          # rest-frame line emission frequency
        dnu         = dnu,         # COMAP frequency channelization
        nuObs       = nuObs,       # median observed line frequency
        Delta_nu    = Delta_nu,    # redshifted frequency band
        tobs        = tobs,        # observing time
        Omega_field = Omega_field, # Observatory field of view
        Tsys_NEFD   = T_sys,       # noise equivalent flux density
        beam_FWHM   = beam_FWHM,   # telescope beam FWHM width
        Nfeeds      = Nfeeds,      # number of feeds in the spectrometer
        Nfield      = Nfield,      # number of fields on the sky
        Mmin        = Mmin,        # minimum DM halo mass
        Mmax        = Mmax         # maximum DM halo mass
        )
    
    return m_co

# Function for correcting redshift errors in LIM cubes that arise because 
# there is an assumption that the luminosity distance is the same for 
# voxels at different redshift    
def redshift_corrected_cube( m , cube=None ):
    
    # Get redshifts of each map slice
    z_bincents = (m.mapinst.nu_rest/m.mapinst.nu_bincents)-1

    # Approx. comoving distance of bin centres
    chi_bincents  = pkp.r_comoving_of_z( z_bincents )
    
    # The evolution fudge factor
    CLT_z_evol  = chi_bincents **2 * (1+z_bincents)**2
    CLT_z_evol /= np.median(CLT_z_evol)
    
    # Make the cubes!
    if cube is None:
        cube = m.maps
    for j in range(len(cube[0,0,:])):
        cube[:,:,j] /= CLT_z_evol[j]

    return cube

# Function for making map cubes
def make_map_cube_COMAP_Fid( m_co, r, noise=True, save_file=None, 
                        overwrite=False, correct_z=True, verbose=True,
                        stdout=None ):
    
    # Load the catalogue file if no catalogue file has been loaded
    if not hasattr(m_co,'catalogue_file'):
        r.add_halos()
        
    # Set the save directory and file
    if save_file is None:
        save_dir  = r.run_dir+'/maps/'
        save_file = save_dir+'cii_cube.npz'
    else:
        save_dir = os.path.dirname(save_file)
        if save_file[-4:] != '.npz':
            save_file += '.npz'
    
    # Report on status of calculation to screen and to a file
    if verbose:
        from peakpatchtools import my_time,timer
        
        # Format output
        def make_output(std_out,output):
            with open(std_out,'a+') as f:
                f.write(output)
            print(output,end='')
        
        # Make sure standard output file exists
        if stdout is None:
            std_out = save_file[:-4]+'.stdout'
        
        # Make output
        stopwatch = timer()
        stopwatch.start()
        output    = f'Calculation start: {stopwatch.dates[0]}\n'
        make_output(std_out,output)
        
    # Read in an intensity cube that was already made
    if (not overwrite) and os.path.exists(save_file):
        m_co.maps = np.load(save_file)['co_cube']*u.Jy/u.sr
        co_cube   = m_co.maps
        
        # Record run time
        if verbose:
            lap    = stopwatch.lap()
            output = ( f'  Map cube loaded: {stopwatch.dates[-1]}\n'
                     + f'         lap time: {lap}\n'                 )
            make_output(std_out,output)
        
        # Noise added map
        if noise:
            m_co.noise_added_map = np.load(save_file)['co_noise']*u.Jy/u.sr
            co_noise             = m_co.noise_added_map
            
            # Record run time
            if verbose:
                lap    = stopwatch.lap()
                output = ( f'Noise cube loaded: {stopwatch.dates[-1]}\n'
                         + f'         lap time: {lap}\n'                 )
                make_output(std_out,output)
        
    # Otherwise, calculate the cube
    else:
        co_cube = m_co.maps
        
        # Record run time
        if verbose:
            lap    = stopwatch.lap()
            output = ( f'    Map cube made: {stopwatch.dates[-1]}\n'
                     + f'         lap time: {lap}\n'                 )
            make_output(std_out,output)
        
        # Noise added map
        if noise:
            coo_noise = m_co.noise_added_map * u.Jy/u.sr
            
            # Record run time
            if verbose:
                lap    = stopwatch.lap()
                output = ( f'  Noise cube made: {stopwatch.dates[-1]}\n'+
                           f'         lap time: {lap}\n'                 )
                make_output(std_out,output)
            
        # Save maps
        if not noise:
            np.savez( save_file, co_cube=m_co.maps.value )
            
            # Record run time
            if verbose:
                lap    = stopwatch.lap()
                output = ( f'   Map cube saved: {stopwatch.dates[-1]}\n'
                         + f'         lap time: {lap}\n'                 )
                make_output(std_out,output)
                
        else:
            np.savez( save_file, co_cube=m_co.maps.value, co_noise=co_noise.value )
            
            # Record run time
            if verbose:
                lap    = stopwatch.lap()
                output = ( f'      Cubes saved: {stopwatch.dates[-1]}\n'
                         + f'         lap time: {lap}\n'                 )
                make_output(std_out,output)
                
    # Do redshift correction
    if correct_z:
        co_cube = redshift_corrected_cube( m_co, co_cube )[:,:,:-1]
        if noise:
            co_noise = redshift_corrected_cube( m_co, co_noise )[:,:,:-1]
        
    # Return cubes
    if not noise:
        return co_cube
    else:
        return co_cube, co_noise

# Function for making map cubes
def add_map_cube_noise_COMAP_Fid( m_co, r, noise_save_file=None, 
                             overwrite=False, correct_z=True, verbose=True,
                             stdout=None ):
    
    # Check that the map cube has already been calculated or read in from file
    if not hasattr( m_co, 'maps'):
        raise AttributeError('You must first calculate map cube or load from file.')
    else:
        co_cube = m_co.maps
    
    # Set the save directory and file
    if noise_save_file is None:
        save_dir        = r.run_dir+'/maps/'
        noise_save_file = save_dir+'co_cube_noise.npz'
    else:
        save_dir = os.path.dirname(noise_save_file)
        if noise_save_file[-4:] != '.npz':
            noise_save_file     += '.npz'
    
    # Report on status of calculation to screen and to a file
    if verbose:
        from peakpatchtools import my_time,timer
        
        # Format output
        def make_output(std_out,output):
            with open(std_out,'a+') as f:
                f.write(output)
            print(output,end='')
        
        # Make sure standard output file exists
        if stdout is None:
            std_out = noise_save_file[:-4]+'.stdout'
        
        # Make output
        stopwatch = timer()
        stopwatch.start()
        output    = f'Calculation start: {stopwatch.dates[0]}\n'
        make_output(std_out,output)
         
    # Read in an intensity cube that was already made
    if (not overwrite) and os.path.exists(noise_save_file):
        
        m_co.noise_added_map = np.load(noise_save_file)['co_noise']*u.Jy/u.sr
        co_noise             = m_co.noise_added_map
            
        # Record run time
        if verbose:
            lap    = stopwatch.lap()
            output = ( f'Noise cube loaded: {stopwatch.dates[-1]}\n'
                     + f'         lap time: {lap}\n'                 )
            make_output(std_out,output)
        
    # Otherwise, calculate the cube
    else:
        co_noise = m_co.noise_added_map * u.Jy/u.sr
            
        # Record run time
        if verbose:
            lap    = stopwatch.lap()
            output = ( f'  Noise cube made: {stopwatch.dates[-1]}\n'+
                       f'         lap time: {lap}\n'                 )
            make_output(std_out,output)
            
        # Save noise map
        np.savez( noise_save_file, co_noise=co_noise.value )
            
        # Record run time
        if verbose:
            lap    = stopwatch.lap()
            output = ( f'      Cubes saved: {stopwatch.dates[-1]}\n'
                     + f'         lap time: {lap}\n'                 )
            make_output(std_out,output)
                
    # Do redshift correction
    if correct_z:
        co_noise = redshift_corrected_cube( m_co, co_noise )[:,:,:-1]
        
    # Return cube
    return co_noise
