from __future__ import absolute_import, print_function
import numpy as np
import scipy
from .tools import * # . tools? 
from . import debug
import time

#@timeme
def load_peakpatch_catalogue(halo_info, filetype='.npz', saveHalos=False, saveFolder='/outputs/'):
    """
    Load peak patch halo catalogue into halos class and cosmology into cosmo class
    
    Slightly modified to work with the lim functions, cosmo class split into separate
    function
    
    When save=True, the function will save the Mcen, Msat, cen_pos, sat_pos output to the saveFolder directory.
    
    Default filetype = '.npz'
    Can take filetype = '.h5'

    Returns
    -------
    halos : class
        Contains all halo information (position, redshift, etc..)
    """
    
    # These are the directories from when saveHalos=True. 
    
    Mcen_file    = saveFolder + 'Mcen.dat'    
    Msat_file    = saveFolder + 'Msat.dat'        
    cen_pos_file = saveFolder + 'cen_pos.dat'   
    sat_pos_file = saveFolder + 'sat_pos.dat'        

    halos       = empty_table()                      # creates empty class to put any halo info into 
    
    if filetype=='.h5':
        # Martine's catalogues use these
        from astropy.cosmology import Planck15 as cosmo15
        h = cosmo15.H(0).value/100
        
        cen_x_fov = 0.  # set zero for now
        cen_y_fov = 0.
        
        if saveHalos:
            print('Saving mass values...')
            mcen       = halo_info['mass_cen'].value
            rm_cenmass = (mcen > 2e13/h)
            mcen = mcen[rm_cenmass]            
            xcen = (halo_info['xpos_cen'].value)[rm_cenmass]
            ycen = (halo_info['ypos_cen'].value)[rm_cenmass]
            zcen = (halo_info['zpos_cen'].value)[rm_cenmass]
            
            with open(Mcen_file, 'w+') as mc:
                np.savetxt(mc, mcen)
            with open(cen_pos_file, 'w+') as cp:
                np.savetxt(cp, (xcen,ycen,zcen))
            
            print('Saving position values...')
            msat       = halo_info['mass_sat'].value
            rm_satmass = (msat > 2e13/h)
            msat = msat[rm_satmass]
            xsat = (halo_info['xpos_sat'].value)[rm_satmass]
            ysat = (halo_info['ypos_sat'].value)[rm_satmass]
            zsat = (halo_info['zpos_sat'].value)[rm_satmass]
            
            with open(Msat_file, 'w+') as ms:
                np.savetxt(ms, msat)
            with open(sat_pos_file, 'w+') as sp:
                np.savetxt(sp, (xsat,ysat,zsat))
            
        else:
            try:
                print('Loading mass files... (takes 3 minutes)') 
                mcen = np.loadtxt(Mcen_file)
                msat = np.loadtxt(Msat_file)
                print('Loading position files...')
                xcen, ycen, zcen = np.loadtxt(cen_pos_file)[0], np.loadtxt(cen_pos_file)[1], np.loadtxt(cen_pos_file)[2]
                xsat, ysat, zsat = np.loadtxt(sat_pos_file)[0], np.loadtxt(sat_pos_file)[1], np.loadtxt(sat_pos_file)[2]
            except IOError:
                print('Need to set saveHalo=True the first run in order to use saveHalo=False subsequently.')
        
        halos.Mcen = mcen
        halos.x_pos_cen = xcen
        halos.y_pos_cen = ycen
        halos.z_pos_cen = zcen
        
        halos.Msat = msat
        halos.x_pos_sat = xsat
        halos.y_pos_sat = ysat
        halos.z_pos_sat = zsat
        
        # both centrals and satellites together --- from here, it doesn't take too long
        halos.M     = np.append(mcen, msat)
        halos.x_pos = np.append(xcen, xsat)
        halos.y_pos = np.append(ycen, ysat)
        halos.z_pos = np.append(zcen, zsat)
        
        halos.chi_cen = np.sqrt(halos.x_pos_cen**2+halos.y_pos_cen**2+halos.z_pos_cen**2)
        halos.chi_sat = np.sqrt(halos.x_pos_sat**2+halos.y_pos_sat**2+halos.z_pos_sat**2)
        halos.chi     = np.append(halos.chi_cen, halos.chi_sat)  
        
        # manually calculate redshift - taken straight from tools.py chi_to_redshift(chi, cosmo)
        zinterp    = np.linspace(0,10,20000)
        dz         = zinterp[1]-zinterp[0]
        chiinterp  = np.cumsum(drdz(zinterp,cosmo15.h,cosmo15.Om0) * dz)
        chiinterp -= chiinterp[0]
        z_of_chi   = scipy.interpolate.interp1d(chiinterp,zinterp)
        
        halos.redshift_cen = z_of_chi(halos.chi_cen)
        halos.redshift_sat = z_of_chi(halos.chi_sat)
        halos.redshift     = z_of_chi(halos.chi)
        
        halos.nhalo_cen = len(halos.Mcen)
        halos.nhalo_sat = len(halos.Msat)
        halos.nhalo     = len(halos.M)
        
        halos.ra_cen    = np.arctan2(-halos.x_pos_cen,halos.z_pos_cen)*180./np.pi - cen_x_fov
        halos.ra_sat    = np.arctan2(-halos.x_pos_sat,halos.z_pos_sat)*180./np.pi - cen_x_fov
        halos.ra        = np.arctan2(-halos.x_pos,halos.z_pos)*180./np.pi - cen_x_fov
        
        halos.dec_cen   = np.arcsin(  halos.y_pos_cen/halos.chi_cen  )*180./np.pi - cen_y_fov
        halos.dec_sat   = np.arcsin(  halos.y_pos_sat/halos.chi_sat  )*180./np.pi - cen_y_fov
        halos.dec       = np.arcsin(  halos.y_pos/halos.chi  )*180./np.pi - cen_y_fov
        
        
        # Using 'ns' and 'sigma8' the same as George's 'cosmo_header'
        params_dict = {'h': h, 'sigma8': 0.82, 'Omega_M': cosmo15.Om0, 'Omega_L': cosmo15.Ode0, 'Omega_B': cosmo15.Ob0, 'ns': 0.96}
        
        # George's 'cosmo_header', in case I need to manually set these values later:
        # params_dict = {'h': 0.7, 'sigma8': 0.82, 'Omega_M': 0.286, 'Omega_L': 0.714, 'Omega_B': 0.047, 'ns': 0.96}

        
    else:
        print('Loading .npz catalogues...')
        params_dict = halo_info['cosmo_header'][()]     
        if debug.verbose: print("\thalo catalogue contains:\n\t\t", halo_info.files)    #.npz version
            
        cen_x_fov  = params_dict.get('cen_x_fov', 0.) #if the halo catalogue is not centered along the z axis
        cen_y_fov  = params_dict.get('cen_y_fov', 0.) #if the halo catalogue is not centered along the z axis
        
        halos.M          = halo_info['M']     # halo mass in Msun 
        halos.x_pos      = halo_info['x']     # halo x position in comoving Mpc 
        halos.y_pos      = halo_info['y']     # halo y position in comoving Mpc 
        halos.z_pos      = halo_info['z']     # halo z position in comoving Mpc 
        
        halos.vx         = halo_info['vx']    # halo x velocity in km/s
        halos.vy         = halo_info['vy']    # halo y velocity in km/s
        halos.vz         = halo_info['vz']    # halo z velocity in km/s
        halos.redshift   = halo_info['zhalo'] # observed redshift incl velocities
        halos.zformation = halo_info['zform'] # formation redshift of halo
        halos.chi        = np.sqrt(halos.x_pos**2+halos.y_pos**2+halos.z_pos**2)  
        halos.nhalo      = len(halos.M)
        halos.ra         = np.arctan2(-halos.x_pos,halos.z_pos)*180./np.pi - cen_x_fov
        halos.dec        = np.arcsin(  halos.y_pos/halos.chi  )*180./np.pi - cen_y_fov

    assert np.max(halos.M) < 1.e17,             "Halos seem too massive"
    # assert np.max(halos.redshift) < 4.,         "need to change max redshift interpolation in tools.py" 
    assert np.max(halos.redshift) < 10,        "need to change max redshift interpolation in tools.py" # changed for George's z=4.77-6.21 sim data because of error in calling m.maps (May 17) - Clara
    if debug.verbose: print('\n\t%d halos loaded' % halos.nhalo)
    
    return halos
    


def load_peakpatch_catalogue_cosmo(halo_info):
    """
    Load peak patch cosmology into cosmo class
    
    Slightly modified to work with the lim functions, halo class split into separate
    function

    Returns
    -------
    cosmo : class
        Contains all cosmology information (Omega_i, sigme_8, etc)
    """
    cosmo      = empty_table()            # creates empty class to put any cosmology info into  

    if debug.verbose: print("\thalo catalogue contains:\n\t\t", halo_info.files)
    
    #get cosmology from halo catalogue
    params_dict    = halo_info['cosmo_header'][()]
    cosmo.Omega_M  = params_dict.get('Omega_M')
    cosmo.Omega_B  = params_dict.get('Omega_B')
    cosmo.Omega_L  = params_dict.get('Omega_L')
    cosmo.h        = params_dict.get('h'      )
    cosmo.ns       = params_dict.get('ns'     )
    cosmo.sigma8   = params_dict.get('sigma8' )

    assert (cosmo.Omega_M + cosmo.Omega_L)==1., "Does not seem to be flat universe cosmology" 

    if debug.verbose: print('\n\t%d halos loaded' % halos.nhalo)

    return cosmo
    


#@timeme
def cull_peakpatch_catalogue(halos, min_mass, max_mass, mapinst, haloType='all'): 
    """
    crops the halo catalogue to only include desired halos
    haloType determines whether you want to use all halos (centrals and satellites), only centrals,
    or only satellites. Options are:
    haloType='all' for both
            ='cen' for only centrals
            ='sat' for only satellites
    """
    
    dir_halos = dir(halos)
    del_attr  = []
    
    if haloType == 'all':
        Mass     = halos.M
        Redshift = halos.redshift
        RA       = halos.ra
        DEC      = halos.dec
        halos.nhalo_all = halos.nhalo
        for attr in dir_halos:      # remove only centrals and satellites attributes
            if (('cen' in attr) or ('sat' in attr)):
                del_attr.append(attr)
    elif haloType == 'cen':
        Mass     = halos.Mcen
        Redshift = halos.redshift_cen
        RA       = halos.ra_cen
        DEC      = halos.dec_cen
        for attr in dir_halos:      # remove attributes that are not centrals
            if 'cen' not in attr:
                del_attr.append(attr)
    elif haloType == 'sat':
        Mass     = halos.Msat
        Redshift = halos.redshift_sat
        RA       = halos.ra_sat
        DEC      = halos.dec_sat
        for attr in dir_halos:      # remove attributes that are not satellites
            if 'sat' not in attr:
                del_attr.append(attr)
    else:
        raise Exception("haloType only takes 'all', 'cen', or 'sat'.")
    
    [dir_halos.remove(attr) for attr in del_attr]
    
    dm = (Mass > min_mass)  * (Redshift    >= mapinst.z_i)\
                            * (np.abs(RA)  <= mapinst.fov_x/2)\
                            * (np.abs(DEC) <= mapinst.fov_y/2)\
                            * (Redshift    <= mapinst.z_f)\
                            * (Mass        <  max_mass)
    
    for i in dir_halos:
        if i[0]=='_': continue
        try:
            setattr(halos,i,getattr(halos,i)[dm])
        except TypeError:
            pass
    
    if haloType == 'cen':
        halos.nhalo    = len(halos.Mcen)
        halos.M        = halos.Mcen
        halos.chi      = halos.chi_cen
        halos.ra       = halos.ra_cen
        halos.dec      = halos.dec_cen
        halos.redshift = halos.redshift_cen
        halos.x_pos    = halos.x_pos_cen
        halos.y_pos    = halos.y_pos_cen
        halos.z_pos    = halos.z_pos_cen
    elif haloType == 'sat':
        halos.nhalo    = len(halos.Msat)
        halos.M        = halos.Msat
        halos.chi      = halos.chi_sat
        halos.ra       = halos.ra_sat
        halos.dec      = halos.dec_sat
        halos.redshift = halos.redshift_sat
        halos.x_pos    = halos.x_pos_sat
        halos.y_pos    = halos.y_pos_sat
        halos.z_pos    = halos.z_pos_sat
    else:
        halos.nhalo = len(halos.M)
    
    if debug.verbose: print('\n\t%d halos remain after mass/map cut' % halos.nhalo)

    return halos
