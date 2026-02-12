import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.units as u
from scipy import special

from lim import lim

matplotlib.rcParams.update({'font.size': 18,'figure.figsize':[8,7]}) 
from scipy.ndimage import gaussian_filter

def own_noise(m, sign, sigm, pix=1.):
    
    '''
    m is the limlam object
    sign is the 3D cube of intensities
    sigm is the noise level
    pix is the pixelization blurring desired. default is 1.
    '''
    
    sm_map = gaussian_filter(sign, [pix, pix, 0])
    noise_sigma  = sigm
    noise_map    = np.random.normal(0, 
                                    noise_sigma.to(u.Jy/u.sr, equivalencies=u.brightness_temperature(m.nuObs)).value,
                                    sign.shape)
    # print(noise_sigma.to(u.Jy/u.sr).value)
    sm_noise_map = sm_map + noise_map
    
    return sm_noise_map

def pix_res(beam_width, map_dim_deg, map_dim_pix):
    
    '''
    
    computes the pixel size corresponding to beaming the signal with a certain angle (beam_width) 
    given the size of the map in degrees (map_dim_deg) and pixels (map_dim_pix) 
    
    beam_width; size of t
    map_dim_pix in pixels
    map_dim_deg in degrees
    
    '''
    
    ang_res = (beam_width.to(u.arcmin))/2.355
    map_dim_deg.to(u.deg)
    
    
    return ang_res/((map_dim_deg/map_dim_pix)*(60*u.arcmin/u.deg))



def halo_centpix(lim_obj, halo_xpos, halo_ypos, halo_zpos):

    '''
    Parameters
    ----------

    halo_xpos: array_like
               list of RA coordinates of halos, from the original distribution of halos

    halo_ypos: array_like
               list of DEC coordinates of halos, from the original distribution of halos


    Returns
    -------

    halo_centpix_x: array_like
                    list of the x-pixel positions of halos, mapped on the lim map


    halo_centpix_x: array_like
                    list of the y-pixel positions of halos, mapped on the lim map

    '''

    map_xs = lim_obj.mapinst.pix_bincents_x
    map_ys = lim_obj.mapinst.pix_bincents_y
    map_zs = (lim_obj.mapinst.nu_rest/lim_obj.mapinst.nu_bincents) - 1
    
    pixcents_x_mesh, halo_xs_mesh = np.meshgrid(map_xs, halo_xpos)
    halo_centpix_x = np.argmin(np.abs(halo_xs_mesh - pixcents_x_mesh), axis=1)
    
    pixcents_y_mesh, halo_ys_mesh = np.meshgrid(map_ys, halo_ypos)
    halo_centpix_y = np.argmin(np.abs(halo_ys_mesh - pixcents_y_mesh), axis=1)
    
    pixcents_z_mesh, halo_zs_mesh = np.meshgrid(map_zs, halo_zpos)
    halo_centpix_z = np.argmin(np.abs(halo_zs_mesh - pixcents_z_mesh), axis=1)
    
    return halo_centpix_x, halo_centpix_y, halo_centpix_z







def halo_map(lim_obj, n, halo_xpos, halo_ypos, halo_zpos):

    halo_centpix_x, halo_centpix_y, halo_centpix_z = halo_centpix(lim_obj, halo_xpos, halo_ypos, halo_zpos)

    halo_mapx = np.linspace(halo_centpix_x - ((n - 1)/2), halo_centpix_x + ((n - 1)/2), n, axis = 1)
    halo_mapy = np.linspace(halo_centpix_y - ((n - 1)/2), halo_centpix_y + ((n - 1)/2), n, axis = 1)
    halo_mapz = halo_centpix_z
    
    npix_x, npix_y = lim_obj.mapinst.npix_x, lim_obj.mapinst.npix_y
    outb_x = halo_mapx >= npix_x
    outb_y = halo_mapy >= npix_y
    
    halo_mapx[outb_x] = None
    halo_mapy[outb_y] = None
    
    return halo_mapx, halo_mapy, halo_mapz


def inbound_halos(lim_obj, n, halo_xpos, halo_ypos, halo_zpos):
    
    halo_mapx, halo_mapy, halo_mapz = halo_map(lim_obj, n, halo_xpos, halo_ypos, halo_zpos)

    inb_x = ~np.isnan(halo_mapx).any(axis = 1)
    inb_y = ~np.isnan(halo_mapy).any(axis = 1)
    
    inb = np.logical_and(inb_x, inb_y)
            
    return halo_mapx, halo_mapy, halo_mapz, inb



def lum_hod(lim_obj, n, halo_xpos, halo_ypos, halo_zpos):

    halo_mapx, halo_mapy, halo_mapz, inb = inbound_halos(lim_obj, n, halo_xpos, halo_ypos, halo_zpos)
    
    inb_mapx = halo_mapx[inb].astype(int)
    inb_mapy = halo_mapy[inb].astype(int)
    inb_mapz = halo_mapz[inb].astype(int)

    nhalos = len(inb_mapx)
    
    pure_map = lim_obj.maps.value
    noisy_map= lim_obj.noise_added_map
    
    grid = [0 for i in range(nhalos)]
    sigs = [0 for i in range(nhalos)]
    noisy= [0 for i in range(nhalos)]
    
    for i in range(nhalos):
    
        grid[i] = np.meshgrid(inb_mapx[i], inb_mapy[i], inb_mapz[i])
        sigs[i] = pure_map[grid[i][0], grid[i][1], grid[i][2]]
        noisy[i]= noisy_map[grid[i][0], grid[i][1], grid[i][2]]
        
    sigs = np.reshape(sigs, (nhalos, n, n))
    noisy= np.reshape(noisy, (nhalos, n, n))
    
    return sigs, noisy, inb

def lum_hod_corr(lim_obj, n, halo_xpos, halo_ypos, halo_zpos):
    
    sig_scale = 1.23

    halo_mapx, halo_mapy, halo_mapz, inb = inbound_halos(lim_obj, n, halo_xpos, halo_ypos, halo_zpos)
    
    inb_mapx = halo_mapx[inb].astype(int)
    inb_mapy = halo_mapy[inb].astype(int)
    inb_mapz = halo_mapz[inb].astype(int)

    nhalos = len(inb_mapx)
    
    pure_map = lim_obj.maps.value
    pure_map = sig_scale*pure_map
    
    noisy_map= own_noise(lim_obj, pure_map, lim_obj.sigma_N)
    
    grid = [0 for i in range(nhalos)]
    sigs = [0 for i in range(nhalos)]
    noisy= [0 for i in range(nhalos)]
    
    for i in range(nhalos):
    
        grid[i] = np.meshgrid(inb_mapx[i], inb_mapy[i], inb_mapz[i])
        sigs[i] = pure_map[grid[i][0], grid[i][1], grid[i][2]]
        noisy[i]= noisy_map[grid[i][0], grid[i][1], grid[i][2]]
        
    sigs = np.reshape(sigs, (nhalos, n, n))
    noisy= np.reshape(noisy, (nhalos, n, n))
    
    return sigs, noisy, inb

# Adding HOD functions:

def n_mh(m_h, logm_min, sigma_logm, alpha, dc):
    
    log_m1p = 1.18*logm_min - 1.28
    log_m0 = 0.76*log_m1p + 2.3
    
    m_1p = 10**(log_m1p)
    m_0 = 10**(log_m0)
    
    
    n_c = (1/2)*(1 + special.erf((np.log10(m_h) - logm_min)/sigma_logm))
    n_s = (n_c)*((m_h - m_0)/m_1p)**alpha
    
    n = dc*(n_c + n_s)
                 
    return n
