from stacking_params import *
from datetime import datetime
plt.rcParams["mathtext.fontset"] = "dejavuserif"


# Calculating luminosities

# 270-lightcone average:

lc_paths = '/home/dongwooc/scratchspace/pprun_hiz_npz/'
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(lc_paths) if isfile(join(lc_paths, f))]
onlyfiles.remove('pksc2npz_5313591.out')
onlyfiles.remove('pksc2npz.sh')
for i in range(len(onlyfiles)):
    onlyfiles[i] = lc_paths+onlyfiles[i]
  
 
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Beginning time =", current_time)
print("------")

for i in range(len(onlyfiles)):
    
    lim_sim.update(catalogue_file = f"{onlyfiles[i]}")
    
    print('Loading', i, 'th lightcone...')
    
    halo_ms = lim_sim.halos.M
    
    mthresh = halo_ms > mass_cut
    
    print('Finished loading', i, 'th lightcone! Now starting the stack computation...')

    map_zs = (lim_sim.mapinst.nu_rest/lim_sim.mapinst.nu_bincents) - 1

    halo_zs = lim_sim.halos.redshift[mthresh]
    good_halo_zs = np.where(np.logical_and(halo_zs >= z_sel - err, halo_zs <= z_sel + err))
    
    halo_xs = lim_sim.halos.ra[mthresh][good_halo_zs]
    halo_ys = lim_sim.halos.dec[mthresh][good_halo_zs]
    halo_zs = halo_zs[good_halo_zs]
    
    n_gal_tot = n_mh(halo_ms[mthresh][good_halo_zs], log_m_min, sigma_log_m, alph, dutyc)
    draws = np.random.poisson(n_gal_tot)

    print('------------------------')
    print(' - The total forecast observing time has been set to', t_obs, '-')
    print(' - The redshift of selected halos is', round(z_sel, 2), 'and accepted halos are in the redshift range [', round(z_sel - err, 2), ',', round(z_sel + err, 2), '], which accounts for', len(halo_xs), 'halos -')
    print(' - The lightcone is in the redshift range z = [', round(np.min(map_zs), 3), ', ', round(np.max(map_zs), 3), '] -') 
    print(' - Stacked map is', n, 'by', n, 'which covers', stack_dim, 'deg by', stack_dim, 'deg -')
    print(' - We beam the stacked map with a width of', beam_width, ', which corresponds to a Gaussian filter of radius', beam_res, 'pixels - ')
    print(' - The poissonian draw has yielded', np.sum(draws), 'LBG galaxies among our', len(halo_xs), 'halos - ')
    print('------------------------')
    
    pure_map, noisy_map, inb = lum_hod(lim_sim, n, halo_xs, halo_ys, halo_zs)
    
    draws = draws[inb]
    map_size = np.ones_like(pure_map)
    draws = np.reshape(draws, (len(draws), 1, 1))
    sized_draws = draws*map_size
    
    pure_map, noisy_map = np.ma.masked_array(pure_map, np.isnan(pure_map)), np.ma.masked_array(noisy_map, np.isnan(noisy_map))

    pure_stack, noisy_stack = np.average(pure_map, axis = 0, weights = sized_draws), np.average(noisy_map, axis = 0, weights = sized_draws)
    
    pure_stack.dump('/mnt/scratch-lustre/horlaville/nuObs270/zdex04/alpha_cii_0-024/alpha_mhi_0-74/Mmin_2-10e10/alpha0_-1_412/gamma0_0_31/stacks/LBGs/stacks_z5-9_tobs2000h_LBGs/sig/sig'+str(i)+'.npy')
    noisy_stack.dump('/mnt/scratch-lustre/horlaville/nuObs270/zdex04/alpha_cii_0-024/alpha_mhi_0-74/Mmin_2-10e10/alpha0_-1_412/gamma0_0_31/stacks/LBGs/stacks_z5-9_tobs2000h_LBGs/for/for'+str(i)+'.npy')
    
    print('Finished loading', i, 'th stack!')
    
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("------")
print("Finished at =", current_time)



