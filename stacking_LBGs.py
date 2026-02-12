from stacking_params import *
plt.rcParams["mathtext.fontset"] = "dejavuserif"


# Calculating luminosities

halo_ms = lim_sim.halos.M

mthresh = halo_ms > mass_cut

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


pure_map, noisy_map = lum(lim_sim, n, halo_xs, halo_ys, halo_zs)
pure_map, noisy_map = np.ma.masked_array(pure_map, np.isnan(pure_map)), np.ma.masked_array(noisy_map, np.isnan(noisy_map))

pure_stack, noisy_stack = np.average(pure_map, axis = 0, weights = draws), np.average(noisy_map, axis = 0, weights = draws)


# Plotting

fig , axes = plt.subplots(nrows = 2, ncols = 1, figsize = (10, 16))

plt.subplot(211)
plt.imshow(gaussian_filter(pure_stack, beam_res), cmap = 'CMRmap', extent = [-stack_dim/2, stack_dim/2, -stack_dim/2, stack_dim/2])
plt.title(r'$Beamed\ Pure\ Signal\ Stacked\ Map$')
plt.ylabel(r'$Dec\ (Degrees)$')
plt.colorbar(label = r'$[C_{II}]\ Luminosity\ (Jy/sr)$')


plt.subplot(212)
plt.imshow(noisy_stack, cmap = 'CMRmap', extent = [-stack_dim/2, stack_dim/2, -stack_dim/2, stack_dim/2])
plt.title(r'$Beamed\ Forecast\ Stacked\ Map$')
plt.xlabel(r'$RA\ (Degrees)$')
plt.ylabel(r'$Dec\ (Degrees)$')
plt.colorbar(label = r'$[C_{II}]\ Luminosity\ (Jy/sr)$')

plt.show()

plt.savefig('Stacking/aug2_lbg_test.png', bbox_inches = 'tight')



    




# forget below

# import sys

# lightcone_paths = [...]

# for i in range(len(lightcone_paths)):
#    sys.argv = ['stacking.py','lightcone_path[i]', 'arg2']
#    execfile('stacking.py')

# stacking would need to maybe write and save the stacked map luminosities?


