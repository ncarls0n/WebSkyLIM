# Martine's original code to load .h5 catalogues

import numpy as np

import healpy as hp

import matplotlib.pyplot as plt

import h5py

from scipy.stats import norm

from scipy.stats import uniform

from cosmology import *

from astropy.io import fits

from astropy.cosmology import Planck15 as cosmo, z_at_value

import astropy.units as u

import os

slice_width = 100

test = False

h = cosmo.H(0).value/100

frac_cens = 0.5



def xyz_to_ang_z(x,y,z):

chi = np.sqrt(x**2+y**2+z**2) # Mpc

redshift = zofchi(chi)

theta, phi = hp.vec2ang(np.column_stack((x,y,z))) # in radians

return (theta,phi,redshift)



def get_uniform(mass, stepsize, nperbin):

idx_keep = []

for i in np.arange(2e13/h, 3e13/h, stepsize):

idx_inbin = np.where(np.logical_and((mass>i), (mass<(i+stepsize))))[0]

if len(idx_inbin) != 0:

idx_keep.append(np.random.choice(idx_inbin, size=nperbin))

idx_keep = np.concatenate(idx_keep)

return idx_keep



# mode = 'cmass'

mode = 'redmagic'

filepath = '/scratch/r/rbond/lokken/peak-patch-runs/june_18/galaxy_catalogue.h5'
# filepath = '/mnt/raid-cita/echung/surp2020/lim_Clara/Martines_galaxy_catalogue.h5'



galcat = h5py.File(filepath, 'r')





mcen = galcat['mass_cen'].value



if mode == 'redmagic' and not os.path.exists('/scratch/r/rbond/lokken/peak-patch-runs/june_18/sats_mhalo_uniform_2e13_3e13_overh_chi_2032_3032.npy'):

msat = galcat['mass_sat'].value

rm_satmass = (msat > 2e13/h)

msat = msat[rm_satmass]

xsat = (galcat['xpos_sat'].value)[rm_satmass]

ysat = (galcat['ypos_sat'].value)[rm_satmass]

zsat = (galcat['zpos_sat'].value)[rm_satmass]


rm_cenmass = (mcen > 2e13/h)

mcen = mcen[rm_cenmass]

xcen = (galcat['xpos_cen'].value)[rm_cenmass]

ycen = (galcat['ypos_cen'].value)[rm_cenmass]

zcen = (galcat['zpos_cen'].value)[rm_cenmass]

galcat.close()



if test == True:

subset_cen = np.random.choice(np.arange(len(mcen)), size=5000)

mcen = mcen[subset_cen]; xcen = xcen[subset_cen]; ycen=ycen[subset_cen]; zcen = zcen[subset_cen]

subset_sat = np.random.choice(np.arange(len(msat)), size=5000)

msat = msat[subset_sat]; xsat = xsat[subset_sat]; ysat=ysat[subset_sat]; zsat = zsat[subset_sat]


chicen = np.sqrt(xcen**2 + ycen**2 + zcen**2)

chisat = np.sqrt(xsat**2 + ysat**2 + zsat**2)

# get a tophat distribution from 2e13/h to 3e13/h for each slice in distance 

minz = 0.1

maxz = 0.8

nbins = int((cosmo.comoving_distance(maxz).value - cosmo.comoving_distance(minz).value) // slice_width)



rm_cen = []

rm_sat = []

for i in range(nbins):

dist_slice_min = cosmo.comoving_distance(minz).value + slice_width*i

dist_slice_max = dist_slice_min + slice_width

cen_bin_cond = np.where(np.logical_and((chicen < dist_slice_max), (chicen > dist_slice_min)))[0]

sat_bin_cond = np.where(np.logical_and((chisat < dist_slice_max), (chisat > dist_slice_min)))[0]

idx_cen = np.arange(len(cen_bin_cond))

idx_sat = np.arange(len(sat_bin_cond))

N = 20

stepsize = (3e13/h-2e13/h)/N

# find the number of halos in the highest mass step, which has the least

nperbin_cen = len(mcen[cen_bin_cond] > (3e13/h-stepsize))

nperbin_sat = int(len(msat[sat_bin_cond] > (3e13/h-stepsize))/2.)

if (nperbin_cen != 0) & (nperbin_sat != 0):

print("Selecting {:d} satellites and {:d} centrals per mass bin of size {:.2e} for slice {:.0f} to {:.0f} Mpc.".format(nperbin_cen, nperbin_sat, stepsize, dist_slice_min, dist_slice_max))

idx_keep_cen = get_uniform(mcen[cen_bin_cond], stepsize, nperbin_cen)

idx_keep_sat = get_uniform(msat[sat_bin_cond], stepsize, nperbin_sat)

addcen = np.zeros((len(idx_keep_cen),4))

addcen[:,0] = xcen[cen_bin_cond][idx_keep_cen]

addcen[:,1] = ycen[cen_bin_cond][idx_keep_cen]

addcen[:,2] = zcen[cen_bin_cond][idx_keep_cen]

addcen[:,3] = mcen[cen_bin_cond][idx_keep_cen]

addsat = np.zeros((len(idx_keep_sat),4))

addsat[:,0] = xsat[sat_bin_cond][idx_keep_sat]

addsat[:,1] = ysat[sat_bin_cond][idx_keep_sat]

addsat[:,2] = zsat[sat_bin_cond][idx_keep_sat]

addsat[:,3] = msat[sat_bin_cond][idx_keep_sat]

rm_cen.append(addcen)

rm_sat.append(addsat)

else:

print("No galaxies in this range. Moving on.")

# highm_cen_cond = np.where(mcen > 3e13/h)[0]

# highm_sat_cond = np.where(msat > 3e13/h)[0]

# highm_cen = np.zeros((len(highm_cen_cond),4))

# highm_sat = np.zeros((len(highm_sat_cond),4))

# highm_cen[:,0] = xcen[highm_cen_cond]

# highm_cen[:,1] = ycen[highm_cen_cond]

# highm_cen[:,2] = zcen[highm_cen_cond]

# highm_cen[:,3] = mcen[highm_cen_cond]

# highm_sat[:,0] = xsat[highm_sat_cond]

# highm_sat[:,1] = ysat[highm_sat_cond]

# highm_sat[:,2] = zsat[highm_sat_cond]

# highm_sat[:,3] = msat[highm_sat_cond]

# rm_cen.append(highm_cen)

# rm_sat.append(highm_sat)

rm_cen = np.concatenate(rm_cen)

rm_sat = np.concatenate(rm_sat)

# break up satellite catalog by distance

chisat = np.sqrt(rm_sat[:,0]**2 + rm_sat[:,1]**2 + rm_sat[:,2]**2)

np.save("/scratch/r/rbond/lokken/peak-patch-runs/june_18/sats_mhalo_uniform_2e13_3e13_overh_chi_lt_1032.npy", rm_sat[chisat<1032])

np.save("/scratch/r/rbond/lokken/peak-patch-runs/june_18/sats_mhalo_uniform_2e13_3e13_overh_chi_1032_2032.npy", rm_sat[(chisat>1032) & (chisat<2032)])

np.save("/scratch/r/rbond/lokken/peak-patch-runs/june_18/sats_mhalo_uniform_2e13_3e13_overh_chi_2032_3032.npy", rm_sat[(chisat>2032) & (chisat<3032)])

np.save("/scratch/r/rbond/lokken/peak-patch-runs/june_18/cens_mhalo_uniform_2e13_3e13_overh.npy", rm_cen)

print("Saved files.")

if test == True:

plt.hist(rm_cen[:,3], bins=100)

plt.savefig("/scratch/r/rbond/lokken/peak-patch-runs/june_18/rm_mock_centrals_mh_hist.pdf")

plt.clf()

plt.hist(rm_sat[:,3], bins=100)

plt.savefig("/scratch/r/rbond/lokken/peak-patch-runs/june_18/rm_mock_sats_mh_hist.pdf")

plt.clf()

plt.hist(np.sqrt(rm_cen[:,0]**2+rm_cen[:,1]**2+rm_cen[:,2]**2), bins=100)

plt.savefig("/scratch/r/rbond/lokken/peak-patch-runs/june_18/rm_mock_centrals_chi_hist.pdf")

plt.clf()

plt.hist(np.sqrt(rm_sat[:,0]**2+rm_sat[:,1]**2+rm_sat[:,2]**2), bins=100)

plt.savefig("/scratch/r/rbond/lokken/peak-patch-runs/june_18/rm_mock_sats_chi_hist.pdf")

elif mode == 'redmagic' and os.path.exists('/scratch/r/rbond/lokken/peak-patch-runs/june_18/sats_mhalo_uniform_2e13_3e13_overh_chi_2032_3032.npy'):

print("this one")

with fits.open("/scratch/r/rbond/lokken/data/large_region_y3_gold_2.2.1_wide_sofcol+deep_mof_run_redmagic_highdens.fit") as hdu:

rmdat = hdu[1].data



ra_min = 0.; ra_max = 7.; dec_min = -3.; dec_max = 4.

# select all redmagic within a box of sq deg, and make sure this is filled at redshifts out to 0.8

in_ra = (rmdat['ra'] > ra_min) & (rmdat['ra'] < ra_max)

in_dec = (rmdat['dec'] > dec_min) & (rmdat['dec'] < dec_max)

box_area = (ra_max - ra_min) * (dec_max - dec_min)

rm_inbox = rmdat[in_ra&in_dec]



rm_cen = np.load("/scratch/r/rbond/lokken/peak-patch-runs/june_18/cens_mhalo_uniform_2e13_3e13_overh.npy")

cenx = rm_cen[:,0]; ceny = rm_cen[:,1]; cenz = rm_cen[:,2]



theta_cen,phi_cen,rshift_cen = xyz_to_ang_z(cenx,ceny,cenz)



minz = 0.1

maxz = 0.8



rm_mocks_all = []



nbins = int((cosmo.comoving_distance(maxz).value - cosmo.comoving_distance(minz).value) // slice_width)

print("Number of distance bins: %d" %nbins)

sat_data1_loaded = False

sat_data2_loaded = False

sat_data3_loaded = False

for i in range(nbins):

dist_slice_min = cosmo.comoving_distance(minz) + slice_width*u.megaparsec*i

dist_slice_max = dist_slice_min + slice_width*u.megaparsec

z_slice_min = z_at_value(cosmo.comoving_distance, dist_slice_min)

z_slice_max = z_at_value(cosmo.comoving_distance, dist_slice_max)

if dist_slice_max.value < 1100:

if sat_data1_loaded == False:

rm_sat = np.load("/scratch/r/rbond/lokken/peak-patch-runs/june_18/sats_mhalo_uniform_2e13_3e13_overh_chi_lt_1032.npy")

satx = rm_sat[:,0]; saty = rm_sat[:,1]; satz = rm_sat[:,2]

theta_sat,phi_sat,rshift_sat = xyz_to_ang_z(satx,saty,satz)

sat_data1_loaded = True

print("Loaded satellite data for chi < 1032.\n")

elif (dist_slice_min.value >= 1032) & (dist_slice_max.value < 2100):

if sat_data2_loaded == False:

rm_sat = np.load("/scratch/r/rbond/lokken/peak-patch-runs/june_18/sats_mhalo_uniform_2e13_3e13_overh_chi_1032_2032.npy")

satx = rm_sat[:,0]; saty = rm_sat[:,1]; satz = rm_sat[:,2]

theta_sat,phi_sat,rshift_sat = xyz_to_ang_z(satx,saty,satz)

sat_data2_loaded = True

print("Loaded satellite data for 1032 < chi < 2032.\n")

print("number of sats in this range:", len(satx))

else:

if sat_data3_loaded == False:

rm_sat = np.load("/scratch/r/rbond/lokken/peak-patch-runs/june_18/sats_mhalo_uniform_2e13_3e13_overh_chi_2032_3032.npy")

sat_data3_loaded = True

satx = rm_sat[:,0]; saty = rm_sat[:,1]; satz = rm_sat[:,2]

theta_sat,phi_sat,rshift_sat = xyz_to_ang_z(satx,saty,satz)

print("Loaded satellite data for 2032 < chi.\n")


# limit my sample to objects which exist in this redshift bin

rm_inbin = np.where(np.logical_and(rm_inbox['zredmagic']>z_slice_min,rm_inbox['zredmagic']<z_slice_max))[0]

nrm_inbox = len(rm_inbin)

scale_up = 41253./box_area

print("{:d} Redmagic galaxies in this box, in redshift slice {:.2} to {:.2}.\nScaling up to sky by multiplying by {:.2}.".format(nrm_inbox, z_slice_min, z_slice_max, scale_up))

nrm_fullsky = scale_up * nrm_inbox

print("Max rshift_sat = ", np.amax(rshift_sat), " and min rshift_sat = ", np.amin(rshift_sat))

sat_inbin = np.where(np.logical_and(rshift_sat>z_slice_min,rshift_sat<z_slice_max))[0]

nsat_inbin = len(sat_inbin)

if frac_cens != 0.:

cen_inbin = np.where(np.logical_and(rshift_cen>z_slice_min,rshift_cen<z_slice_max))[0]

ncen_inbin = len(cen_inbin)

print("Selecting {:d} simulated satellite galaxies from {:d} total.".format(int(nrm_fullsky*(1-frac_cens)), nsat_inbin))

print("Selecting {:d} simulated central galaxies from {:d} total.".format(int(nrm_fullsky*(frac_cens)), ncen_inbin))

idx_choice_sat = np.random.choice(range(nsat_inbin), size=int(nrm_fullsky*(1-frac_cens)))

idx_choice_cen = np.random.choice(range(ncen_inbin), size=int(nrm_fullsky*frac_cens))

rm_mocks_sat = rm_sat[sat_inbin][idx_choice_sat]

rm_mocks_cen = rm_cen[cen_inbin][idx_choice_cen]

rm_mocks_all.append(np.concatenate((rm_mocks_sat, rm_mocks_cen)))



else:

print("Selecting {:d} simulated satellite galaxies from {:d} total.".format(int(nrm_fullsky), nsat_inbin))

idx_choice = np.random.choice(range(nsat_inbin), size=int(nrm_fullsky))

rm_mocks_all.append(rm_sat[sat_inbin][idx_choice])

rm_mocks_all = np.concatenate(rm_mocks_all)

np.save("/scratch/r/rbond/lokken/peak-patch-runs/june_18/rm_mock_galaxies_mhalo_gt_2e13_50percent_cens.npy", rm_mocks_all)







elif mode == 'cmass':

# First things first, sort all halos by mass from smallest to largest

print("Sorting all centrals by halo mass.")

M_sort = np.argsort(mcen)

xcen = (galcat['xpos_cen'].value)[M_sort]

ycen = (galcat['ypos_cen'].value)[M_sort]

zcen = (galcat['zpos_cen'].value)[M_sort]

mcen = mcen[M_sort]

galcat.close()

if test == True:

subset = np.random.choice(np.arange(len(mcen)), size=10000)

mcen = mcen[subset]; xcen = xcen[subset]; ycen=ycen[subset]; zcen = zcen[subset]

chi = np.sqrt(xcen**2+ycen**2+zcen**2)

logM = np.log10(mcen)



# Prepare subsample of centrals for the mock CMASS gal selection

mean_mass = 12.78 # mean log(M) of halos in CMASS sample

sigma_mass = 0.37



# cut the halo catalog to logM range that's similar to CMASS sample 

mean_minus_minmass = mean_mass - np.amin(logM)

first_cut = (logM<mean_mass+mean_minus_minmass)

print("Doing first mass cut.")

M_firstcut = mcen[first_cut]

x_firstcut = xcen[first_cut]

y_firstcut = ycen[first_cut]

z_firstcut = zcen[first_cut]

chi_firstcut = chi[first_cut]

logM_firstcut = logM[first_cut]

indices = np.arange(len(M_firstcut))

# the following is a list of the indices of the mass and position arrays that I want to keep after a mass cut

indices_keep = []

# make a Gaussian probability density function centered around the mean

Nhalo = 100

stepsize = (np.amax(logM_firstcut)-np.amin(logM_firstcut))/Nhalo

print("Making a Gaussian probability density function centered around logM = 12.78, with stepsize logM = {:.2}.".format(stepsize))

normal_masses = np.arange(np.amin(logM_firstcut),np.amax(logM_firstcut),stepsize)

gaussian_logM = norm.pdf(normal_masses, mean_mass, sigma_mass)

for i in range(len(normal_masses)-1):

n = int(stepsize * gaussian_logM[i] * 10**7) # number of masses to include in this range

# randomly sample the actual halo sample in that range

indices_in_range = indices[(logM_firstcut<normal_masses[i+1])&(logM_firstcut>normal_masses[i])]

idx_keep = np.random.choice(indices_in_range,size=n)

for index in idx_keep:

indices_keep.append(index)



# update halo properties again

M_mock_cm = M_firstcut[indices_keep]

x_mock_cm = x_firstcut[indices_keep]

y_mock_cm = y_firstcut[indices_keep]

z_mock_cm = z_firstcut[indices_keep]

chi_mock_cm = chi_firstcut[indices_keep]



# Plot this mass distribution



# plt.hist(np.log10(M_cmass), bins=500)

# plt.title("Mock CMASS Mhalo dist, all z")

# plt.xlabel(r'log(M/$M_\odot$)')

# plt.ylabel('Counts')

# plt.savefig("/scratch/r/rbond/lokken/peak-patch-runs/june_18/cmass_mock_centrals_logMh_hist.pdf")



# Plot the current redshift distribution

# plt.hist(chi_mock_cm, bins=100)

# plt.title("Mock CMASS gal redshift dist")

# plt.xlabel("$\chi$")

# plt.ylabel("Counts")

# plt.savefig("/scratch/r/rbond/lokken/peak-patch-runs/june_18/cmass_mock_centrals_chi_hist.pdf")





# convert to mass, comoving distance, radial velocity, redshfit, RA and DEc



with fits.open("/scratch/r/rbond/lokken/data/large_region_galaxy_DR12v5_CMASS_South.fits") as hdu:

cmdat = hdu[1].data

ra_min = 0.; ra_max = 7.; dec_min = -3.; dec_max = 4.



# select all cmass within a box, and make sure this is filled at redshifts out to 0.8

in_ra = (cmdat['ra'] > ra_min) & (cmdat['ra'] < ra_max)

in_dec = (cmdat['dec'] > dec_min) & (cmdat['dec'] < dec_max)

box_area = (ra_max - ra_min) * (dec_max - dec_min)



cm_inbox = cmdat[in_ra&in_dec]

theta_mock_cm,phi_mock_cm,rshift_mock_cm = xyz_to_ang_z(x_mock_cm,y_mock_cm,z_mock_cm)



minz = 0.1

maxz = 0.8



cm_mocks_all = []



nbins = int((cosmo.comoving_distance(maxz).value - cosmo.comoving_distance(minz).value) // slice_width)

print("Number of distance bins: %d" %nbins)

for i in range(nbins):

dist_slice_min = cosmo.comoving_distance(minz) + slice_width*u.megaparsec*i

dist_slice_max = dist_slice_min + slice_width*u.megaparsec

z_slice_min = z_at_value(cosmo.comoving_distance, dist_slice_min)

z_slice_max = z_at_value(cosmo.comoving_distance, dist_slice_max)

# limit my sample to objects which exist in this redshift bin

cm_inbin = np.where(np.logical_and(cm_inbox['Z']>z_slice_min,cm_inbox['Z']<z_slice_max))[0]

ncm_inbox = len(cm_inbin)

scale_up = 41253./box_area

print("{:d} CMASS galaxies in this box, in redshift slice {:.2} to {:.2}.\nScaling up to full sky by multiplying by {:.2}.".format(ncm_inbox, z_slice_min, z_slice_max, scale_up))



ncm_fullsky = int(scale_up * ncm_inbox)

mock_cm_inbin = np.where(np.logical_and(rshift_mock_cm>z_slice_min,rshift_mock_cm<z_slice_max))[0]

n_mocks_inbin = len(mock_cm_inbin)

print("Selecting {:d} simulated central galaxies from {:d} total.".format(ncm_fullsky, n_mocks_inbin))

idx_choice = np.random.choice(range(n_mocks_inbin), size=ncm_fullsky)

mocks_inbin = np.zeros((ncm_fullsky, 4))

mocks_inbin[:,0] = x_mock_cm[mock_cm_inbin][idx_choice]

mocks_inbin[:,1] = y_mock_cm[mock_cm_inbin][idx_choice]

mocks_inbin[:,2] = z_mock_cm[mock_cm_inbin][idx_choice]

mocks_inbin[:,3] = M_mock_cm[mock_cm_inbin][idx_choice]

cm_mocks_all.append(mocks_inbin)



cm_mocks_all = np.concatenate(cm_mocks_all)

np.save('/scratch/r/rbond/lokken/peak-patch-runs/june_18/cm_mock_galaxies_mhalo_cmdist_allcens.npy', cm_mocks_all)

