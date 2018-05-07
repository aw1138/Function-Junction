#!/usr/bin/env python

import sys
import numpy as np
import glob
import pandas as pd
from astropy.io import fits
from astropy import units as u
import itertools
from astropy.time import Time

# quick decription of inputs
# blue file path, red file path, save name

# import blue filters
b_fil = pd.read_table('../filters/B.fil', sep='  ', names=['wavelength', 'weight'], engine='python')
b_fil.index.name = 'B_Band'
ux_fil = pd.read_table('../filters/UX.fil', sep='  ', names=['wavelength', 'weight'], engine='python')
ux_fil.index.name = 'U_Band'
# bjp_fil = pd.read_table('../filters/BJ+.fil', sep='  ', names=['wavelength', 'weight'], engine='python')
# bjm_fil = pd.read_table('../filters/BJ-.fil', sep='  ', names=['wavelength', 'weight'], engine='python')

#import red filters
i_fil = pd.read_table('../filters/I.fil', sep='  ', names=['wavelength', 'weight'], engine='python')
i_fil.index.name = 'I_Band'
r_fil = pd.read_table('../filters/R.fil', sep='  ', names=['wavelength', 'weight'], engine='python')
r_fil.index.name = 'R_Band'

#import middle filter
v_fil = pd.read_table('../filters/V.fil', sep='  ', names=['wavelength', 'weight'], engine='python')
v_fil.index.name = 'V_Band'

# import a table of the system errors
extra_err = pd.read_table('HPOL_Sys_Err_Aislynn.txt', sep=' ') 

# create lists of all the files matching the pattern provided
files_blue = glob.glob(sys.argv[1])
files_red = glob.glob(sys.argv[2])

def polarimetry_calc(tup):
	
	"""
	A function that will calculate the angle and amount of polarized light
	given a false filter and a path to spectra data.

	Parameters:
	-----------

	tup: tuple, list of tuples
		A tuple or list of tuples of the form (filter, path). These are 
		used to obtain all relevant information that's needed.
	"""
	fil, path = tup
	# tuple to make the map function work later
	
	# load in the table and its data
	table = fits.open(path)
	table_data = table[1].data

	# get the filter-wavelength weight
	weights = np.interp(table_data[0]['wavelength'], fil['wavelength'],
							 fil['weight'], left=0, right=0)
	# define the flux
	flux = table[0].data

	# get the flux weight for q
	qweight = table_data[0]['q'] * flux * weights
	# sum/integrate it
	qsum = np.trapz(qweight, table_data[0]['wavelength'])

	# get the weight for just the flux
	fweight = flux * weights
	# sum it
	fsum = np.trapz(fweight, table_data[0]['wavelength'])

	# find q
	polq = qsum / fsum

	# get the flux weight for u
	uweight = table_data[0]['u'] * flux * weights
	# sum it
	usum = np.trapz(uweight, table_data[0]['wavelength'])

	# find u
	polu = usum / fsum

	# find the total polarized light
	pol_tot = np.sqrt(polq**2 + polu**2)

	# find the position angle of the light, convert to deg
	pos_ang = 0.5 * np.arctan2(polu, polq) * u.rad
	pos_ang = pos_ang.to(u.deg).value

	# grab the date to sort by later
	date = table[0].header['DATE-OBS']

	# define the error from the og table
	err = table_data[0]['error']

	# words
	s = pd.DataFrame(weights * flux * err, columns=['pol_err'])
	s = s.loc[(s!=0).any(axis=1)]

	# get poisson error
	fish = np.sum(s) / np.sum(flux * weights) / np.sqrt(s.size)

	# turn the start and end times for sys errs into Time objs
	start = Time(extra_err['Start_Date'].values.astype(str), scale='utc')
	end = Time(extra_err['End_Date'].values.astype(str), scale='utc')

	# turn the date of obs into a time object and find which err it falls into
	time = Time(date, scale='utc')
	errmatch = (time >= start) & (time <= end)

	if all(errmatch) is False:
		t = pd.Series(start-time).append(pd.Series(end-time))
		errmatch = np.abs(t).values.argmin()
		if errmatch > 19:
			errmatch -= 19
	# get the err for the polization
	errpol = np.sqrt(fish**2 + extra_err[fil.index.name][errmatch]**2)

	# get the error for the position angle
	errang = 90 / np.pi * errpol / pol_tot

	table.close()

	return pol_tot, pos_ang, date, errpol.values[0], errang.values[0]


def det_vdates(path, datekey='DATE-OBS'):
	"""
	A function to be used to obtain the date associated with a fits file.
	Meant to be used in conjunction with the map built-in function.

	Parameters:
	-----------

	path: string
		The path to the files to get dates from
	datekey: string
		The key in the header of the file for the UTC or other needed date.
	"""

	table = fits.open(path)
	
	# grab the date of the obs and the path to the data
	tbl = {}
	tbl['date'] = [table[0].header[datekey]]
	tbl['path'] = [path]
	
	# turn them into a DataFrame
	tbl = pd.DataFrame(tbl)

	# close the file to conserve memory
	table.close()
	
	return tbl


def v_band_pol(tup, fil=v_fil):
	"""
	A function to do polarimetry in only the V band.

	Parameters:
	-----------

	tup: tuple of strings
		A tuple of the form (red files path, blue files path).

	fil: dictionary or DataFrame
		The filter the calculations will be done for. Default is the V filter.
	"""

	# split the tuple into the paths
	path_r, path_b = tup
	
	# load in the red table and its data
	table_r = fits.open(path_r)
	data_r = table_r[1].data
	# grab the data needed, put into its own dataframe
	r_hold = {'wavelength': data_r[0]['wavelength'], 'q': data_r[0]['q'] , 
			  'u': data_r[0]['u'] , 'error': data_r[0]['error']}
	data_r = pd.DataFrame(r_hold)

	# load in the blue table and its data
	table_b = fits.open(path_b)
	data_b = table_b[1].data
	# grab the data needed, put into its own dataframe
	b_hold = {'wavelength': data_b[0]['wavelength'], 'q': data_b[0]['q'] , 
			  'u': data_b[0]['u'] , 'error': data_b[0]['error']}
	data_b = pd.DataFrame(b_hold)

	# get the fluxes
	flux_r = table_r[0].data
	flux_b = table_b[0].data
	# normalize one of the fluxes to match the other
	flux_b = flux_b * np.mean(flux_r) / np.mean(flux_b)

	# make the fluxes their own columns in the table
	data_r['flux'] = flux_r
	data_b['flux'] = flux_b

	# concatinate the data
	# cannot do byteswap solution to merge bc will fuck up data
	data = pd.concat([data_b, data_r], join='outer')
	
	# get the filter wavelength weights
	vweights = np.interp(data['wavelength'], v_fil['wavelength'], 
							v_fil['weight'], left=0, right=0)

	# get the weight for q
	qweight = data['q'] * data['flux'] * vweights
	# sum them
	qsum = np.trapz(qweight, data['wavelength'])

	# get the flux weight
	fweight = data['flux'] * vweights
	# sum them
	fsum = np.trapz(fweight, data['wavelength'])

	# get total q
	polq = qsum / fsum

	# get u weight
	uweight = data['u'] * data['flux'] * vweights
	# sum them
	usum = np.trapz(uweight, data['wavelength'])

	# get total u
	polu = usum / fsum

	# get total polarized light
	pol_tot = np.sqrt(polq**2 + polu**2)

	# get the position angle
	pos_ang = 0.5 * np.arctan2(polu, polq) * u.rad
	pos_ang = pos_ang.to(u.deg).value
	
	# get the date for later
	date = table_r[0].header['DATE-OBS']

	# define the error from the og table
	err = data['error']
	flux = data['flux']
	# words
	s = pd.DataFrame(vweights * flux * err, columns=['pol_err'])
	s = s.loc[(s!=0).any(axis=1)]

	# get poisson error
	fish = np.sum(s) / np.sum(flux * vweights) / np.sqrt(s.size)

	# turn the start and end times for sys errs into Time objs
	start = Time(extra_err['Start_Date'].values.astype(str), scale='utc')
	end = Time(extra_err['End_Date'].values.astype(str), scale='utc')

	# turn the date of obs into a time object and find which err it falls into
	time = Time(date, scale='utc')
	errmatch = (time >= start) & (time <= end)

	if all(errmatch) is False:
		t = pd.Series(start-time).append(pd.Series(end-time))
		errmatch = np.abs(t).values.argmin()
		if errmatch > 19:
			errmatch -= 19

	# get the err for the polarization
	errpol = np.sqrt(fish**2 + extra_err['V_Band'][errmatch]**2)

	# get the error for the position angle
	errang = 90 / np.pi * errpol / pol_tot
	
	table_r.close()
	table_b.close()

	return pol_tot, pos_ang, date, errpol.values[0], errang.values[0]

######################################################################
#################### Acually creating the tables #####################
######################################################################

## get red filter mags

# create empty placeholder dictionary
rmags = {}

# create a list of tuples with (i_fil, filename)
prod_i = itertools.product([i_fil], files_red)
# map the function to the list of tuples
rmags['Ipol'] = [i for i,j,k,l,m in list(map(polarimetry_calc, prod_i))]

# for some reason you have to define it again or it doesn't work
prod_i = itertools.product([i_fil], files_red)
rmags['Iangle'] = [j for i,j,k,l,m in list(map(polarimetry_calc, prod_i))]

# create a list of tuples with (i_fil, filename)
prod_i = itertools.product([i_fil], files_red)
# map the function to the list of tuples
rmags['Ipol_err'] = [l for i,j,k,l,m in list(map(polarimetry_calc, prod_i))]

# for some reason you have to define it again or it doesn't work
prod_i = itertools.product([i_fil], files_red)
rmags['Iangle_err'] = [m for i,j,k,l,m in list(map(polarimetry_calc, prod_i))]

# create a list of tuples with (r_fil, filename)
prod_r = itertools.product([r_fil], files_red)
rmags['Rpol'] = [i for i,j,k,l,m in list(map(polarimetry_calc, prod_r))]

prod_r = itertools.product([r_fil], files_red)
rmags['Rangle'] = [j for i,j,k,l,m in list(map(polarimetry_calc, prod_r))]

# create a list of tuples with (r_fil, filename)
prod_r = itertools.product([r_fil], files_red)
rmags['Rpol_err'] = [l for i,j,k,l,m in list(map(polarimetry_calc, prod_r))]

prod_r = itertools.product([r_fil], files_red)
rmags['Rangle_err'] = [m for i,j,k,l,m in list(map(polarimetry_calc, prod_r))]

prod_r = itertools.product([r_fil], files_red)
rmags['date'] = [k for i,j,k,l,m in list(map(polarimetry_calc, prod_r))]

# convert to DataFrame
rmags = pd.DataFrame(rmags)

## get blue filter mags

# create empty placeholder dictionary
bmags = {}

# create a list of tuples with (b_fil, filename)
prod_b = itertools.product([b_fil], files_blue)
bmags['Bpol'] = [i for i,j,k,l,m in list(map(polarimetry_calc, prod_b))]

# for some reason you have to define it again or it doesn't work
prod_b = itertools.product([b_fil], files_blue)
bmags['Bangle'] = [j for i,j,k,l,m in list(map(polarimetry_calc, prod_b))]

prod_b = itertools.product([b_fil], files_blue)
bmags['Bpol_err'] = [l for i,j,k,l,m in list(map(polarimetry_calc, prod_b))]

prod_b = itertools.product([b_fil], files_blue)
bmags['Bangle_err'] = [m for i,j,k,l,m in list(map(polarimetry_calc, prod_b))]

# create a list of tuples with (u_fil, filename)
prod_u = itertools.product([ux_fil], files_blue)
bmags['Upol'] = [i for i,j,k,l,m in list(map(polarimetry_calc, prod_u))]

prod_u = itertools.product([ux_fil], files_blue)
bmags['Uangle'] = [j for i,j,k,l,m in list(map(polarimetry_calc, prod_u))]

prod_u = itertools.product([ux_fil], files_blue)
bmags['Upol_err'] = [l for i,j,k,l,m in list(map(polarimetry_calc, prod_u))]

prod_u = itertools.product([ux_fil], files_blue)
bmags['Uangle_err'] = [m for i,j,k,l,m in list(map(polarimetry_calc, prod_u))]

# create a list of tuples with (bjp_fil, filename)
# prod_bjp = itertools.product([bjp_fil], files_blue)
# bmags['BJ+pol'] = [i for i,j,k,l,m in list(map(polarimetry_calc, prod_bjp))]

# prod_bjp = itertools.product([bjp_fil], files_blue)
# bmags['BJ+angle'] = [j for i,j,k,l,m in list(map(polarimetry_calc, prod_bjp))]

# # create a list of tuples with (bjp_fil, filename)
# prod_bjm = itertools.product([bjm_fil], files_blue)
# bmags['BJ-pol'] = [i for i,j,k,l,m in list(map(polarimetry_calc, prod_bjm))]

# prod_bjm = itertools.product([bjm_fil], files_blue)
# bmags['BJ-angle'] = [j for i,j,k,l,m in list(map(polarimetry_calc, prod_bjm))]

prod_b = itertools.product([b_fil], files_blue)
bmags['date'] = [k for i,j,k,l,m in list(map(polarimetry_calc, prod_b))]

# turn the placeholder into a pandas DataFrame
bmags = pd.DataFrame(bmags)

## get the v filter mags

# get date and path for blue and red separately
bobs = pd.concat(list(map(det_vdates, files_blue)))
robs = pd.concat(list(map(det_vdates, files_red)))

# combine, take only dates where both obs happened
v_dates = pd.merge(bobs, robs, how='inner', on=['date'], suffixes=['_b','_r'])

# zip the paths for those dates together into a tuple
file_tup = list(zip(v_dates['path_r'], v_dates['path_b']))

# create a placeholder table
vmags = {}

# fill it with the values from the function
vmags['Vpol'] = [i for i,j,k,l,m in list(map(v_band_pol, file_tup))]
vmags['Vangle'] = [j for i,j,k,l,m in list(map(v_band_pol, file_tup))]
vmags['date'] = [k for i,j,k,l,m in list(map(v_band_pol, file_tup))]
vmags['Vpol_err'] = [l for i,j,k,l,m in list(map(v_band_pol, file_tup))]
vmags['Vangle_err'] = [m for i,j,k,l,m in list(map(v_band_pol, file_tup))]

# convert placeholder to a DataFrame
vmags = pd.DataFrame(vmags)

# merge r and b tables
final = pd.merge(rmags, bmags, how='outer', on=['date'])
final = final.merge(vmags, how='outer', on=['date'])
final = final.sort_values(by=['date']).reset_index(drop=True)

# save the merged table
final.to_csv(sys.argv[3])