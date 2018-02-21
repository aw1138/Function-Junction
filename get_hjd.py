#!/usr/bin/python

import sys
import numpy as np
import glob
import pandas as pd
from astropy.io import fits
from astropy import time, coordinates as coord, units as u
from astropy.units import imperial
from astropy.units import cds

# input:
# path to observation tables (b, r)
# path to observation table (if no second table, make sure save
# 	does not currentl exist)
# name of table to save w/ hjd and mjd

try:
	
	# grab every file with the pattern specified by the user
	files_blue = glob.glob(sys.argv[1])
	files_red = glob.glob(sys.argv[2])

	# create some empty dictionaries to store information
	bobs = {}
	bobs['dates_blue'] = []
	bobs['mjd_blue'] = []

	robs = {}
	robs['dates_red'] = []
	robs['mjd_red'] = [] 

	# for every file specified by user:
	for path1, path2 in zip(files_blue, files_red):
	    
	    # open the file
	    table1 = fits.open(path1)
	    table2 = fits.open(path2)
	    
	    # get the actual utc date for the observations from the header
	    bobs['dates_blue'].append(table1[0].header['DATE-OBS'])
	    robs['dates_red'].append(table2[0].header['DATE-OBS'])
	    
	    # get the mjd for the observations from the header
	    bobs['mjd_blue'].append(table1[0].header['MJD-OBS'])
	    robs['mjd_red'].append(table2[0].header['MJD-OBS'])
	    
	    # close the file
	    table1.close()
	    table2.close()

	bobs = pd.DataFrame(bobs) # blue observation dates, mjd
	robs = pd.DataFrame(robs) # red observation dates, mjd

	# get the mjd for dates where there is both blue and red obs

	dates = pd.merge(bobs, robs, how='inner', left_on=['dates_blue'], right_on=['dates_red'])
	dates = dates[['dates_blue', 'dates_red', 'mjd_blue', 'mjd_red']]

	# get which mjd is greater for every observation
	big_blue = dates['mjd_blue'] > dates['mjd_red']
	big_red = dates['mjd_red'] > dates['mjd_blue']

	# where blue is later obs take that mjd
	mjd = dates['mjd_blue'].where(big_blue).rename(index=str, columns={'mjd_blue':'mjd'})
	utc = dates['dates_blue'].where(big_blue).rename(index=str, columns={'dates_blue':'utc'})

	# where red is later obs take that mjd
	temp = dates['mjd_red'].where(big_red).rename(index=str, columns={'mjd_red':'mjd'}).dropna()
	temp_u = dates['dates_red'].where(big_red).rename(index=str, columns={'dates_red':'utc'}).dropna()
	#replace missing mjd slots in blue list with the ones in the red
	#this gives us the bigger mjd for every date in a table
	mjd = mjd.fillna(temp)
	utc = utc.fillna(temp_u)

	# getting mjd for dates where there is only one observation

	dates_not = pd.merge(bobs, robs, how='outer', left_on=['dates_blue'], 
	                     right_on=['dates_red'], left_index=True, indicator=True)
	dates_not = dates_not.query('_merge != "both"')

	# rename the columns so we can create one merged table
	mjd_not = dates_not['mjd_blue'].rename(index=str, columns={'mjd_blue':'mjd'})
	utc_not = dates_not['dates_blue'].rename(index=str, columns={'dates_blue':'utc'})

	temp_not = dates_not['mjd_red'].rename(index=str, columns={'mjd_red':'mjd'}).dropna()
	utc_temp = dates_not['dates_red'].rename(index=str, columns={'dates_red':'utc'}).dropna()

	mjd_not = mjd_not.fillna(temp_not).reset_index(drop=True) + (30 * u.min).to(cds.MJD).value
	utc_not = utc_not.fillna(utc_temp).reset_index(drop=True)

	# combine single obs date table with the double obs date table
	mjd = mjd.append(mjd_not).sort_values().reset_index(drop=True)
	utc = utc.append(utc_not).sort_values().reset_index(drop=True)

	# get the hjd times from the mjd times

	# specify the star location in the sky
	v367cyg = coord.SkyCoord("20:47:59.5849", "+39:17:15.723", unit=(u.hourangle, u.deg), frame='icrs')
	# specify observatory location
	loc = coord.EarthLocation.from_geodetic(43.0776*u.degree, 89.6717*u.degree, height=1188*imperial.ft)
	# turn mjd table into an astropy time object
	times = time.Time(mjd, format='mjd', scale='utc', location=loc)
	# astropy calcs the difference between hjd and mjd
	ltt_helio = times.light_travel_time(v367cyg, 'heliocentric')
	# add difference to mjd values to get hjd 
	hjd = mjd + ltt_helio.value

	# put mjd and hjd into one, shared table
	jd = {}
	jd['mjd'] = mjd
	jd['hjd'] = hjd
	jd['utc'] = utc
	jd = pd.DataFrame(jd)

	# save the full table
	jd.to_csv(sys.argv[3])


except ValueError:

		# pull all red observation files
	files = glob.glob(sys.argv[1])

	obs = {}
	obs['dates'] = []
	obs['mjd'] = [] 

		# for every file specified by user:
	for path in files:
	    
	    # open the file
	    table = fits.open(path)
	    
	    # get the actual utc date for the observations from the header
	    obs['dates'].append(table[0].header['DATE-OBS'])
	    
	    # get the mjd for the observations from the header
	    obs['mjd'].append(table[0].header['MJD-OBS'])
	    
	    # close the file
	    table.close()

	obs = pd.DataFrame(obs) # red observation dates, mjd

	mjd = obs['mjd'] + (30 * u.min).to(cds.MJD).value
	utc = obs['dates']

	
	# specify the star location in the sky
	v367cyg = coord.SkyCoord("20:47:59.5849", "+39:17:15.723", unit=(u.hourangle, u.deg), frame='icrs')
	# specify observatory location
	loc = coord.EarthLocation.from_geodetic(43.0776*u.degree, 89.6717*u.degree, height=1188*imperial.ft)
	# turn mjd table into an astropy time object
	times = time.Time(mjd, format='mjd', scale='utc', location=loc)
	# astropy calcs the difference between hjd and mjd
	ltt_helio = times.light_travel_time(v367cyg, 'heliocentric')
	# add difference to mjd values to get hjd 
	hjd = mjd + ltt_helio.value

	# put mjd and hjd into one, shared table
	jd = {}
	jd['mjd'] = mjd
	jd['hjd'] = hjd
	jd['utc'] = utc
	jd = pd.DataFrame(jd)

	# save the full table
	jd.to_csv(sys.argv[2])


# else:

