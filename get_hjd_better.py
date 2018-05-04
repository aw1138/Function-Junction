#!/usr/bin/env python

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

## complete except for adding in flags

# define a function that will grab times from fits table headers
def get_times(path, datekey='DATE-OBS', mjdkey='MJD-OBS'):

	# open the file
    table = fits.open(path)
    
    # get the actual utc date for the observations from the header
    date = table[0].header[datekey]
    
    # get the mjd for the observations from the header
    mjd = table[0].header[mjdkey]
    
    # close the file
    table.close()

    return date, mjd



try:
	
	# grab every file with the pattern specified by the user
	files_blue = glob.glob(sys.argv[1])
	files_red = glob.glob(sys.argv[2])

	# create some empty dictionaries to store information
	bobs = {}
	robs = {}

	# for every file in path given, get utc dates and mjd

	bobs['date'] = [i for i,j in list(map(get_times, files_blue))]
	bobs['mjd_blue'] = [j for i,j in list(map(get_times, files_blue))]

	robs['date'] = [i for i,j in list(map(get_times, files_red))]
	robs['mjd_red'] = [j for i,j in list(map(get_times, files_red))]


	bobs = pd.DataFrame(bobs) # blue observation dates, mjd
	robs = pd.DataFrame(robs) # red observation dates, mjd

	# get the mjd for dates where there is both blue and red obs

	dates = pd.merge(bobs, robs, how='inner', on=['date'])

	# get which mjd is greater for every observation
	big_blue = dates['mjd_blue'] > dates['mjd_red']
	big_red = dates['mjd_red'] > dates['mjd_blue']

	# where blue is later obs take that mjd
	mjd = dates['mjd_blue'].where(big_blue).rename(index=str, columns={'mjd_blue':'mjd'})
	utc = dates['date'].where(big_blue)

	# where red is later obs take that mjd
	temp = dates['mjd_red'].where(big_red).rename(index=str, columns={'mjd_red':'mjd'}).dropna()
	temp_u = dates['date'].where(big_red).dropna()
	#replace missing mjd slots in blue list with the ones in the red
	#this gives us the bigger mjd for every date in a table
	mjd = mjd.fillna(temp)
	utc = utc.fillna(temp_u)

	# getting mjd for dates where there is only one observation

	dates_not = pd.merge(bobs, robs, how='outer', on=['date'], 
							left_index=True, indicator=True)
	dates_not = dates_not.query('_merge != "both"')

	# rename the columns so we can create one merged table
	mjd_not = dates_not['mjd_blue'].rename(index=str, columns={'mjd_blue':'mjd'})
	utc_not = dates_not['date']

	temp_not = dates_not['mjd_red'].rename(index=str, columns={'mjd_red':'mjd'}).dropna()
	utc_temp = dates_not['date'].dropna()

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
	jd['date'] = utc
	jd = pd.DataFrame(jd)

	# save the full table
	jd.to_csv(sys.argv[3])


except ValueError:

	# this could use some fixing to make sure not an issue

	# pull all red observation files
	files = glob.glob(sys.argv[1])

	obs = {}

	# for every file in path get the utc dates and mjd
	obs['dates'] = [i for i,j in list(map(get_times, files_blue))]
	obs['mjd'] = [j for i,j in list(map(get_times, files_blue))]

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


else:

	import traceback

	traceback.print_exc()