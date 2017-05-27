import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from astropy.table import Table

def getindex(d1, d2, r):
    """
    A function, created by (name: not me), to obtain the index of the star that matches those put into it.
    """
    t = cKDTree(d2)
    d, idx = t.query(d1, distance_upper_bound=r)
    return idx


def matchLists(tol, ra1, dec1, ra2, dec2, filename):
    """
    A function that matches the indexes of objects in one set of data to the same objects' indexes in a second set of data. Main credit: (not me). Edited/revised by A. Wallach.

    Attributes:
	tol: float
	    The tolerance within which the match must be. How close the RA and DEC must be between sets in order to be matched.
	ra1: numpy array
	    RA of the objects in the first set of data.
	dec1: numpy array
	    Declination of the objects in the first set of data.
	ra2: numpy array
	     RA of the objects in the second set of data.
	dec2: numpy array
	    Declination of the objects in the second set of data.
	filename: string
	    The name of the file that will be saved as a csv file.

    This function takes a name for the final file, tolerance, RAs, and Declinations, and creates a pandas DataFrame object (table) containing the index of the objects in the first set of data as index1 and the matches to those objects as index2. If there is no match, the function returns an index value of -1. The table is then written as a .csv file, and imported back in using pandas.
    """
    d1,d2 = np.empty((ra1.size, 2)), np.empty((ra2.size, 2))
    d1[:,0],d1[:,1],d2[:,0],d2[:,1] = ra1,dec1,ra2,dec2
    out2 = getindex(d1,d2,tol)
    out2[out2==ra2.size] = -1
    
    out1 = np.arange(ra1.size)
    table = [out1, out2]
    nms = ('index1','index2')
    fmt = {'index1':'%10d', 'index2':'%10d'}
    t = Table(table, names=nms)
    t.to_pandas().to_csv(filename+'.csv', index=False)
    final = pd.read_csv(filename+'.csv')
    return final


def find_diff_mag(magnitude1_1, magnitude2_1, idx1, idx2, 
                   magnitude1_2=None, magnitude2_2=None,
                   magnitude1_3=None, magnitude2_3=None,
                   right_ascen1=None, right_ascen2=None,
                   declination1=None, declination2=None):
	"""
	A function to find the difference in magnitude between two observation epics.
    
	Attributes:
		magnitude1_1: list or numpy array
			magnitude in epic one in the first chosen band
		magnitude2_1: list or numpy array
			magnitude in epic two in the first chosen band
		idx1: list or numpy array
			index of the stars in epic one
		idx2: list or numpy array
			indices of the stars in epic two that match the stars in idx1
		magnitude1_2: list or numpy array
			magnitude in epic one in a second chosen band
		magnitude2_2: list or numpy array
			magnitude in epic two in a second chosen band
		magnitude1_3: list or numpy array
			magnitude in epic one in a third chosen band
		magnitude2_3: list or numpy array
			magnitude in epic two in a third chosen band
		right_ascen1: list or numpy array
			ra of the stars in epic one
		right_ascen2: list or numpy array
			ra of the stars in epic two
		declination1: list or numpy array
			dec of the stars in epic one
		declination2: list or numpy array
			dec of the stars in epic two
    
	This function creates an table with columns made from the input data and the difference between each star's magnitude(from one epic to the next) in the bands given as well as the difference in ra and dec (from one epic to the next). It returns the newly created table as a pandas DataFrame object.

	This function can take the magnitudes of up to three different bands of light. Only one is required for the function to run.
	"""
	## create emtpy dictionary with columns for the data
	diffapexmags = {}
	diffapexmags['diffmag_1'] = []
	diffapexmags['mag1_1'] = []
	diffapexmags['mag2_1'] = []
	diffapexmags['diffmag_2'] = []
	diffapexmags['mag1_2'] = []
	diffapexmags['mag2_2'] = []
	diffapexmags['diffmag_3'] = []
	diffapexmags['mag1_3'] = []
	diffapexmags['mag2_3'] = []
	diffapexmags['diff_ra'] = []
	diffapexmags['ra1'] = []
	diffapexmags['ra2'] = []
	diffapexmags['diff_dec'] = []
	diffapexmags['dec1'] = []
	diffapexmags['dec2'] = []
	for idx1, idx2 in zip(idx1, idx2):
	## remove non-matched stars
		if idx2 == -1:
			diffmag_1 = np.nan
			mag1_1 = np.nan
			mag2_1 = np.nan
			diffmag_2 = np.nan
			mag1_2 = np.nan
			mag2_2 = np.nan
			diffmag_3 = np.nan
			mag1_3 = np.nan
			mag2_3 = np.nan
			ra1 = np.nan
			ra2 = np.nan
			dec1 = np.nan
			dec2 = np.nan
	## find the differences, fill all columns
		if idx2 != -1:
			mag1_1 = magnitude1_1[idx1]
			mag2_1 = magnitude2_1[idx2]
			diffmag_1 = mag1_1 - mag2_1
			diffapexmags['diffmag_1'].append(diffmag_1)
			diffapexmags['mag1_1'].append(mag1_1)
			diffapexmags['mag2_1'].append(mag2_1)

			if magnitude1_2 == None:
				mag1_2 = np.nan
				mag2_2 = np.nan
				diffmag_2 = np.nan
			else:
				mag1_2 = magnitude1_2[idx1]
				mag2_2 = magnitude2_2[idx2]
				diffmag_2 = mag1_2 - mag2_2
			diffapexmags['diffmag_2'].append(diffmag_2)
			diffapexmags['mag1_2'].append(mag1_2)
			diffapexmags['mag2_2'].append(mag2_2)
		    
			if magnitude1_3 == None:
				mag1_3 = np.nan
				mag2_3 = np.nan
				diffmag_3 = np.nan
			else:
				mag1_3 = np.nan
				mag2_3 = np.nan
				diffmag_3 = np.nan
			diffapexmags['diffmag_3'].append(diffmag_3)
			diffapexmags['mag1_3'].append(mag1_3)
			diffapexmags['mag2_3'].append(mag2_3)

			if right_ascen1 == None:
				ra1 = np.nan
				ra2 = np.nan
				diff_ra = np.nan
			else:
				ra1 = right_ascen1[idx1]
				ra2 = right_ascen2[idx2]
				diff_ra = ra1 - ra2
			diffapexmags['diff_ra'].append(diff_ra)
			diffapexmags['ra1'].append(ra1)
			diffapexmags['ra2'].append(ra2)

			if declination1 == None:
				dec1 = np.nan
				dec2 = np.nan
				diff_dec = np.nan
			else:
				dec1 = declination1[idx1]
				dec2 = declination2[idx2]
				diff_dec = dec1 - dec2
			diffapexmags['diff_dec'].append(diff_dec)
			diffapexmags['dec1'].append(dec1)
			diffapexmags['dec2'].append(dec2)

	# convert dictonary into pandas DataFrame object
	tablever = pd.DataFrame.from_dict(diffapexmags)
	# return the DataFrame
	return tablever

