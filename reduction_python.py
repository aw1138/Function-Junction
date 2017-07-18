def reduce_data(data_path, flat_path, dark_path, bias_path, 		reduced_image_folder, over_write=False):
	
	"""
	A function to simply reduce astronomical data in Python. Needed 
	packages: matplotlib, numpy, astropy, scipy, glob, and os.

	This function currently has no way to differentiate images with 
	different exposure lengths from the darks and/or different filters 
	from the flats. MAKE SURE YOU PUT IN THE CORRECT IMAGES.

	Attributes
		
		data_path: string
			The path to your science images.
		flat_path: string
			The path to your flat images for this filter.
		dark_path: string
			The path to your dark images for this exposure 
			length.
		bias_path: string 
			The path to your bias images.
		reduced_image_folder: string
			How you want your reduced images saved.
		over_write: boolean, optional
			Whether you want the files saved to the reduced image
			folder overwritten by the system. The default is
			False.
	"""

	import matplotlib.pyplot as plt
	import numpy as np
	from astropy.io import fits
	from matplotlib.colors import LogNorm
	import glob
	import os.path
	from scipy import stats

	flat_data = []
	for i,v in enumerate(glob.glob(flat_path)):
	## Idk how to get the number of items in this path
	## Is there a way to do it without a for loop?
		flat = fits.getdata(v)
		flat = np.float64(flat)
		flat_data.append(flat)
	flat_data = np.array(flat_data)

	dark_data = []
	for i,v in enumerate(glob.glob(dark_path)):
		dark = fits.getdata(v)
		dark = np.float64(dark)
		dark_data.append(dark)
	dark_data = np.array(dark_data)

	bias_data = []
	for i,v in enumerate(glob.glob(bias_path)):
		bias = fits.getdata(v)
		bias = np.float64(bias)
		bias_data.append(bias)
	bias_data = np.array(bias_data)

	image_data = []
	for i,v in enumerate(glob.glob(data_path)):
		image = fits.getdata(v)
		image = np.float64(image)
		image_data.append(image)	
	image_data= np.array(image_data)

	
	master_dark = np.mean(dark_data, axis=0)
	## this will probs need to change to do pixel by pixel
	
	master_bias = np.mean(bias_data, axis=0)

	flat_sub_bias = flat_data - master_bias
	master_flat = np.median(flat_sub_bias, axis=0)
	master_flat /= stats.mode(master_flat, axis=None)[0]

	reduced = (image_data - master_dark)/(master_flat)

	for idx,val in enumerate(reduced):
		image = fits.PrimaryHDU(val)
		image_name = os.path.basename(glob.glob(data_path)[idx])
		path = os.path.join(reduced_image_folder, image_name)
		image.writeto(path, overwrite=over_write)

	return plt.imshow(reduced[0], cmap='gray', norm=LogNorm(vmin=1000, 					vmax=1800), origin='lower left')
 
