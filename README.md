# Function-Junction
A repository to dump files I make while working.

All functions are written for Python 3 unless otherwise specified.

diffmag: A set of functions (getindex, matchLists, find_diff_mag) that can be used to find variability in magnitude of stars in a data set. The first two functions match between two sets of data based on ra and dec while the last finds the variability.

get_hjd_better: A script that will obtain Heliocentric Julian Dates for your observations from the MJD keyword in fits file headers. The input is the red and blue file paths and the output file name. It can alternatively be only one file path and the output file name.

mags_from_spec_better: A script that will find the polarization and position angle in all filters for multiple days worth of polarization spectra. All it takes as input is the path string to your blue range spectra, the path string to your red range spectra, and the output file name. It will save a csv table filled with your values and their errors.

reduction_python: This function (reduce_data) will take the paths to your data and calibrations and do basic dark, bias, and flat removal using Python.
