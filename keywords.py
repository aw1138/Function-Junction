def get_keyword(folders, images, keyword, write_name=None):
    """
    
    A function to get any chosen single keyword from the header of a fits image 
    or set of fits images. It will extract the value of the chosen keyword, and
    then returns the table written to the file as a pandas DataFrame object. It
    writes the table to a file with the chosen name if write_name has been 
    filled in.
    
    Attributes
    
        folders: list of strings
            The folder(s) the images are in. Must be iterable set of strings
            with each string being a path to the folder the images are in.
            If there is more than one folder to pull images from, there will be 
            more than one item in the set.
        images: string
            The image or image pattern to search for and pull header information
            from. Can include wild cards (*), etc.
        keyword: string
            The keyword in the header of the fits image(s) to be pulled. 
        write_name: optional, string
            The name of the file the output table is to be saved to. The default 
            is None.
    
    """
    import glob
    import os.path
    from astropy.io import fits
    import pandas as pd
    
    dates_jd = []
        
    for i,v in enumerate(folders):
        file_path = os.path.join(v, images)

        for i,v in enumerate(glob.glob(file_path)):
            image = fits.open(v)
            date = image[0].header[keyword]
            image.close(v)
            dates_jd.append(date)

    dates_final = pd.DataFrame(dates_jd, columns=[keyword])
    
    if write_name is not None:
        dates_final.to_csv(write_name, index=False)
    
    return dates_final