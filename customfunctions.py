import numpy as np

def extract_mapextent_data(file_name, n_header):
    """Extract map extent and data.

    Parameters:
        file_name (str): path to the map file in ASCII format
        n_header (int) :- number of header lines 
    
    Returns: 
        map_extent (tuple): map extent with left, right, bottom, and 
                            top coordinates in UTM
        data (numpdy array): map data  
    
    """
    # Read file header
    header = {}
    row = 1
    with open(file_name, 'rt') as file_header:
        for line in file_header:
            if row <= n_header:
                line = line.split(" ", 1)
                header[line[0]] = float(line[1])
            else:
                break
            row = row + 1
    
    # Compute map extent
    left = header['xllcorner']
    right = header['xllcorner'] +  header['ncols'] * header['cellsize']
    bottom = header['yllcorner']
    top = header['yllcorner'] + header['nrows'] * header['cellsize']
    map_extent = (left, right, bottom, top)
    
    # Extract data values from the file and replace no data value with NaN
    data = np.loadtxt(file_name, skiprows = n_header, dtype ='float')
    data[data == header['NODATA_value']] = np.nan
    
    
    return map_extent, data
