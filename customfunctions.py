import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


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
    right = header['xllcorner'] + header['ncols'] * header['cellsize']
    bottom = header['yllcorner']
    top = header['yllcorner'] + header['nrows'] * header['cellsize']
    map_extent = (left, right, bottom, top)

    # Extract data values from the file and replace no data value with NaN
    data = np.loadtxt(file_name, skiprows=n_header, dtype='float')
    data[data == header['NODATA_value']] = np.nan

    return map_extent, data


def scoops3doutput_plots(run_id, dem_filepath, scoops3doutput_folderpath):
    """Plot maps to visualize DEM & output from Scoops 3D with custom colorbars

    Parameters:
        run_id (str): ID of the Scoops3D run
        dem_filepath (str): path to the DEM file in ASCII format
        scoops3doutput_folderpath (str): path to the folder with scoops 3D 
        output files
        
    Returns: 
        Saves 6 Plots in the scoops3doutput_folderpath defined in parameters.
        *dem_map (pdf): DEM
        *fos3d_map (pdf): Minimum 3D Factor of safety for each DEM cell 
        *fosvol_map (pdf): Volume of the critical surface for each DEM cell
        *searchgrid_map (pdf): Location of horizontal search space relative 
        to the DEM
        *numcols_map (pdf): Number of columns associated with the critical 
        surface at each DEM cell
        *critcheck_map (pdf): Check on the volumes associated with the critical 
        surfaces
        *boundcheck_map (pdf): Check on the search-lattice boundaries
        
    
    """

    ##########################################################################
    # Define file paths
    ##########################################################################

    # File path for DEM
    dem_filepath = dem_filepath

    # File path for Run Outputs
    fos3d_filepath = str(scoops3doutput_folderpath) + \
        str(run_id) + '_fos3d_out.asc'
    fosvol_filepath = str(scoops3doutput_folderpath) + \
        str(run_id) + '_fosvol_out.asc'
    critcheck_filepath = str(scoops3doutput_folderpath) + \
        str(run_id) + '_critcheck_out.asc'
    numcols_filepath = str(scoops3doutput_folderpath) + \
        str(run_id) + '_numcols_out.asc'
    searchgrid_filepath = str(scoops3doutput_folderpath) + \
        str(run_id) + '_searchgrid_out.asc'
    boundcheck_filepath = str(scoops3doutput_folderpath) + \
        str(run_id) + '_boundcheck_out.asc'

    ##########################################################################
    # Define colorscales
    ##########################################################################

    #  colorbar for DEM (dem)
    dem_cmap = mpl.cm.get_cmap("viridis").copy()
    dem_bounds = [0, 100, 250, 500, 750, 1000, 1250, 1500]
    dem_norm = mpl.colors.BoundaryNorm(dem_bounds, dem_cmap.N, extend='both')

    # colorbar for Factor of Safety (fos3d)
    fos3d_cmap = (mpl.colors.ListedColormap(['#f03b20', '#fd8d3c', '#feb24c',
                                             '#fed976', '#ffffb2', '#ffffcc',
                                             '#c7e9b4', '#7fcdbb', '#41b6c4',
                                             '#2c7fb8']).with_extremes(
                                                 over='#253494',
                                                 under='#bd0026'))
    fos3d_bounds = [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3]
    fos3d_norm = mpl.colors.BoundaryNorm(fos3d_bounds, fos3d_cmap.N)

    # color bar for failure volume (fosvol)
    fosvol_cmap = mpl.cm.YlGnBu.copy().reversed()
    fosvol_bounds = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    fosvol_norm = mpl.colors.BoundaryNorm(fosvol_bounds, fosvol_cmap.N,
                                          extend='both')

    # color bar for horizontal search space relative to DEM grid (searchgrid)
    searchgrid_cmap = (mpl.colors.ListedColormap(['red', 'orange',
                                                  'green', 'grey']))
    searchgrid_bounds = [-25, -1.5, 0, 2.5, 25]
    searchgrid_norm = mpl.colors.BoundaryNorm(searchgrid_bounds,
                                              searchgrid_cmap.N)

    # color bar for showing number of active columns (numcols)
    numcols_cmap = mpl.cm.cividis.copy()
    numcols_bounds = [100, 200, 300, 400, 500]
    numcols_norm = mpl.colors.BoundaryNorm(numcols_bounds, numcols_cmap.N,
                                           extend='both')

    # color bar for search lattice boudnary check file (boundcheck)
    boundcheck_cmap = mpl.cm.cividis.copy()
    boundcheck_bounds = [-0.1, 0.1]
    boundcheck_norm = mpl.colors.BoundaryNorm(boundcheck_bounds,
                                              boundcheck_cmap.N, extend='both')

    # colorbar for critical size check file (critcheck)
    critcheck_cmap = (mpl.colors.ListedColormap(['black', 'grey',
                                                 'blue', 'red']))
    critcheck_bounds = [-1.5, -0.5, 0.5, 1.5, 2.5]
    critcheck_norm = mpl.colors.BoundaryNorm(critcheck_bounds,
                                             critcheck_cmap.N)

    ##########################################################################
    # Create plots
    ##########################################################################
    # Read DEM map extent and Data
    dem_map_extent, dem = extract_mapextent_data(dem_filepath, 6)
    # Plot Elevation Map with Contours 100m interval
    fig, ax = plt.subplots(1)
    fig.set_size_inches(15, 10)
    ax.set_title('Elevation Map')
    dem_map = ax.imshow(dem,
                        origin='upper',
                        extent=dem_map_extent,
                        cmap=dem_cmap,
                        norm=dem_norm)
    ax.contour(dem, origin='upper', colors='black', extent=dem_map_extent,
               levels=list(range(0, int(np.nanmax(dem)), 100)))
    fig.colorbar(dem_map,
                 extend='both',
                 label='Elevation (m)')
    plt.ticklabel_format(style='plain')
    plt.savefig(str(scoops3doutput_folderpath) +
                '/{n}_dem_map.pdf'.format(n=run_id))

    # Factor of Safety
    # Extract FOS map extent and data
    fos3d_map_extent, fos3d = extract_mapextent_data(fos3d_filepath, 6)
    # Plot Factor of safety map
    fig, ax = plt.subplots(1)
    fig.set_size_inches(15, 10)
    ax.set_title('Factor of Safety Map', fontsize=16)
    fos3d_map = ax.imshow(fos3d,
                          origin='upper',
                          extent=fos3d_map_extent,
                          cmap=fos3d_cmap,
                          norm=fos3d_norm)
    ax.contour(dem, origin='upper', colors='black', extent=dem_map_extent,
               levels=list(range(0, int(np.nanmax(dem)), 100)))
    plt.ticklabel_format(style='plain')
    fig.colorbar(fos3d_map, extend='both',
                 label='Factor of Safety')
    plt.savefig(str(scoops3doutput_folderpath) +
                '/{n}_fos3d_map.pdf'.format(n=run_id))

    # Failure volume
    # Extract failure volume map extent and data
    fosvol_map_extent, fosvol = extract_mapextent_data(fosvol_filepath, 6)
    # Plot Failure volume map
    fig, ax = plt.subplots(1)
    fig.set_size_inches(15, 10)
    ax.set_title('Failure Volume Map', fontsize=16)
    fosvol_map = ax.imshow(fosvol/1000000000,
                           origin='upper',
                           extent=fosvol_map_extent,
                           cmap=fosvol_cmap,
                           norm=fosvol_norm)
    ax.contour(dem, origin='upper', colors='black', extent=dem_map_extent,
               levels=list(range(0, int(np.nanmax(dem)), 100)))
    plt.ticklabel_format(style='plain')
    fig.colorbar(fosvol_map,
                 extend='both',
                 label='Volume of critical failures $(km^3)$')
    plt.savefig(str(scoops3doutput_folderpath) +
                '/{n}_fosvol_map.pdf'.format(n=run_id))

    # Search Quality Evaluation - Horizontal Search Space
    # Extract search grid map extent and data
    searchgrid_map_extent, searchgrid = extract_mapextent_data(
        searchgrid_filepath, 6)
    # Plot search grid map
    fig, ax = plt.subplots(1)
    fig.set_size_inches(15, 10)
    ax.set_title('Horizontal Search Space Map', fontsize=16)
    searchgrid_map = ax.imshow(searchgrid,
                               origin='upper',
                               extent=searchgrid_map_extent,
                               cmap=searchgrid_cmap,
                               norm=searchgrid_norm)
    ax.set_xlim([472750, 478250])
    ax.set_ylim([6577750, 6583250])
    ax.contour(dem, origin='upper', colors='black', extent=dem_map_extent,
               levels=list(range(0, int(np.nanmax(dem)), 100)))
    plt.ticklabel_format(style='plain') 
    fig.colorbar(searchgrid_map,
                 ticks=[-11, -1, 2, 22]).set_ticklabels(
                     ['No search \nlattice \nnodes \nabove, \noutside DEM',
                     'No search \nlattice \nnodes \nabove, \inside DEM',
                     'Search lattice \nnodes above, \ninside DEM',
                     'Search lattice \nnodes above, \noutside DEM'])
    plt.savefig(str(scoops3doutput_folderpath) +
                '/{n}_searchgrid_map.pdf'.format(n=run_id))

    # Search Quality Evaluation - Number of active columns
    # Extract active columns map extent and data
    numcols_map_extent, numcols = extract_mapextent_data(numcols_filepath, 6)
    # Plot active columns map
    fig, ax = plt.subplots(1)
    fig.set_size_inches(15, 10)
    ax.set_title('Number of active columns intersected by critical surfaces',
                 fontsize=16)
    numcols_map = ax.imshow(numcols,
                            origin='upper',
                            extent=numcols_map_extent,
                            cmap=numcols_cmap,
                            norm=numcols_norm)
    ax.contour(dem, origin='upper', colors='black', extent=dem_map_extent,
               levels=list(range(0, int(np.nanmax(dem)), 100)))
    plt.ticklabel_format(style='plain')
    fig.colorbar(numcols_map,
                 extend='both')

    plt.savefig(str(scoops3doutput_folderpath) +
                '/{n}_numcols_map.pdf'.format(n=run_id))

    # Search Quality Evaluation - Critical size check
    # Extract critical size check map extent and data
    critcheck_map_extent, critcheck = extract_mapextent_data(
        critcheck_filepath, 6)
    # Plot critical size check map
    fig, ax = plt.subplots(1)
    fig.set_size_inches(15, 10)
    ax.set_title('Critical size check', fontsize=16)
    critcheck_map = ax.imshow(critcheck,
                              origin='upper',
                              extent=critcheck_map_extent,
                              cmap=critcheck_cmap,
                              norm=critcheck_norm)
    ax.contour(dem, origin='upper', colors='black', extent=dem_map_extent,
               levels=list(range(0, int(np.nanmax(dem)), 100)))
    plt.ticklabel_format(style='plain')
    fig.colorbar(critcheck_map,
                 ticks=[-1, 0, 1, 2]).set_ticklabels(
                     ['DEM cell was \nnot included \nin any trial surface',
                      'Size was \nnot \nrestricted',
                      'Volume was \nless than \nvmin + tol',
                      'Volume was \ngreater than \nvmax - tol'])

    plt.savefig(str(scoops3doutput_folderpath) +
                '/{n}_critcheck_map.pdf'.format(n=run_id))

    # Search Quality Evaluation - Boundary check
    # Extract boundary check map extent and data
    boundcheck_map_extent, boundcheck = extract_mapextent_data(
        boundcheck_filepath, 6)
    # Plot boundary check map
    fig, ax = plt.subplots(1)
    fig.set_size_inches(15, 10)
    ax.set_title('Search lattice boundary check', fontsize=16)
    boundcheck_map = ax.imshow(boundcheck,
                               origin='upper',
                               extent=boundcheck_map_extent,
                               cmap=boundcheck_cmap,
                               norm=boundcheck_norm)
    ax.contour(dem, origin='upper', colors='black', extent=dem_map_extent,
               levels=list(range(0, int(np.nanmax(dem)), 100)))
    plt.ticklabel_format(style='plain')
    fig.colorbar(boundcheck_map, extend='both', ticks=[0])

    plt.savefig(str(scoops3doutput_folderpath) +
                '/{n}_boundcheck_map.pdf'.format(n=run_id))

    return None
