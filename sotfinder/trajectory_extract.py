import numpy as np
from numpy.linalg import inv,norm
import pandas as pd
from glob import glob
import multiprocessing as mp
from functools import partial
from itertools import compress

from photutils.aperture import CircularAperture
from starextractor import parse_image

from .astroalign import find_transform_tree,matrix_transform
from .ransac import sequential_ransac
from .classes import Trace,Line,Stars

def arc_length_acceleration(x, y, t):
    """
    Calculates the tangential acceleration along the motion trajectory of a space object.

    This function is used to validate the feasibility of identified motion trajectories. It operates
    under the assumption that an object's tangential acceleration should not exhibit significant 
    variations over time. By examining the time series of tangential acceleration, the function
    helps to determine if the detected motion trajectory is plausible or not. Large variations in
    acceleration may indicate an invalid trajectory.

    Usage:
        >>> tangential_accel = arc_length_acceleration(x, y, t)
    Inputs:
        x -> [array-like] x pixel coordinates of the space object over time.
        y -> [array-like] y pixel coordinates of the object over time.
        t -> [array-like] Timestamps corresponding to x and y coordinates.
    Outputs:
        tangential_accel -> [array-like] Tangential acceleration values for each time moments.
    """

    # Calculate the time intervals between observations
    delta_t = np.diff(t)
    
    # Calculate the velocity components in x and y directions for each interval
    vx = np.diff(x) / delta_t
    vy = np.diff(y) / delta_t
    
    # Compute the speed (magnitude of velocity)
    speed = np.sqrt(vx**2 + vy**2)
    
    # Calculate the change in speed (delta speed)
    delta_speed = np.diff(speed)
    
    # Compute the tangential acceleration
    delta_time_for_acceleration = (delta_t[:-1] + delta_t[1:]) / 2
    tangential_acceleration = delta_speed / delta_time_for_acceleration
    
    # Prepend a NaN value to make the array size consistent
    tangential_acceleration = np.insert(tangential_acceleration, 0, np.nan)
    
    return tangential_acceleration

def iqr_outliers(x, scale=4.0):
    """
    Identifies outliers in a dataset using the Interquartile Range (IQR) method.

    Usage:
        outliers, is_outlier = iqr_outliers(x)
    Inputs:
        x -> [array-like] Data array in which to find outliers.
        scale -> [float] Scaling factor for the IQR to adjust sensitivity. A lower value makes the method more sensitive to outliers.
    Outputs:
        outliers -> [array-like] Outliers.
        is_outlier -> [array-like,bool] A boolean array marking the outliers in the original array.
    """

    # Compute the first and third quartiles (25th and 75th percentiles)
    Q1 = np.nanpercentile(x, 25)
    Q3 = np.nanpercentile(x, 75)

    # Calculate the IQR and determine outlier thresholds
    IQR = Q3 - Q1
    lower_bound = Q1 - scale * IQR
    upper_bound = Q3 + scale * IQR

    # Identify outliers
    is_outlier = (x < lower_bound) | (x > upper_bound)
    outliers = x[is_outlier]

    return outliers, is_outlier

def parse_astroimage_mp(params, arg):
    """
    Function designed for multi-core parallel processing of astronomical images. It parses an individual 
    image file and extracts relevant data, including star spots and their characteristics, for further analysis.

    The function performs the following steps for each image:
        1. Parse the image to extract star spots.
        2. Extract centroids, brightness, and SNR of these star spots.
        3. Use the brightest star spots to build asterisms, calculate their invariants, and construct a 2D tree.

    This function is typically called in parallel across multiple processor cores to significantly 
    improve the efficiency of processing a large set of astronomical images.

    Usage:
        >>> # In a multiprocessing context
        >>> with mp.Pool(processes=number_of_cores) as pool:
        >>>     results = pool.map(partial(parse_astroimage_mp, params), image_paths)

    Inputs:
        params -> [tuple]: Parameters for image parsing (max control points, FWHM).
        arg -> [str] -> Path of the astronomical image file.

    Outputs:
        result -> [tuple] Result of multi-core parallel processing.
    """

    mcp, fwhm, phot = params
    imagefile = arg
    
    try:
        # Parse the image and extract necessary data
        data = parse_image(imagefile)
        res, offset = data.res, data._offset

        # Extract star spots
        sources = data.find_source(phot=phot,max_control_points=mcp,fwhm=fwhm)

        # Build asterisms and calculate invariants from brightest star spots
        xy, asterisms, kdtree = sources.invariantfeatures(max_control_points=15)

        # Collect additional data for analysis
        pixels_camera, image, brightness, snr = sources.xy, data.image, sources.brightness, sources.snr
        camera = (xy, asterisms, kdtree)

        result = (res, offset, pixels_camera, image, brightness, snr, camera)
    except Exception as e:
        # Handling any errors that occur during image parsing
        print(f"Error processing image {imagefile}: {e}")
        result = None
    
    return result

def find_moving(image_dir, mcp=200, phot='aperture',scale=4):
    """
    Analyzes a series of astronomical images to identify moving space objects and their trajectories.
    It employs parallel processing and data manipulation with Pandas DataFrame for efficiency.

    The function performs the following steps:
        1. Reads and sorts astronomical images from a directory.
        2. Sets up parallel processing to analyze each image.
        3. Extracts star spots and their properties from each image.
        4. Identifies potential moving objects based on brightness variations and validates them.

    Usage:
        traces = find_moving('/path/to/image/directory/')
    Inputs:
        image_dir -> [str]: Directory containing a series of image files.
        mcp -> [int,optional,default=200] Maximum control points for detecting moving objects, used in 'find_source'. Higher values assist in detecting dimmer objects.
        phot -> [str,optional,default='aperture'] Method of photometry for sources, used in 'find_source'. Avaliable options include 'aperture', 'psf' and 'dao'. 
        scale -> [int,optional,default=4] Scale factor for the IQR method used in outlier detection.
    Outputs:
        Trace -> An object containing structured information about the trajectories, brightness, and signal-to-noise ratio of moving objects.
    """

    # Loading and sorting image files from the specified directory
    image_list = sorted(glob(image_dir + '*'))
    n = len(image_list)  # Total number of images in the directory
    mid = n // 2  # Index of the intermediate frame used as a reference

    # Processing the intermediate image to determine FWHM (Full Width at Half Maximum)
    data_mid = parse_image(image_list[mid])
    fwhm = data_mid.fwhm  # Full Width at Half Maximum in pixels

    # Set up multiprocessing for parallel processing of images
    number_of_cores = mp.cpu_count() - 1  # Number of available processor cores
    params = (mcp, fwhm, phot)  # Parameters for image parsing
    func = partial(parse_astroimage_mp, params)  # Function to be executed in parallel

    # Process all images in parallel and gather results
    with mp.Pool(number_of_cores) as pool:
        reslist = pool.map(func, image_list)

    # Extract results of parallel computation for the intermediate frame the first frame
    res,offset,pixels_camera_mid,image_mid,brightness_mid,snr_mid,camera_mid = reslist[mid]
    res,offset,pixels_camera0,image0,brightness0,sne0,camera0 = reslist[0]

    # Initialize DataFrame and data lists
    df = pd.DataFrame(columns=['frame','transf'])
    xy_candis, transf_matrices, xys, fluxes, snrs, images = [], [], [], [], [], []
    # Dictionary to accumulate trace information for final output
    trace_info = {'res':res,'offset':offset,'fwhm':fwhm,'_ref_ind':mid}

    # Process each image in the series
    for i in range(n):
        imagefile = image_list[i]
        # Extract the result of parallel computation for the current frame
        res,offset,pixels_camera,image,brightness,snr,camera = reslist[i]

        # Append extracted data to respective lists
        xys.append(pixels_camera) # Pixel coordinates of detected sources
        fluxes.append(brightness) # Brightness values of sources
        snrs.append(snr) # Signal-to-noise ratios of sources
        images.append(image) # Image of background subtracted

        # Calculate the affine transformation and identify candidate coordinates for moving objects
        if i == mid:
            # Calculate affine transformation from the intermediate frame to the first frame
            transf,(pixels_camera_mid_match,pixels_camera0_match),_s,_d = find_transform_tree(camera_mid,camera0)
            # Transform coordinates from the intermediate frame to the first frame
            pixels_camera_affine = matrix_transform(pixels_camera_mid,transf.params)
            # Create circular apertures centered on the affine coordinates
            apertures_affine = CircularAperture(pixels_camera_affine+offset, r=fwhm)
            # Perform photometry within these apertures on the first frame
            brightness_affine = apertures_affine.do_photometry(image0)[0] 

            # Identifies potential candidates(outliers) for moving objects using the Interquartile Range (IQR) method.
            # In this context, outliers are defined as sources whose brightness in the local frame significantly
            # differs from their brightness after affine transformation to the reference frame.
            _,flags = iqr_outliers(np.abs(brightness_mid - brightness_affine),scale)
            # Initialize transformation matrix as identity for the intermediate frame
            transf_matrix = np.eye(3)
            xy_candi = pixels_camera_mid[flags] # Pixel coordinates of candidates for the intermediate frame
        else:
            # For other frames, calculate transformation to the reference frame
            transf,(pixels_camera_match, pixels_camera_mid_match),_s,_d = find_transform_tree(camera,camera_mid)
            # # Transform coordinates to the reference frame
            pixels_camera_affine = matrix_transform(pixels_camera,transf.params)
            # Circular apertures and photometry similar to the above
            apertures_affine = CircularAperture(pixels_camera_affine+offset,r=fwhm)
            brightness_affine = apertures_affine.do_photometry(image_mid)[0] 
            _,flags = iqr_outliers(np.abs(brightness - brightness_affine),scale)
            transf_matrix = transf.params # Affine matrix from the current frame to the reference frame
            xy_candi = pixels_camera_affine[flags] # Pixel coordinates of candidates for the current frame
            
        # Update DataFrame with image file name and transformation matrix
        df.loc[i] = [imagefile,transf_matrix]

        # Handle cases with no candidate moving objects
        xy_candis.append(np.full((1, 2), np.nan) if not flags.any() else xy_candi) # Assign NaN if no candidates found
        transf_matrices.append(transf_matrix)

    # Convert list of transformation matrices and images to a NumPy array for efficient processing
    transf_matrices = np.array(transf_matrices)
    images = np.array(images)

    # Concatenate candidate coordinates, remove NaN values, and separate into x and y
    xy_array = np.concatenate(xy_candis) 
    xy_dropna = xy_array[~np.isnan(xy_array).all(axis=1)]
    x,y = xy_dropna.T

    # Identify trajectories of moving objects using the sequential RANSAC algorithm
    lines = sequential_ransac(x,y) 
    m = len(lines.keys()) # Number of detected trajectories

    if m != 0:
        iis = [] # List to store indices of images with valid trajectories
        for key in lines.keys(): 
            xy_affine = lines[key] # Coordinates of the moving object in the reference frame
            flags,xy_inline = [],[] # List to store flags for valid frames and list to store in-trajectory coordinates

            for xy_candi in xy_candis:
                # Check if candidate points are on the trajectory
                in_line = np.isin(xy_candi,xy_affine).all(axis=1)
                if in_line.sum() == 1:
                    flag = True # Valid frame if exactly one candidate is on the trajectory
                    xy_inline.append(xy_candi[in_line]) # Collect valid candidate point
                else:    
                    flag = False # Invalid frame if none or multiple candidates on the trajectory
                flags.append(flag)

            flags = np.array(flags)  
            xy_inline = np.concatenate(xy_inline)

            # Validate trajectories by checking for consistent tangential acceleration
            x,y = xy_inline.T
            t = np.where(flags)[0] # Time points corresponding to valid frames
            tan_acc = arc_length_acceleration(x,y,t) # Calculate tangential acceleration
            _,tan_acc_flags = iqr_outliers(tan_acc,scale)

            if tan_acc_flags.any(): 
                continue # Skip trajectory if it has irregular acceleration

            # Update DataFrame with trajectory data
            df.loc[flags,'xy_affine_line'+key] = pd.Series(xy_inline.tolist()).values

            # Process each valid frame for the trajectory
            df_xy,df_flux,df_snr = [],[],[] # Lists for coordinates, flux, and SNR
            transf_valid = list(compress(transf_matrices, flags))
            xys_valid= list(compress(xys, flags))
            fluxes_valid = list(compress(fluxes, flags))
            snrs_valid = list(compress(snrs, flags))

            inds = [] # Indices of the space object among sources in the valid frames
            for j in range(len(xy_inline)):
                # Transform affined coordinates back to their original frame for each valid point on the trajectory
                df_xy_j = matrix_transform(xy_inline[j],inv(transf_valid[j]))[0]
                # Find the closest point in the original frame to the transformed point
                ind = np.argmin(norm(df_xy_j - xys_valid[j],axis=1))
                inds.append(ind) # Store the index of the closest point among sources
                # Append coordinates, flux, and SNR of the space object to the lists
                df_xy.append(xys_valid[j][ind])
                df_flux.append(fluxes_valid[j][ind])
                df_snr.append(snrs_valid[j][ind])

            # Update DataFrame with local frame data
            df.loc[flags, f'xy_local_line{key}'] = pd.Series(df_xy).values
            df.loc[flags, f'flux_line{key}'] = pd.Series(df_flux).values
            df.loc[flags, f'snr_line{key}'] = pd.Series(df_snr).values

            # Store trajectory information
            line_info = {'xy_affine':xy_inline,'xy_local':df_xy,'flux':df_flux,'snr':df_snr,'_flags':flags}
            trace_info[f'line{key}'] = Line(line_info)

            # Create an array 'ii' with the same shape as 'flags'. Initialize it with -9, which
            # represents invalid or non-existent indices. This array will be used to track the
            # indices of space objects among sources extracted in each frame that are part of the detected trajectories.
            # Valid indices will be updated later where the trajectory is present.
            ii = np.full_like(flags,-9,dtype=int)
            ii[flags] = inds
            iis.append(ii)

        # Creates an array of shape (n, m) where 'n' is the number of frames and 'm' is the number of trajectories.
        # Each element in this array represents the indices of multiple space objects in a specific frame for specific trajectories.
        iis = np.stack(iis).T

        # For each frame (each row in the 'iis' array), filter out indices marked as -9.
        # -9 indicates an invalid or non-existent index, implying that the trajectory is not present in that frame.
        iis = [kk[kk != -9] for kk in iis]

        # Remove trajectory points from sources
        for i in range(n):
            xys[i] = np.delete(xys[i],iis[i],axis=0)
            fluxes[i] = np.delete(fluxes[i],iis[i])
            snrs[i] = np.delete(snrs[i],iis[i])

    # Combine star data into DataFrame
    df = df.join(pd.DataFrame({'xy_local_star':xys,'flux_star':fluxes,'snr_star':snrs}))

    # Prepare final trace information
    stars_info = {'xy_local':xys,'flux':fluxes,'snr':snrs,'_n_mid':len(snrs[mid])}
    trace_info['stars'] = Stars(stars_info)
    trace_info.update({'df':df,'frames':image_list,'transf':transf_matrices,'_num_frames':n,'_num_traces':m,'_images':images,'_xy_candis':xy_dropna})

    return Trace(trace_info)