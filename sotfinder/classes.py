import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path
import h5py  # Importing the h5py library for HDF5 file handling

def from_hdf5(file_hdf5):
    """
    Read data from an HDF5 file and construct a Trace object containing various astronomical data.
    This function opens the specified HDF5 file and reads various data related to astronomical observations,
    such as resolution, Full Width Half Maximum (FWHM), offsets, frames, reference frame indice, stars, and traces of space objects.

    Usage:
        >>> from sotfinder import from_hdf5
        >>> trace = from_hdf5('path_to_hdf5_file.hdf5')
    Inputs:
        file_hdf5 -> [str] Path to the HDF5 file.
    Outputs:
        Trace -> A Trace object containing structured data extracted from the HDF5 file.

    The function reads the following key data from the HDF5 file:
        - 'res': Resolution of the observation image.
        - 'fwhm': Full Width Half Maximum (FWHM) values of the image.
        - 'offset': The coordinates of the image center relative to the original coordinate origin (the center of the lower-left pixel).
        - 'ref_ind': Index of the reference frame.
        - 'transf': Affine transformation matrix to the reference frame.
        - 'frames': image frames.
        - 'stars/xy_local': Local pixel coordinates of stars in each frame.
        - 'stars/flux': Flux of stars.
        - 'stars/snr': Signal-to-noise ratio of stars.
        - 'lines': Data representing the motion trails of space objects in the observation.

    For each line (space object's motion trail):
        - 'xy_affine': Affined pixel coordinates in the reference frame.
        - 'xy_local': Local pixel coordinates in each frame.
        - 'flux': Flux values of the space object on trail.
        - 'snr': Signal-to-noise ratio.
        - 'flags': Indicators for the presence of the space object in the image frame.

    The function closes the HDF5 file after reading all necessary data.
    """

    # Open the HDF5 file for reading
    with h5py.File(file_hdf5, 'r') as fin:

        # Extracting basic information
        res = tuple(fin['res'][:])
        fwhm = fin['fwhm'][()]
        offset = fin['offset'][:]
        ref_ind = fin['ref_ind'][()]

        # Extracting frame names, transformation matrix, and images
        frames = list(fin['frames'][:])
        transf = fin['transf'][:]
        images = fin['images'][:]

        # Preparing trace information
        trace_info = {
            'res': res,
            'fwhm': fwhm,
            'offset': offset,
            '_ref_ind': ref_ind,
            'frames': frames,
            'transf': transf,
            '_images': images
        }

        # Extracting star information
        xy_local_star = [ele.reshape(len(ele)//2, 2) for ele in fin['stars/xy_local'][:]]
        # 'xy_local_star' contains the reshaped local coordinates of stars in each frame

        flux_star = fin['stars/flux'][:].tolist()
        # 'flux_star' holds the flux values of stars, converted to a list for easier manipulation

        snr_star = fin['stars/snr'][:].tolist()
        # 'snr_star' contains the signal-to-noise ratio of stars, also converted to a list

        # Organizing star information into a dictionary for the Trace object
        stars_info = {
            'xy_local': xy_local_star,
            'flux': flux_star,
            'snr': snr_star,
            '_n_mid': len(snr_star[ref_ind]) # The number of stars in the reference frame
        }
        # Adding the structured star information to 'trace_info' for further processing
        trace_info['stars'] = Stars(stars_info)

        # Processing lines representing space objects' motion trails
        lines = fin['lines']
        for key in lines.keys():
            # Extracting affine and local coordinates, flux, SNR, and flags for each trace
            xy_affine = lines[f'{key}/xy_affine'][:]
            xy_local = lines[f'{key}/xy_local'][:]
            flux = lines[f'{key}/flux'][:]
            snr = lines[f'{key}/snr'][:]
            _flags = lines[f'{key}/_flags'][:]

            # Structuring line information into a dictionary
            line_info = {
                'xy_affine': xy_affine, # Affined pixel coordinates in the reference frame
                'xy_local': xy_local,   # Local pixel coordinates in each frame
                'flux': flux,           # Flux values for the trace
                'snr': snr,             # Signal-to-noise ratio for trace
                '_flags': _flags         # Flags indicating presence of trace in the image frame
            }
            # Adding the structured line information to the trace info
            trace_info[key] = Line(line_info)

        # Update trace information with the number of frames and traces
        trace_info.update({
            '_num_frames': len(frames), # Total number of frames processed
            '_num_traces': len(lines.keys()) # Total number of identified traces
        })

    # Return the structured Trace object with all gathered information
    return Trace(trace_info)

class Trace(object):
    """
    Class representing a collection of astronomical observation data. It includes information about
    resolution, stars, trajectories of space objects, and other relevant astronomical data.

    Attributes:
        res (tuple): Resolution of the image frame. 
        Typically stored as a tuple (width, height), such as (1024,1024), representing the number of pixels in each dimension of the image.
        fwhm (float): Full Width at Half Maximum of the Gaussian kernel for the image.
        offset (tuple): Pixel coordinates of the image center relative to the lower-left pixel center, such as (511.5,511.5).
        frames (list): The collection of image frames in the observation. Each frame captures a view of the sky at a specific time.
        transf (2D numpy array): Transformation matrices. 
        Used for image alignment or to map coordinates from one frame to another. Essential for comparing and combining data across multiple frames.
        stars (Stars): Information about the stars in the observation.
        In form of a complex object about identified stars, such as their positions, brightness, and signal-to-noise ratios.
        line0,line1,.. (Line): The motion trails of space objects.  It contains affined/local coordinates, brightness, signal-to-noise ratios, and other relevant data.
        _ref_ind (int): Index of the reference frame in the series of observations.
        _num_frames (int): Total number of frames in the observation.
        _num_traces (int): Number of identified traces or trajectories.
        _images (3D numpy array): The background-subtracted grayscale image.
        _xy_candis (2D numpy array): Candidate coordinates for moving objects.
        It is crucial for further analysis, particularly when fine-tuning the RANSAC algorithm for trajectory detection.

    Methods:
        savefile: Saves the data into an HDF5 file, preserving its structure.
        show: Displays images based on the astronomical data.
    """
    def __init__(self, info):
        """
        Initializes a Trace instance with information from an astronomical observation.
        """
        for key in info.keys():
            setattr(self, key, info[key]) # Efficiently setting attributes from the info dictionary

    def __repr__(self):
        """
        Returns a detailed string representation of the Trace object.
        """
        return '<Trace object: RES = {} FWHM ≈ {:.2f} REF_INDEX = {:d} NUM_FRAMES = {:d} NUM_TRACES = {:d}>'.format(
            self.res, self.fwhm, self._ref_ind, self._num_frames, self._num_traces)

    def savehdf5(self, file_out):
        """
        Saves the Trace data to an HDF5 file. This method ensures that the complex
        structure of the Trace data is preserved and stored in an organized manner.

        Inputs:
            file_out -> [str] Path of the output HDF5 file. Must end with '.hdf5'.

        Raises:
            ValueError -> If the file extension is not '.hdf5'.
        """
        # Ensuring the file format is HDF5
        if not file_out.endswith('.hdf5'):
            raise ValueError("Output file format must be '.hdf5'.")

        # Create directories if they don't exist
        Path(file_out).parent.mkdir(parents=True, exist_ok=True)

        # Define a variable-length data type for storing arrays of float64
        vlen_dtype = h5py.special_dtype(vlen=np.dtype('float64'))

        # Opening the HDF5 file for writing
        with h5py.File(file_out, 'w') as fout:
            # Writing basic information like resolution, FWHM, and offset
            fout.create_dataset("res", data=self.res)
            fout.create_dataset("fwhm", data=self.fwhm)
            fout.create_dataset("offset", data=self.offset)

            # Saving star information using variable-length data type
            stars_grp = fout.create_group("stars")
            stars = self.stars

            num_frames = self._num_frames
            xy_local_star = stars_grp.create_dataset('xy_local', (num_frames,), dtype=vlen_dtype)
            flux_star = stars_grp.create_dataset('flux', (num_frames,), dtype=vlen_dtype)
            snr_star = stars_grp.create_dataset('snr', (num_frames,), dtype=vlen_dtype)

            for i in range(num_frames):
                xy_local_star[i] = stars.xy_local[i].flatten()
                flux_star[i] = stars.flux[i]
                snr_star[i] = stars.snr[i]

            # Saving trace information for each identified line (motion trail)
            lines_grp = fout.create_group("lines")
            for key, value in self.__dict__.items():
                if 'line' in key:
                    subgrp = lines_grp.create_group(key)
                    for sub_key, sub_val in value.__dict__.items():
                        subgrp.create_dataset(sub_key, data=sub_val)

            # Create datasets for other Trace attributes
            fout.create_dataset("frames", data=self.frames)
            fout.create_dataset("transf", data=self.transf)
            fout.create_dataset("ref_ind", data=self._ref_ind)
            fout.create_dataset("images", data=self._images)

    def show(self, mode='trace', fig_out=None):
        """
        Visualizes the astronomical data based on the specified mode.

        Inputs:
            mode -> [str,optional,default='trace'] Determines the visualization mode. Available Options includ 'trace' and 'image', where
                        'trace' mode plots the trajectories of moving objects.
                        'image' mode displays the positions of moving objects and stars in a specific frame.
            fig_out -> [str, optional, default=None): Path to save the output figure. If None, the figure is displayed inline.
        Outputs:
            If fig_out is specified, saves the figure to the given path. Otherwise, displays the figure on screen. 
        Raises:
            ValueError: If an unknown mode is specified.
        """

        if mode == 'trace':
            self._plot_trajectories(fig_out)
        elif mode == 'image':
            self._plot_image(fig_out)
        else:
            raise ValueError(f"Unknown mode: {mode}")

    def _plot_image(self, fig_out=None):
        """
        Plots both moving space objects and stars in a specific frame.

        Inputs:
            fig_out -> [str, optional, default=None): Path to save the output figure. If None, the figure is displayed inline.
        Outputs:
            If fig_out is specified, saves the figure to the given path. Otherwise, displays the figure on screen.    
        """
        '''
        images = self._images
        xy = self.xy_local + offset
        xy_star = self.xy_local + offset
        if fig_out is None:
            plot_kwargs = {'mark':(xy,'s','red')}
        else:
            plot_kwargs = {'mark':(xy,'s','red'),'figname':fig_out}
        show_image(image,origin='lower',**plot_kwargs) 
        '''
        pass        

    def _plot_trajectories(self, fig_out=None):
        """
        This method is called by the 'show' method to create a scatter plot of the space objects' trajectories.
        Each point in the trajectories is plotted as a small dot, and the trajectories are overlaid on a grid for reference.

        The plotting area is adjusted to be 1.5 times larger than the resolution of the images
        to ensure that trajectories extending beyond the reference frame image are fully visible.

        Inputs:
            fig_out -> [str, optional, default=None] Path to save the output figure. If None, the figure is displayed 
            on screen.
        Outputs:
            If fig_out is specified, saves the figure to the given path. Otherwise, displays the figure on screen.
        """
        fig, ax = plt.subplots(dpi=300)

        with plt.style.context(('seaborn-darkgrid')):
            # Iterate over each trajectory and plot it
            for key, line in self.__dict__.items():
                if key.startswith('line'):
                    x, y = line.xy_affine.T
                    plt.scatter(x, y, marker='o',s=5,label=f'Trajectory {key}')            

        ax.set_aspect('equal', 'box')
        plt.xlabel("x")
        plt.ylabel("y")

        # Extend the plot's display range to accommodate trajectories that extend beyond the reference frame's bounds
        height,width = res = self.res
        h_half,w_half = res[0]/2, res[1]/2
        xb, yb = w_half * 1.5, h_half * 1.5
        plt.xlim([-xb, xb])
        plt.ylim([-yb, yb])

        # Add boundary of the frame
        ax.add_patch(Rectangle([-w_half,-h_half], width, height,fill=False,lw=1,color='k',ls='dashdot'))

        if fig_out is None:
            plt.show()
        else:
            Path(fig_out).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(fig_out, bbox_inches='tight')

class Line(object):
    """
    Represents the trajectory of a moving space object as observed in the astronomical data.

    Attributes:
        xy_affine -> [array-like] Coordinates of the trajectory in the reference frame.
        xy_local -> [array-like] Coordinates of the trajectory in local frame.
        flux -> [array-like] Flux values associated with points on the trajectory.
        snr -> [array-like] Signal-to-noise ratio for each point in the trajectory.
        _flags -> [array-like,bool] Flags indicating the presence of the moving object in each image frame.
    """
    def __init__(self,info): 
        """
        Initializes a Line instance with trajectory information.

        Inputs:
            info -> [dict] A dictionary containing trajectory data including coordinates, flux, SNR, and presence flags.
        """ 

        for key, value in info.items():
            setattr(self, key, np.array(value))  # Convert lists to NumPy arrays for efficient handling    

    def __repr__(self):
        """
        String representation of the Line object, showing key information.
        """
        return f'<Line object: NUM_POINTS = {len(self.xy_affine)}>'

class Stars(object):
    """
    Represents the stars observed in the astronomical images.

    Attributes:
        xy_local [list of 2d array] Coordinates of stars in each local frame. 
        The list length corresponds to the number of frames, and each 2D array has 
        shape (k, 2), where 'k' is the number of stars in that frame.
        flux [list of 1d array] Flux values of stars across all frames.
        snr [list of 1d array] Signal-to-noise ratio for stars across all frames.
        _n_mid (int): Number of stars in the reference frame.
    """
    def __init__(self,info): 
        """
        Initializes a Stars instance with information about the stars.

        Inputs:
            info -> [dict] A dictionary containing data about the stars, including their coordinates, flux, and SNR.
        """
        for key, value in info.items():
            setattr(self, key, value)  # Assigning attributes from the info dictionary    

    def __repr__(self):
        """
        String representation of the Stars object, summarizing key data.
        """
        return f'<Stars object: NUM_STARS = {self._n_mid}>'
