# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 13:23:53 2019

@author: casimp
"""

import numpy as np
import csv
import matplotlib.pyplot as plt
from cpex.transformation import strain_transformation

class Extract():
    
    def __init__(self):
        pass
        
    def extract_grains(self, data='elastic', idx=1, grain_idx=None):
        """
        Extracts data (stress, strain etc.) for either all grains, or a 
        specified grain at a given (orthogonal) orientation or component
        (where data='time', 'frame' etc. indexing does not work.
        
        Parameters
        ----------
        data: str
            The data label, either 'stress', 'strain', 'elastic' (strain),
            'back stress', 'rot', 'time', 'frame'
        idx: int
            The orientation (referenced via an idx) of the defined data
            e.g. data='stress', idx=1 => sigma_yy
        grain_idx: int
            The index of the grain (note GRAIN-1 => idx=0)
            
        Returns
        -------
        order: array
            Dataset
        """
        if idx == None and grain_idx != None:
            idx = np.s_[:, grain_idx]
        elif idx == None and grain_idx == None:
            idx = np.s_[:, :]
        elif idx != None and grain_idx == None:
            idx = np.s_[idx, :]
        else:
            idx = np.s_[idx, grain_idx]
        
        d = {'strain':self.e,
             'stress':self.s,
             'elastic':self.elastic,
             'back stress':self.b_stress,
             'rot':self.rot - self.rot[:,:, 0][:, :, None],
             'time':self.t,
             'frame':np.arange(self.num_frames)}
        
        if data not in ['time', 'frame', 'rot']:
            ex = d[data][idx]
        else:
            ex = d[data]
            
        return ex

    
    def extract_neighbours_idx(self, grain_idx, frame=0):
        """
        Extracts the indinces of all grains ordered with respect to position 
        away from a given grain (index).
        
        Grains move a small amount during deformation, the frame can be defined
        to explicity interrogtae neightbours at a given load level/time.
        
        Parameters
        ----------
        grain_idx: int
            The index of the grain to search around
        frame: int, None
            The frame to reference (default = 0). If None extracts ordered 
            inidices for all frames.
            
        Returns
        -------
        order: list
            Grain indices ordered by euclidean distance from selected grain
        """
        
        if frame == None:
            frame = np.s_[:]
        dims = self.dims[:, :, frame]
        rel_dims = dims - dims[:, grain_idx, None] # Keeps correct dimensionality
        euc_dist = np.sum(rel_dims ** 2, axis=0)**0.5
        
        order = np.argsort(euc_dist, axis=0)
        
        return order[1:]
    
    
    def extract_neighbours(self, grain_idx, data='strain', idx=1, frame=-1,
                           cmean=False, dimframe='simple'):
        """
        Extracts data (stress, strain etc.) for all grains, with data being 
        ordered with respect to position away from a given grain (index).

        Calls extract_grains and extrains_neighbours_idx methods.
        
        Parameters
        ----------
        grain_idx: int
            The index of the grain to search around
        data: str
            The data label, either 'stress', 'strain', 'elastic' (strain),
            'back stress'
        idx: int
            The orientation (referenced via an idx) of the defined data
            e.g. data='stress', idx=1 => sigma_yy
        frame: int, None
            The frame to reference (default = 0). If None extracts ordered 
            data for all frames.
        cmean: bool
            Compute a rolling, cumulative mean
        dimframe: str, int, None
            If frame is not None then the neighbour ordering is done on same
            frame. If frame is None then the dimensions are taken from the 
            final frame or a specified frame (int) unless dimframes==None, in
            which case neighbour ordering is done for each frame . Warning 
            this is slow!
            
        Returns
        -------
        order: array
            Ordered dataset
        """
        if frame==None:
            frame=np.s_[:]
            
        
        ex = self.extract_grains(data=data, idx=idx)
        
        if frame == np.s_[:] and dimframe !='simple':
            # Not ideal!!!
            order = self.extract_neighbours_idx(grain_idx, None)
            ex_ordered = np.column_stack([i[j] for i, j in zip(np.split(ex, ex.shape[1], axis=1), 
                                          np.split(order, order.shape[1], axis=1))]).squeeze()
        
        elif frame == np.s_[:] and isinstance(dimframe, int):
            order = self.extract_neighbours_idx(grain_idx, dimframe)
            ex_ordered = ex[order]
        
        else:
            dimframe = frame if frame != np.s_[:] else -1
            order = self.extract_neighbours_idx(grain_idx, dimframe)
            # print(order)
            ex_ordered = ex[order]
            
        if cmean:
            ex_csum = np.cumsum(ex_ordered, axis=0) 
            ex_cmean = ex_csum / np.arange(1, ex_csum.shape[0] + 1)[:, None]
        
            return ex_cmean[..., frame]
        
        return ex_ordered[..., frame]
        
        
    def plot_neighbours(self, grain_idx, data='plastic', idx=1, frame=-1,
                        cmean=True, ):
        """
        Plots data (stress, strain etc.) for all n grains, with data being 
        ordered with respect to position away from a given grain (index).
    
        
        Parameters
        ----------
        grain_idx: int
            The index of the grain to search around
        data: str
            The data to plot either 'stress', 'strain', 'elastic' (strain), 
            'back stress'
        idx: int
            The orientation (referenced via an idx) of the defined data
            e.g. data='stress', idx=1 => sigma_yy
        frame: int, None
            The frame to reference (default = 0). If None extracts ordered 
            data for all frames.
        cmean: bool
            Compute a rolling, cumulative mean
        """
        assert frame != None, "Can't study response across all frames."
        
        ex_ordered = self.extract_neighbours(grain_idx, data=data, 
                                             idx=idx, frame=frame,
                                             cmean=cmean)
        
        # Tinkering with axis labels
        x = 'nth nearest neighbour'
        y = 'cumulative mean {} (window=n)'.format(data) if cmean else data
        
        # Plotting
        plt.plot(np.arange(1, np.size(ex_ordered) +1), ex_ordered, label=grain_idx)
        plt.legend()
        plt.ylabel(y)
        plt.xlabel(x)
    
    def extract_lattice(self, data='lattice', family='311', 
                        grain_idx=None, plane_idx=None):
        """
        Routine to extract information about some or all (default) grains for a 
        specified lattice plane.

        Parameters:
        -----------
        data: str
            Either 'lattice' or 'phi'
        family: str
            The lattice plane family to assess
        grain_idx: int, [int,...], None
            If None then all grains of this family to be extracted else
            the individual grain (or list of grains)
        plane_idx: int, [int,...], None
            If None then all planes of this family/grain combination to be 
            extracted else the individual planes (or list of planes)
            
        Returns:
        --------
        data: array
            Lattice strains (or phi) for given family (and potentially 
            grain/plane specification)
        """
        
        
        if plane_idx == None and grain_idx != None:
            idx = np.s_[:, grain_idx]
        elif plane_idx == None and grain_idx == None:
            idx = np.s_[:, :]
        elif plane_idx != None and grain_idx == None:
            idx = np.s_[plane_idx, :]
        else:
            idx = np.s_[plane_idx, grain_idx]
        
        lattice = self.lattice_strain[family][idx]
        phi = self.lattice_phi[family]
        
        d = {'phi':phi,'lattice':lattice}
        
        return d[data]
    
    def extract_phi_idx(self, family='311', phi=0, window=10, frame=0):
        """
        Allows for selection of the index of lattice planes wityh a defined 
        orientation with resepect to the y axis (nominally the loading axis).
        A 2D array of indices with be returned if a frame is specified, the
        elemtns in the array will be structured:
            
            [[grain_idx, plane_idx],
            [grain_idx, plane_idx],
            ...]
        
        If None is passed as the frame variable then the rotation of
        the grain during loading/dwell etc. is being considered - a 2D array 
        is returned with each element being structured as follows:
            
            [[grain_idx, frame_idx, plane_idx],
            [grain_idx, frame_idx, plane_idx],
            ...]
            
        ** In addition to the list of indices an equivalent boolean array is 
        returned in each case. **
        
        Parameters
        ----------
        family: str
            The index of the grain to search around
        phi: float
            The data to extractm either 'stress', 'strain', 'elastic' (strain),
            'back stress'
        window: float
            The orientation (referenced via an idx) of the defined data
            e.g. data='stress', idx=1 => sigma_yy
        frame: int, None
            The frame to reference (default = 0). If None extracts ordered 
            data for all frames.
            
        Returns
        -------
        va: array (bool)
            Boolean array of the same dimension as the lattice strain array -
            elements are True if they are within the window, else False)
        select: array (int)
            A list of the grains/plane indices for all grains that lie within 
            specified orientation/window combination.
            
        """
        if frame == None:
            frame = np.s_[:]
            
        phi_ = 180 * self.lattice_phi[family][:, frame] / np.pi
        
        phi_ -= 90
        phi -= 90
        w = window / 2
        p0, p1 = phi - w, phi + w
        
        s0 = np.logical_and(phi_ > np.min(p0), phi_ < np.max(p1))
        s1 = np.logical_and(-phi_ > np.min(p0), -phi_ < np.max(p1))
        select = np.logical_or(s0, s1)
        
        va = np.argwhere(select)
        return va, select

    
    def plot_phi(self, y='lattice', family='200', frame=-1, idx=0, 
                         alpha=0.1, restrict_z=False, restrict_range = [70, 110]):
        """
        For a given lattice family (and frame) plots the variation in the 
        *resolved* lattice strain (or back stress) with respect to the angle 
        the planes make to the loading axis (phi). Can be restricted across 
        a smaller z_rot if required. N.b. rotations of grains defined as
        (x_rot, phi, z_rot).
        
        Parameters
        ----------
        y: str
            The data to plot on the y axis. This is typically lattice strain
            but it is also possible to plot wrt. back stress.
        family: str
            The lattice plane family to assess
        frame: int
            The frame to extract data from (default = 0).
        idx: int
            The compnent (referenced via an idx) of the defined data. Only 
            valid for back stress (for fcc, idx = 0-11)
        alpha: float
            Plotting data transparency
        restrict_z: bool
            Restrict data extraction/plotting across one angular range. Can be 
            used to normalise the amount of data wrt. phi
        restrict_range: [float, float]
            Range across which to limit z rotations.
            
        """
        lattice = self.lattice_strain
        
        y_ = {'lattice': lattice[family],
              'back stress': self.b_stress[idx]}[y]
        try:
            y_tensor = self.lattice_tensor[family]
            tens = True
        except KeyError:
            print('Tensor not available')
            tens=False
        
        if y == 'back stress':
            x = self.rot[1]
        else:
            x = self.lattice_phi[family]
            
        rot = self.lattice_rot[family]

        
        if restrict_z == True and y == 'lattice':
            r0, r1 = restrict_range
            t_z = rot[:, :, 2]* 180 / np.pi
            va = np.logical_and(t_z > r0, t_z < r1)
            vaf = np.zeros_like(rot[:, :, 2], dtype='bool')
            vaf[:, frame, :] += True
            va = np.logical_and(va, vaf)
        else:
            va = np.s_[:, frame]

        plt.plot(x[va].flatten(), y_[va].flatten(), '.', alpha=alpha)
        if y == 'lattice' and tens:
            st = strain_transformation(np.linspace(0, np.pi, 1001), *y_tensor[:, frame])
            plt.plot(np.linspace(0, np.pi, 1001), st, 'r')
        x = 'lattice rot (phi)' if y == 'lattice' else 'grain rot (phi)'
        plt.xlabel(x)
        plt.ylabel(y)
    
    
    def plot_grains(self, y='elastic', x='stress', x_mean=True, 
             y_mean=False, x_idx=1, y_idx=1, grain_idx=None, alpha=0.2,
             color='k', mcolor='r'):
        """
        The plot_grain method is very general plotting routing 
        and any grain (not lattice) specific vaues can be plotted on 
        either axis.

        - Define data to plot on either axis i.e. y='stress', x='strain'
        - Specify whether the data on given axis is the mean response of all grains
        - Where relevant, the index of that data must be specified 
        i.e. for y='stress', y_idx = 1 for sigma_yy
        
        While general a limited number of x, y combinations will, 
        unsurprisingly, not work.
        
        Parameters
        ----------
        y, x: str, str
            The data (label), either 'stress', 'strain', 'elastic' (strain),
            'back stress', 'rot', 'time', 'frame' to plot on x/y axis
        x_mean, y_mean: bool, bool
            Whether to take the mean (across all grains) of the data on the
            x/y axis
        x_idx, y_idx: int, int
            Component/orientation of the specified data to plot
            e.g. x='stress', idx=1 => sigma_xx
        grain_idx: [int, ...]
            List on grains (indices) to plot (if None, all grains plotted)
        alpha, color: float, str
            Plotting options for the grain specific lines
        mcolor:
            The color of the grain average (across x and y) line
        """
        # If necessary put grain_idx into list for fancy indexing
        if isinstance(grain_idx, int):
            grain_idx = [grain_idx,]
            
        # Time and frame can't be averaged
        if x in ['time', 'frame']:
            x_mean = False
        if y in ['time', 'frame']:
            y_mean = False
        
        # Data extraction
        x_ = self.extract_grains(data=x, idx=x_idx, grain_idx=grain_idx)
        y_ = self.extract_grains(data=y, idx=y_idx, grain_idx=grain_idx)

        # Saving x, y locations (?)
        csvfile = open('strain_grain.csv', 'w', newline='')
        obj = csv.writer(csvfile)
        for val in np.transpose(x_):
            obj.writerow(val)
        csvfile.close()
        
        csvfile = open('stress_grain.csv', 'w', newline='')
        obj = csv.writer(csvfile)
        for val in np.transpose(y_):
            obj.writerow(val)
        csvfile.close()
		
        # Calculate mean of arrays
        xm = np.nanmean(x_, axis=0) if x not in ['time', 'frame'] else x_
        ym = np.nanmean(y_, axis=0) if y not in ['time', 'frame'] else y_

        x__ = xm if x_mean else x_.T
        y__ = ym if y_mean else y_.T
        
        # Tinkering with axis labels
        x = '{} (idx={})'.format(x, x_idx) if x not in ['time', 'frame'] else x
        y = '{} (idx={})'.format(y, y_idx) if y not in ['time', 'frame'] else y
        x = 'mean {}'.format(x) if x_mean else x
        y = 'mean {}'.format(y) if y_mean else y
        
        # Plotting
        plt.plot(np.squeeze(x__), np.squeeze(y__), color=color, alpha=alpha)
        if (not y_mean or not x_mean) and (grain_idx == None or len(grain_idx) != 1):
            plt.plot(xm, ym, color=mcolor, label='Mean response')
            plt.legend()
        plt.ylabel(y)
        plt.xlabel(x)
        
    
    def plot_lattice(self, family='200', phi=0, window=10, lat_ax='x', 
                            ax2='stress', ax2_idx=1, ax2_mean=True,  
                            alpha=0.2, color='k', mcolor='r',
                            plot_select=True, phi_frame=0):
        
        """
        The lattice strains for a given family are plotted if they lie at (or 
        close to) an angle, phi (with the loading axis). The angular tolerance 
        / azimuthal window is defined by the user (window). For XRD, a window 
        of 10deg is often used.

        Parameters:
        -----------
        family: str
            The lattice plane family to assess
        phi: float
            Angle at which to extract the lattice plane strains
        window: float
            Azimuthal tolerance (absolute window width) for lattice data 
            extraction
        lat_ax: str
            Axis to plot the lattice data on, either 'x' or 'y'
        ax2: str
            The data to plot against the lattice strain. Either 'stress', 
            'strain', 'elastic' (strain), 'back stress'
        ax2_idx: int
            Component/orientation of the specified second axis data to plot
            e.g. ax2='stress', ax2_idx=1 => sigma_xx
        ax2_mean: bool
            Whether to take the mean (across all grains) of the data on the
            second axis
        alpha, color: float, str
            Plotting options for the grain specific lines
        mcolor:
            The color of the grain average (across x and y) line
        plot_select: bool
            If plot_select is True the individual lattice planes will be 
            plotted in addition to the mean result, when False just the mean
            response
        phi_frame: int
            The frame to define the grains that lie within the aimuthal
            window (default = 0).
        """
        
        ax2_mean = False if ax2 in ['time', 'frame'] else ax2_mean

        d = self.extract_grains(data=ax2, idx=ax2_idx, grain_idx=None)

        valid, select = self.extract_phi_idx(family=family, phi=phi,window=window, frame=phi_frame)
        if ax2 in ['time', 'frame']:
            d, dm = d, d
            
        else:
            
            d = np.nanmean(d, axis=0) if ax2_mean else d[valid[:,0]].T
            dm = d if ax2_mean else np.nanmean(d, axis=1)
            
        lattice = self.extract_lattice(family=family)
        
        lattice = lattice[valid[:,0], :, valid[:,1]].T
        
        x_ = lattice if lat_ax == 'x' else d
        y_ = lattice if lat_ax != 'x' else d

        
        assert np.sum(select) > 0, 'Phi window too small for {} - no grains/planes selected'.format(family)
        if plot_select:
            plt.plot(x_, y_, 'k', alpha=alpha)
            
        x_ = np.nanmean(lattice, axis=1) if lat_ax == 'x' else dm
        y_ = np.nanmean(lattice, axis=1) if lat_ax != 'x' else dm

        plt.plot(x_, y_, label=family, color=mcolor)
        
        
        ax2 = '{} (idx={})'.format(ax2, ax2_idx) if ax2 not in ['time', 'frame'] else ax2
        ax2 = ax2 if not ax2_mean else 'mean {}'.format(ax2)
        xlabel = ax2 if lat_ax != 'x' else 'lattice'
        ylabel = ax2 if lat_ax == 'x' else 'lattice'   
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
    def extract_lattice_map(self, family='200', az_bins=19):
        """
        Average the lattice strains data across a defined number of bins 
        (i.e.  azimuthally integrate), return 2D array of lattice strains against frame
        for the specified family.

        Parameters:
        -----------
        family: str
            The lattice plane family to assess
        az_bins: int
            Number of bins to extract lattice strains across

            
        Returns:
        --------
        bins: list
            List of the phi bins that data has been extracted at
        data: array
            Lattice strains for given family averaged across a user 
            defined (az_bins) number of azimuthally arrayed bins
        """
        
        phi_steps = az_bins + 1
        arr1 = np.moveaxis(self.lattice_strain[family], 1, 2)
        arr1 = arr1.reshape((-1, arr1.shape[-1]))
        
        arr2 = np.moveaxis(self.lattice_phi[family], 1, 2)
        arr2 = arr2.reshape((-1, arr2.shape[-1]))
        arr2[arr2 > np.pi/2] -= np.pi # -90 to 90
        
        bins = np.linspace(-90, 90, phi_steps)
        e_phi = np.nan * np.ones((phi_steps - 1, self.num_frames))
        
        for idx, i in enumerate(bins[:-1]):
            va = np.logical_and(arr2 < bins[idx + 1] * np.pi / 180, arr2 > bins[idx] * np.pi / 180)
            try:
                e_phi[idx] = np.sum(arr1 * va, axis=0) / np.nansum(va, axis=0)
            except ZeroDivisionError:
                pass
            
        return (bins[:-1]+bins[1:])/2, e_phi
        
    def plot_lattice_map(self, family='200', az_bins=19, ax2='time',
                                ax2_idx=1):

        """
        Plot 2D map of the azimtuhally arrayed lattice strains as a function of 
        second variable such as time or frame. Also works with macro stress
        or strain although obvious issues may arise if there is creep dwells.


        Parameters:
        -----------
        family: str
            The lattice plane family to assess
        az_bins: int
            Number of bins to extract lattice strains across
        ax2: str
            The data to plot against the lattice strain. Either 'stress', 
            'strain', 'elastic' (strain), 'back stress'
        ax2_idx: int
            Component/orientation of the specified second axis data to plot
            e.g. ax2='stress', ax2_idx=1 => sigma_xx

            
        Returns:
        --------
        bins: list
            List of the phi bins that data has been extracted at
        data: array
            Lattice strains for given family averaged across a user 
            defined (az_bins) number of azimuthally arrayed bins
        """
        
        bin_c, e_phi = self.extract_lattice_map(family=family, az_bins=az_bins)
        
        d = self.extract_grains(data=ax2, idx=ax2_idx, grain_idx=None)
        ax2_mean = False if ax2 in ['time', 'frame'] else True
        if ax2_mean:
            d = np.nanmean(d, axis=0)
        
        time, phi = np.meshgrid(d, bin_c)
        plt.contourf(time, phi, e_phi)
        plt.colorbar()

        ax2 = 'mean {} (idx={})'.format(ax2, ax2_idx) if ax2 not in ['time', 'frame'] else ax2

        plt.xlabel(ax2)
        plt.ylabel('phi (reflected at 0$^o$)')
            
    
    def plot_lattice_all(self, phi=0, window=10, lat_ax='x', ax2='stress', 
                         ax2_idx=1, ax2_mean=True, phi_frame=0):
        """
        The lattice strains for a ALL families are plotted if they lie at (or 
        close to) an angle, phi (with the loading axis). The angular tolerance 
        / azimuthal window is defined by the user (window). For XRD, a window 
        of 10deg is often used.

        Parameters:
        -----------
        phi: float
            Angle at which to extract the lattice plane strains
        window: float
            Azimuthal tolerance (absolute window width) for lattice data 
            extraction
        lat_ax: str
            Axis to plot the lattice data on, either 'x' or 'y'
        ax2: str
            The data to plot against the lattice strain. Either 'stress', 
            'strain', 'elastic' (strain), 'back stress'
        ax2_idx: int
            Component/orientation of the specified second axis data to plot
            e.g. ax2='stress', ax2_idx=1 => sigma_xx
        ax2_mean: bool
            Whether to take the mean (across all grains) of the data on the
            second axis
        phi_frame: int
            The frame to define the grains that lie within the aimuthal
            window (default = 0).
        """
        for family in self.lattice_list:
            try:
                self.plot_lattice(family=family, lat_ax=lat_ax, ax2=ax2, ax2_idx=ax2_idx, phi=phi, 
                         window=window, phi_frame=phi_frame, plot_select=False, mcolor=None, ax2_mean=ax2_mean)
            except AssertionError:
                print('Phi window too small for {} - no grains/planes selected'.format(family))
        plt.legend(self.lattice_list)
            

    def plot_back_lattice(self, family='200', phi=0, window=10,
                          back_ax='y', b_idx=1, ax2='stress', ax2_idx=1, 
                          alpha=0.2, color='k', mcolor='r',
                          plot_select=True, phi_frame=0):
        
        """
        Plot a component of back stress for a specified family of lattice 
        planes at a defined azimuthal angle. Plot against any other extracted 
        stress, strain, time etc. component.
        
        Parameters:
        -----------
        family: str
            The lattice plane family to assess
        phi: float
            Angle at which to extract the lattice plane strains
        window: float
            Azimuthal tolerance (absolute window width) for lattice data 
            extraction    
        back_ax: str
            Axis to plot the lattice data on, either 'x' or 'y'
        back_idx: int
            Component of the back stress to plot (for fcc 0-11)
        ax2: str
            The data to plot against the lattice strain. Either 'stress', 
            'strain', 'elastic' (strain), 'back stress'
        ax2_idx: int
            Component/orientation of the specified second axis data to plot
            e.g. ax2='stress', ax2_idx=1 => sigma_xx
        alpha, color: float, str
            Plotting options for the grain specific lines
        mcolor:
            The color of the grain average (across x and y) line
        plot_select: bool
            If plot_select is True the individual lattice planes will be 
            plotted in addition to the mean result, when False just the mean
            response
        phi_frame: int
            The frame to define the grains that lie within the aimuthal
            window (default = 0).
        """
        
        back = self.extract_grains(data='back stress', idx=b_idx, grain_idx=None)
        
        d = self.extract_grains(data=ax2, idx=ax2_idx, grain_idx=None)
        d = d if ax2 in ['time', 'frame'] else np.nanmean(d, axis=0)

        
        valid, select = self.extract_phi_idx(family=family, phi=phi,window=window, frame=phi_frame)
        
        # back = back[valid[:,0], :, valid[:,1]].T
        v = np.unique(valid[:,0])
        back = back[v, :].T
        
        x_ = back if back_ax == 'x' else d
        y_ = back if back_ax != 'x' else d
        
        assert np.sum(select) > 0, 'Phi window too small for {} - no grains/planes selected'.format(family)
        if plot_select:
            plt.plot(x_, y_, 'k', alpha=alpha)
            
        ax2 = 'mean {} (idx={})'.format(ax2, ax2_idx) if ax2 not in ['time', 'frame'] else ax2
        
        xlabel = ax2 if back_ax != 'x' else 'back stress'
        ylabel = ax2 if back_ax == 'x' else 'back stress'   
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
            
    def plot_active_slip(self, family='200', phi=0, window=10, 
                         back_ax='y', b_active=2, ax2='stress', ax2_idx=1, 
                         alpha=0.2, color='k', mcolor='r',
                         plot_select=True, phi_frame=0):
        
        """
        Plot the number of active slip systems for every plane for a specified 
        family of lattice planes at a defined azimuthal angle (angle wrt y axis). 
        Plotting is a function of time, frame, stress strain etc. The 
        activation of a slip system is taken to occur when the absolute back
        stress associated with that system (i.e. back stress component) rises 
        above a user define value
        
        Parameters:
        -----------
        family: str
            The lattice plane family to assess
        phi: float
            Angle at which to extract the lattice plane strains
        window: float
            Azimuthal tolerance (absolute window width) for lattice data 
            extraction    
        back_ax: str
            Axis to plot the lattice data on, either 'x' or 'y'
        b_active: int
            Component of the back stress to plot (for fcc 0-11)
        ax2: str
            The data to plot against the lattice strain. Either 'stress', 
            'strain', 'elastic' (strain), 'back stress'
        ax2_idx: int
            Component/orientation of the specified second axis data to plot
            e.g. ax2='stress', ax2_idx=1 => sigma_xx
        alpha, color: float, str
            Plotting options for the grain specific lines
        mcolor:
            The color of the grain average (across x and y) line
        plot_select: bool
            If plot_select is True the individual lattice planes will be 
            plotted in addition to the mean result, when False just the mean
            response
        phi_frame: int
            The frame to define the grains that lie within the aimuthal
            window (default = 0).
        """
        
        back = self.extract_grains(data='back stress', idx=None, grain_idx=None)
        back_bool = np.abs(back) > b_active
        
        d = self.extract_grains(data=ax2, idx=ax2_idx, grain_idx=None)
        d = d if ax2 in ['time', 'frame'] else np.nanmean(d, axis=0)

        
        valid, select = self.extract_phi_idx(family=family, phi=phi,window=window, frame=phi_frame)
        
        # back = back[valid[:,0], :, valid[:,1]].T
        v = np.unique(valid[:,0])
        back_active = np.sum(back_bool, axis=0)[v, :].T
        
        x_ = back_active if back_ax == 'x' else d
        y_ = back_active if back_ax != 'x' else d
        
        assert np.sum(select) > 0, 'Phi window too small for {} - no grains/planes selected'.format(family)
        if plot_select:
            plt.plot(x_, y_, 'k', alpha=alpha)
            
        x_ = np.nanmean(back_active, axis=1) if back_ax == 'x' else d
        y_ = np.nanmean(back_active, axis=1) if back_ax != 'x' else d

        plt.plot(x_, y_, label=family, color=mcolor)
        
        ax2 = 'mean {} (idx={})'.format(ax2, ax2_idx) if ax2 not in ['time', 'frame'] else ax2
        xlabel = ax2 if back_ax != 'x' else 'Active slip systems'
        ylabel = ax2 if back_ax == 'x' else 'Active slip systems'   
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
    def plot_active_slip_all(self, phi=0, window=10, back_ax='y', b_active = 2,
                          ax2='stress', ax2_idx=1, phi_frame=0):
        """
        Plot the plane averaged number of active slip systems for all families 
        of lattice planes at a defined azimuthal angle (angle wrt y axis). 
        Plotting is a function of time, frame, stress strain etc. The 
        activation of a slip system is taken to occur when the absolute back
        stress associated with that system (i.e. back stress component) rises 
        above a user define value
        
        Parameters:
        -----------
        phi: float
            Angle at which to extract the lattice plane strains
        window: float
            Azimuthal tolerance (absolute window width) for lattice data 
            extraction    
        back_ax: str
            Axis to plot the lattice data on, either 'x' or 'y'
        b_active: int
            Component of the back stress to plot (for fcc 0-11)
        ax2: str
            The data to plot against the lattice strain. Either 'stress', 
            'strain', 'elastic' (strain), 'back stress'
        ax2_idx: int
            Component/orientation of the specified second axis data to plot
            e.g. ax2='stress', ax2_idx=1 => sigma_xx
        phi_frame: int
            The frame to define the grains that lie within the aimuthal
            window (default = 0).
        """
        for family in self.lattice_list:
            try:
                self.plot_active_slip(family=family, back_ax=back_ax, ax2=ax2, ax2_idx=ax2_idx, phi=phi, 
                         window=window, frame=phi_frame, plot_select=False, mcolor=None)
            except AssertionError:
                print('Phi window too small for {} - no grains/planes selected'.format(family))
        plt.legend(self.lattice_list)


