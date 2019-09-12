# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 19:34:24 2019

@author: cs17809
"""
import numpy as np
from scipy.optimize import curve_fit

def trans_matrix(r0, r1, r2):
    """ Creates a multidimensional transofrmation array from 3 rotational arrays:

    Args:
        phi (float, ndarray): Azimuthal angle (rad)
        p[0] (float, ndarray): e_xx
        p[1] (float, ndarray): e_yy
        p[2] (float, ndarray): e_xy

    Returns:
        float, ndarray: Stress/strain wrt. azimuthal angle(s)
    """            
    zrot = np.array([np.stack([np.cos(r0), np.sin(r0), np.zeros_like(r0)]),
           np.stack([-np.sin(r0), np.cos(r0), np.zeros_like(r0)]), 
           np.stack([np.zeros_like(r0),np.zeros_like(r0),np.ones_like(r0)])])
    
    xrot = np.array([np.stack([np.ones_like(r1),np.zeros_like(r1),np.zeros_like(r1)]),
           np.stack([np.zeros_like(r1), np.cos(r1),np.sin(r1)]),
           np.stack([np.zeros_like(r1),-np.sin(r1), np.cos(r1)])])
    
    zrot2 = np.array([np.stack([np.cos(r2),np.sin(r2),np.zeros_like(r2)]),
            np.stack([-np.sin(r2),np.cos(r2),np.zeros_like(r2)]),
            np.stack([np.zeros_like(r2),np.zeros_like(r2),np.ones_like(r2)])])
    
    zrot = np.moveaxis(np.moveaxis(zrot, 0, -1), 0, -1)#.swapaxes(-1, -2)
    xrot = np.moveaxis(np.moveaxis(xrot, 0, -1), 0, -1)# .swapaxes(-1, -2)
    zrot2 = np.moveaxis(np.moveaxis(zrot2, 0, -1), 0, -1)#.swapaxes(-1, -2)
            
    R = np.matmul(np.matmul(zrot2,xrot),zrot)
    
    
    return R
            

def strain_transformation(phi, *p):
    """ Stress/strain (normal) transformation

    Args:
        phi (float, ndarray): Azimuthal angle (rad)
        p[0] (float, ndarray): e_xx
        p[1] (float, ndarray): e_yy
        p[2] (float, ndarray): e_xy

    Returns:
        float, ndarray: Stress/strain wrt. azimuthal angle(s)
    """
    average = (p[0] + p[1]) / 2
    radius = (p[0] - p[1]) / 2
    return average + np.cos(2 * phi) * radius + p[2] * np.sin(2 * phi)
    
