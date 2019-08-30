import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

from odbAccess import *
import pickle
from abaqusConstants import *
import time
from cpex.data_io import cpex_to_hdf5
from cpex.lattice import Lattice

    

class ScrapeODB(Lattice):
    def __init__(self, fpath, N=12):
        """
        Args:
            fpath (str): Path to an odb file
            N (int): Number of slip systems
        """
        self.fpath = fpath
        self.odb = openOdb(path=self.fpath)    
        self.frames= self.odb.steps['Loading'].frames
        self.instances = self.odb.rootAssembly.instances['DREAM-1'.upper()]
        self.num_frames = len(self.frames)
        self.N = N
        
        ### Need to calc number of grains
        self.num_grains = 1 # HOW?

        # Create empy arrays to store data
        d_shape = (len(self.num_grains), len(self.num_frames))
        
        s = np.empty((6,) + d_shape)
        e = np.empty((6,) + d_shape)
        lat = np.empty((6,) + d_shape)
        rot = np.empty((3,) + d_shape)
        dims = np.empty((3,) + d_shape)
        v = np.empty(d_shape)

        
        # Open each frame in turn and extract data
        for fidx, frame in enumerate(self.frames):            
            sc = scrape_frame(frame, self.num_grains, self.N)

            s[:, :, fidx] = sc[0]
            e[:, :, fidx] = sc[1]
            lat[:, :, fidx] = sc[2]
            rot[:, :, fidx] = sc[3]
            dims[:, :, fidx] = sc[4]
            v[:, fidx] = sc[5]

        self.s = s
        self.e = e
        self.lat = lat        
        self.dims = dims
        self.rot = rot
        self.v = v  
        


def scrape_frame(frame, num_grains, N):
    
    lat_SDV_nums = range((N * 12 + 4), (N * 12 + 4) + 6 + 1)
    rot_SDV_nums = range((N * 34 + 3), (N * 34 + 3) + 3 + 1)
    
    s = np.empty((6,len(num_grains)))
    e = np.empty((6,len(num_grains)))  
    
    lat = np.empty((6,len(num_grains)))
    rot = np.empty((3,len(num_grains)))
    
    v = np.empty((len(num_grains)))
    dims = np.empty((3,len(num_grains)))
    
    data_fo = frame.fieldOutputs

    # Abaqus
    stress_ = data_fo['S']
    strain_ = data_fo['LE']
    Volume_ = data_fo['IVOL']
    coords_ = data_fo['COORD']
    
    # SDVs
    latticeSDV_ = [data_fo['SDV{}'.format(i)] for i in lat_SDV_nums]
    rotSDV_ = [data_fo['SDV{}'.format(i)] for i in rot_SDV_nums]
    
    for idx, grain in enumerate(range(1, num_grains + 1)):
        
        myInstance2=instances.elementSets[grain] # Select grain
        numElements=len(myInstance2.elements) # Find all elements in grain
        
        # Macro, elastic-plastic response
        stress = stress_.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values
        strain = strain_.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values
        Volume = Volume_.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values
        coords = coords_.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values
        
        # Lattice results exx, eyy, ezz, exy, exz, eyz
        latSDV = [i.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values for i in latticeSDV_]
        rotSDV = [i.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values for i in rotSDV_]
        
        # Iterate over total number of elements and sum the strain/unit volume
        s_v = np.zeros((6,))
        e_v = np.zeros((6,))
        lat_v = np.zeros((6,))
        rot_v = np.zeros((3,))
        coords_v = np.zeros((3,))
        vv = 0

        for ip in range(numElements*8):
            vv_ = Volume[ip].data
            
            s_v += np.array([stress[ip].data[i] * vv_ for i in range(1, 7)])
            e_v += np.array([strain[ip].data[i] * vv_ for i in range(1, 7)])
            coords_v += np.array([coords[ip].data[i] * vv_ for i in range(1, 4)])

            lat_v += np.array([i[ip].data * vv_ for idx, i in enumerate(latSDV)])
            rot_v += np.array([i[ip].data * vv_ for idx, i in enumerate(rotSDV)])
            
            vv += Volume[ip].data
        

        # Unpack values
        s[:,idx] = s_v / vv
        e[:,idx] = e_v  / vv
        lat[:,idx] = lat_v / vv            
        dims[:,idx] = coords_v / vv
        rot[:,idx] = rot_v / vv
        v[idx] = vv
        
    return  (s, e, lat, dims, rot, v)

