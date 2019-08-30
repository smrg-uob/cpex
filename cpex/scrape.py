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
        
        s = [np.empty(d_shape),]*6
        e = [np.empty(d_shape),]*6
        lat = [np.empty(d_shape),]*6
        rot = [np.empty(d_shape),]*3
        dims = [np.empty(d_shape),]*3
        v = np.empty(d_shape)

        
        # Open each frame in turn and extract data
        for fidx, frame in enumerate(self.frames):
            
            sc = scrape_frame(frame, self.num_grains, self.N)
            
            for i in range(6):
                s[i][fidx, :] = sc[0][i]
                e[i][fidx, :] = sc[1][i]
                lat[i][fidx, :] = sc[2][i]
                
            for i in range(3):
                dims[i][fidx, :] = sc[3][i]
                rot[i][fidx, :] = sc[4][i]
            
            v[fidx, :] = sc[5]

            
            
        self.sxx, self.syy, self.szz, self.sxy, self.syz, self.sxz = s
        self.exx, self.eyy, self.ezz, self.exy, self.eyz, self.exz = e
        self.latxx, self.latyy, self.latzz, self.latxy, self.latyz, self.latxz = lat        
        self.x, self.y, self.z = dims
        self.phi1, self.phi2, self.phi3, rot
        self.v = v  
        


    

def scrape_frame(frame, num_grains, N):
    
    lat_SDV_nums = range((N * 12 + 4), (N * 12 + 4) + 6)
    rot_SDV_nums = range((N * 34 + 3), (N * 34 + 3) + 3)
    
    s = [np.empty((len(num_grains))),]*6    
    e = [np.empty((len(num_grains))),]*6  
    
    lat = [np.empty((len(num_grains))),]*6
    rot = [np.empty((len(num_grains))),]*3
    
    v = np.empty((len(num_grains)))
    dims = [np.empty((len(num_grains))),]*3
    
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
        latticeSDV = [i.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values for i in latticeSDV_]
        rotSDV = [i.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values for i in rotSDV_]
        
        # Iterate over total number of elements and sum the strain/unit volume
        s_v = [0, ] * 6
        e_v = [0, ] * 6
        lattice_v = [0, ] * 6
        rot_v = [0, ] * 3
        coords_v = [0, ] * 3
        vv = 0

        for ip in range(numElements*8):
            vv_ = Volume[ip].data
            
            s_v = [s_v[i-1] + stress[ip].data[i] * vv_ for i in range(1, 7)]
            e_v = [e_v[i-1] + strain[ip].data[i] * vv_ for i in range(1, 7)]
            coords_v = [coords_v[i-1] + coords[ip].data[i] * vv_ for i in range(1, 4)]
            
            
            lattice_v = [lattice_v[idx] + i[ip].data * vv_ for idx, i in enumerate(latticeSDV)]
            rot_v = [rot_v[idx] + i[ip].data * vv_ for idx, i in enumerate(rotSDV)]
            
            vv += Volume[ip].data
        
        s_v = [s_v[i] / vv for i in range(6)]
        e_v = [e_v[i] / vv for i in range(6)]
        lattice_v = [lattice_v[i] / vv for i in range(6)]
        rot_v = [rot_v[i] / vv for i in range(3)]
        coords_v = [coords_v[i] / vv for i in range(3)]
        vv = 0

        # Unpack values
        for i in range(6):
            s[i][idx] = s_v[i]
            e[i][idx] = e_v[i]
            lat[i][idx] = lattice_v[i]
            
        for i in range(3):
            dims[i][idx] = coords_v[i]
            rot[i][idx] = rot_v[i]
        
        v[idx] = vv
        
    return  (s, e, lat, dims, rot, v)

