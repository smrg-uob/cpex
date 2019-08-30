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
        sxx, syy, sxy = [np.empty(d_shape),]*3 
        exx, eyy, exy = [np.empty(d_shape),]*3 
        ro1, ro2, ro3 = [np.empty(d_shape),]*3 

        v = np.empty(d_shape)
        eSDV149 = np.empty(d_shape) # what is this?
        
        # Open each frame in turn and extract data
        for fidx, frame in enumerate(self.frames):
            
            sc = scrape_frame(frame, self.num_grains, self.N)
            
            sxx[fidx,:], syy[fidx,:], szz[fidx,:] = sc[0], sc[1], sc[2]
            sxy[fidx,:], syz[fidx,:], sxz[fidx,:] = sc[3], sc[4], sc[5]
            
     ############ Finish       
#            exx[fidx,:], eyy[fidx,:], exy[fidx,:] = sc[3], sc[4], sc[5]
#            v[fidx,:] = sc[6]
#            eSDV149[fidx,:] = sc[7]
#            ro1[fidx,:], ro2[fidx,:], ro3[fidx,:] = sc[8], sc[9], sc[10]
            
        self.sxx, self.syy, self.sxy = sxx, syy, sxy
        self.exx, self.eyy, self.exy = sxx, syy, sxy
        self.v = v
        self.eSDV149 = eSDV149
        self.ro1, self.ro2, self.ro3 = ro1, ro2, ro3     
        


    

def scrape_frame(frame, num_grains, N):
    
    lat_SDV_nums = range((N * 12 + 4), (N * 12 + 4) + 6)
    rot_SDV_nums = range((N * 34 + 3), (N * 34 + 3) + 3)
    
    sxx, syy, szz, sxy, syz, sxz = [np.empty((len(num_grains))),]*6    
    exx, eyy, ezz, exy, eyz, exz = [np.empty((len(num_grains))),]*6  
    
    latxx, latyy, latzz, latxy, latyz, latxz = [np.empty((len(num_grains))),]*6
    phi1, phi2, phi3 = [np.empty((len(num_grains))),]*3
    
    v = np.empty((len(num_grains)))
    x, y, z = [np.empty((len(num_grains))),]*3
    
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

           
        sxx[idx], syy[idx], szz[idx] = s_v[0], s_v[1], s_v[2]
        sxy[idx], syz[idx], sxz[idx] = s_v[3], s_v[4], s_v[5]
        
        exx[idx], eyy[idx], ezz[idx] = e_v[0], e_v[1], e_v[2]
        exy[idx], eyz[idx], exz[idx] = e_v[3], e_v[4], e_v[5]
        
        latxx[idx], latyy[idx], latzz[idx] = lattice_v[0], lattice_v[1], lattice_v[2]
        latxy[idx], latyz[idx], latxz[idx] = lattice_v[3], lattice_v[4], lattice_v[5]
        
        x[idx], y[idx], z[idx] = coords_v
        phi1[idx], phi2[idx], phi3[idx] = rot_v
        
        v[idx] = vv
        
    return  (sxx, syy, szz, sxy, syz, sxz, exx, eyy, ezz, exy, eyz, exz, 
             latxx, latyy, latzz, latxy, latyz, latxz, 
             v, x, y, z, phi1, phi2, phi3)

