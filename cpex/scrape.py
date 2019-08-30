import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

from odbAccess import *
import pickle
from abaqusConstants import *
import time
from cpex.data_io import cpex_to_hdf5


class Lattice():
    def extract_lattice(self):
        pass
    
    def calc_lattice_tensor(self):
        pass
    

class ScrapeODB(Lattice):
    def __init__(self, fpath):
        self.fpath = fpath
        self.odb = openOdb(path=self.fpath)    
        self.frames= self.odb.steps['Loading'].frames
        self.instances = self.odb.rootAssembly.instances['DREAM-1'.upper()]
        self.num_frames = len(self.frames) 
        
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
            
            sc = scrape_frame(frame, self.num_grains)
            
            sxx[fidx,:], syy[fidx,:], sxy[fidx,:] = sc[0], sc[1], sc[2]
            exx[fidx,:], eyy[fidx,:], exy[fidx,:] = sc[3], sc[4], sc[5]
            v[fidx,:] = sc[6]
            eSDV149[fidx,:] = sc[7]
            ro1[fidx,:], ro2[fidx,:], ro3[fidx,:] = sc[8], sc[9], sc[10]
            
        self.sxx, self.syy, self.sxy = sxx, syy, sxy
        self.exx, self.eyy, self.exy = sxx, syy, sxy
        self.v = v
        self.eSDV149 = eSDV149
        self.ro1, self.ro2, self.ro3 = ro1, ro2, ro3
        

    def save_to_hdf5(self, fpath, overwrite=False):
        """ Save data back to hdf5 file format.

        Saves analyzed information and the raw, scraped data.

        Args:
            fpath (str): Abs. path for new file
            overwrite (bool): Overwrite file if it already exists
        """
        cpex_to_hdf5(fpath, self, overwrite)        
        


    

def scrape_frame(frame, num_grains):
    
    sxx, syy, sxy = [np.empty((len(num_grains))),]*3    
    exx, eyy, exy = [np.empty((len(num_grains))),]*3  
    ro1, ro2, ro3 = [np.empty((len(num_grains))),]*3
    
    v = np.empty((len(num_grains)))
    eSDV149 = np.empty((len(num_grains))) # what is this?
    
    data_fo = frame.fieldOutputs

    stress_ = data_fo['S']
    strain_ = data_fo['LE']
    Volume_ = data_fo['IVOL']
    SDV149_ = data_fo['SDV149']
    
    for idx, grain in enumerate(range(1, num_grains + 1)):
        
        myInstance2=instances.elementSets[grain] # Select grain
        numElements=len(myInstance2.elements) # Find all elements in grain
        
        stress = stress_.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values
        strain = strain_.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values
        Volume = Volume_.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values
        SDV149 = SDV149_.getSubset(region=myInstance2,position=INTEGRATION_POINT,elementType='C3D8').values
        
        
        # Iterate over total number of elements and sum the strain/unit volume
        syy_vyy, eyy_vyy, vyy, eSDV149_vyy = 0, 0, 0, 0

        for ip in range(numElements*8):
            vyy_ = Volume[ip].data
            
            syy_vyy += stress[ip].data[1] * vyy_
            eyy_vyy += strain[ip].data[1] * vyy_
            eSDV149_vyy += SDV149[ip].data * vyy_
            vyy += Volume[ip].data
            
        syy[idx] = syy_vyy
        eyy[idx] = eyy_vyy
        v[idx] = vyy
        eSDV149[idx] = eSDV149_vyy
        
    return  sxx, syy, sxy, exx, eyy, exy, v, eSDV149, ro1, ro2, ro3