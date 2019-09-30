# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:48:44 2019

@author: casimp
"""
import os
import numpy as np
from scipy.optimize import curve_fit
from collections import OrderedDict

from cpex.nvec import nvec_extract
from cpex.transformation import trans_matrix, strain_transformation
from cprex.extract import Extract


class Load(Extract):
    
    def __init__(self, fpath, calc=True, lattice_list = ['111', '200', '220', '311']):
        """
        Initialises the cpex Load class. Takes in a .npz file, pulls out the 
        data and then performs initial calculations to find and resolve the
        resolves strains for all lattice planes for lattice familieis specified
        in the lattice list. 
        
        Parameters
        ----------
        fpath: str
            Path to the .npz file
        calc: bool
            Auto run/call the lattice plane extraction and strain resolve.
            Also fit an in-plane (x-y) tensor through the phi resolve lattice
            strain data.
        lattice_list: [str, str, str...]
            List of all the lattice plane families of interest (for fcc this 
            normally would be ['111', '200', '220', '311'])
        
        """
        data = np.load(fpath)
        self.e = data['e']
        self.s = data['s']
        self.elastic = data['lat']
        self.dims = data['dims']
        self.rot = data['rot']
        self.v = data['v']
        self.N = data['N']
        self.num_frames = data['num_frames']
        self.num_grains = data['num_grains']
        try:
            self.t = data['time']
        except KeyError:
            print('time not saved in file, creating zero array')
            self.t = np.zeros((self.num_frames, ))
            
        try:
            self.b_stress = data['b_stress']
        except KeyError:
            print('back stress not saved in file, creating zero array')
            d_shape = (self.num_grains, self.num_frames)
            self.b_stress = np.zeros((12,) + d_shape)
            
        self.rot[:, :,0] = self.rot[:, :,1]
        self.lattice_list = lattice_list
        self.lattice_nvecs = [nvec_extract(*[int(i) for i in hkl]) for hkl in self.lattice_list]
        
        if calc:
            print('Calculating lattice rotations and strains...')
            self.calc_lattice_rot()
            self.calc_lattice_strain()
            self.calc_lattice_tensor()
            
    
    def calc_lattice_rot(self):
        """
        Extracts all angles for FCC grain families. 
         Needs to be generalised to all structures.
        """
        r0, r1, r2 = self.rot[0], self.rot[1], self.rot[2]
        total_rot = trans_matrix(r0, r1, r2)
        total_rot = np.transpose(total_rot, (0, 1, 3, 2))
   
    
        angles = []
        for nvec in self.lattice_nvecs:
    
            rot_v1 = np.matmul(total_rot, nvec.T) # total rot matrix
            yax=np.array([[0,1,0]]).T
            angle = np.arccos(rot_v1/(np.linalg.norm(yax)*np.linalg.norm(rot_v1, axis=-2))[:, :, np.newaxis,:] )
            angles.append(angle)
        
        self.lattice_rot = OrderedDict(zip(self.lattice_list, angles))
        self.lattice_phi = OrderedDict(zip(self.lattice_list, [i[..., 1, :] for i in angles]))


    def calc_lattice_strain(self):
    
        ens = []
        for nvec in self.lattice_nvecs:
            exx, eyy, ezz, exy, exz, eyz =self.elastic
            eT = np.array([[exx, exy, exz],
                              [-exy, eyy, eyz],
                              [-exz, -eyz, ezz]])
    
            eT = np.moveaxis(np.moveaxis(eT, 0, -1), 0, -1) #[:, :, np.newaxis, :, :]
            r0, r1, r2 = self.rot[0], self.rot[1], self.rot[2]
            
            
            r = trans_matrix(r0, r1, r2)
            
            eTrot = np.matmul(np.matmul(r, eT), np.transpose(r, (0,1,3,2)))
            eTrot = eTrot[..., np.newaxis, :, :]
            nvec = nvec[:, np.newaxis, :]
            
            en = nvec@eTrot@np.transpose(nvec, (0, 2, 1))
            en = en[:, :, :, 0, 0]
            
            ens.append(en)
        
        self.lattice_strain = dict(zip(self.lattice_list, ens))
             
        
        
    def calc_lattice_tensor(self):
        
        tensors, tensors_err = [], []

        for e_lat, phi, rot in zip(self.lattice_strain.values(),
                                   self.lattice_phi.values(),
                                   self.lattice_rot.values()):
        
            e_tensor, e_tensor_err = np.zeros((3, self.num_frames)), np.zeros((3, self.num_frames))
            for idx in range(self.num_frames):
                popt, pcov = curve_fit(strain_transformation, phi[:,idx].flatten(), e_lat[:,idx].flatten(), p0=[0, 0, 0])
                e_tensor[:, idx] = popt
                e_tensor_err[:, idx] = np.sqrt(np.diag(pcov))
            
            tensors.append(e_tensor)
            tensors_err.append(e_tensor_err)
            
        self.lattice_tensor = OrderedDict(zip(self.lattice_list, tensors))
        self.lattice_tensor_err = OrderedDict(zip(self.lattice_list, tensors_err))
        
        
if __name__ == '__main__':
    folder = os.path.join(os.path.dirname(__file__), r'data') # should be sub [0]
    fpath = os.path.join(folder, r'test_cpex.npz')
    data = Load(fpath)
    data.plot(x='strain', y='stress')

