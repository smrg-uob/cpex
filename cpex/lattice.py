# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:48:44 2019

@author: cs17809
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class Load():
    
    def __init__(self, fpath):
        data = np.load(fpath)
        self.e = data['e']
        self.s = data['s']
        self.lat = data['lat']
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
        self.rot[:, :,0] = self.rot[:, :,1]
        
        
    
    def plot(self, y='stress', x='frame', idx=1, alpha=0.2):
        
        #Need tp average over volume.
        
        dy = {'strain':self.e,
             'stress':self.s,
             'lattice':self.lat,
             'rot':self.rot - self.rot[:,:, 0][:, :, None]}
        
        dx = {'time':self.t,
             'stress':np.nanmean(self.s[idx], axis=0),
             'strain':np.nanmean(self.e[idx], axis=0),
             'frame':np.arange(self.num_frames),
             'lattice':np.nanmean(self.lat[idx], axis=0)}
        
        
        plt.plot(dx[x], dy[y][idx].T, color='k', alpha=alpha)
        plt.plot(dx[x],np.nanmean(dy[y][idx], axis=0), color='r')
        plt.ylabel(y)
        plt.xlabel(x)


    def plot_subset(self, grain_list, y='stress', x='frame', idx=1, alpha=0.2, 
                    random=None):
        
        
        if random != None:
            grain_list = np.random.choice(self.num_grains, random)
        else:
            grain_list = [i-1 for i in grain_list]
        
        dy = {'strain':self.e[:,grain_list],
             'stress':self.s[:,grain_list],
             'lattice':self.lat[:,grain_list],
             'rot':(self.rot - self.rot[:,:, 0][:, :, None])[:,grain_list]}
        
        dx = {'time':self.t,
             'stress':np.nanmean(self.s[idx], axis=0),
             'strain':np.nanmean(self.e[idx], axis=0),
             'frame':np.arange(self.num_frames),
             'lattice':np.nanmean(self.lat[idx], axis=0)}
        
        
        plt.plot(dx[x], dy[y][idx].T, color='k', alpha=alpha)
        plt.plot(dx[x],np.nanmean(dy[y][idx], axis=0), color='r')
        plt.ylabel(y)
        plt.xlabel(x)
    


    def extract_lattice_rot(self):
        """
        Extracts all angles for FCC grain families. 
         Needs to be generalised to all structures.
        """
        r0, r1, r2 = self.rot[0], self.rot[1], self.rot[2]
        total_rot = trans_matrix(r0, r1, r2)
        total_rot = np.transpose(total_rot, (0, 1, 3, 2))
        
        nvec_111 = (1 / np.sqrt(3)) * np.array([[1, 1, 1],
                                              [1, 1, -1],
                                              [1, -1, 1],
                                              [-1, 1, 1]])
    
        nvec_200 = np.array([[1, 0, 0],
                             [0, 1, 0],
                             [0, 0, 1]])
    
        nvec_220 = (2 / np.sqrt(8)) * np.array([[1, 1, 0],
                                                  [1, 0, 1],
                                                  [0, 1, 1],
                                                  [1, -1, 0],
                                                  [1, 0, -1],
                                                  [0, -1, 1]])
        
        nvec_311 = (1 / np.sqrt(11)) * np.array([[3, 1, 1],
                                          [1, 3, 1],
                                          [1, 1, 3],
                                          [-3, 1, 1],
                                          [1, -3, 1],
                                          [1, 1, -3]])
    



        angles = []
        for nvec in [nvec_111, nvec_200, nvec_220, nvec_311]:
    
            rot_v1 = np.matmul(total_rot, nvec.T) # total rot matrix
            yax=np.array([[0,1,0]]).T
            angle = np.arccos(rot_v1/(np.linalg.norm(yax)*np.linalg.norm(rot_v1, axis=-2))[:, :, np.newaxis,:] )
            angles.append(angle)
        
        self.rot_111, self.rot_200, self.rot_220, self.rot_311 = angles
        self.phi_111, self.phi_200, self.phi_220, self.phi_311 = self.rot_111[..., 1, :], self.rot_200[..., 1, :], self.rot_220[..., 1, :], self.rot_311[..., 1, :]
        


    def calc_lattice_strain(self):
        #r0, r1, r2 = self.phi_200[:, :, 0], self.phi_200[:, :, 1], self.phi_200[:, :, 2]

        nvec_111 = (1 / np.sqrt(3)) * np.array([[1, 1, 1],
                          [1, 1, -1],
                          [1, -1, 1],
                          [-1, 1, 1]])

        nvec_200 = np.array([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]])

        nvec_220 = (2 / np.sqrt(8)) * np.array([[1, 1, 0],
                                                  [1, 0, 1],
                                                  [0, 1, 1],
                                                  [1, -1, 0],
                                                  [1, 0, -1],
                                                  [0, -1, 1]])
        
        nvec_311 = (1 / np.sqrt(11)) * np.array([[3, 1, 1],
                                          [1, 3, 1],
                                          [1, 1, 3],
                                          [-3, 1, 1],
                                          [1, -3, 1],
                                          [1, 1, -3]])

        ens = []
        for nvec in [nvec_111, nvec_200, nvec_220, nvec_311]:
            exx, eyy, ezz, exy, exz, eyz =self.lat
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
            
        self.e_111, self.e_200, self.e_220, self.e_311 = ens
             
        
        
    def calc_lattice_tensor(self):
        
        tensors, tensors_err = [], []
        for e_lat, phi in zip([self.e_111, self.e_200, self.e_220, self.e_311],
                              [self.phi_111, self.phi_200, self.phi_220, self.phi_311]):
            e_tensor, e_tensor_err = np.zeros((3, self.num_frames)), np.zeros((3, self.num_frames))
            for idx in range(self.num_frames):
                popt, pcov = curve_fit(strain_transformation, phi[:,idx].flatten(), e_lat[:,idx].flatten(), p0=[0, 0, 0])
                e_tensor[:, idx] = popt
                e_tensor_err[:, idx] = np.sqrt(np.diag(pcov))
            
            tensors.append(e_tensor)
            tensors_err.append(e_tensor_err)
            
        self.e_111_tensor, self.e_200_tensor, self.e_220_tensor, self.e_311_tensor = tensors
        self.e_111_tensor_err, self.e_200_tensor_err, self.e_220_tensor_err, self.e_311_tensor_err = tensors_err
        
        
    def plot_lattice_phi(self, lattice='200', frame=-1, alpha=0.1):
        
        y = {'111':self.e_111,
             '200':self.e_200,
             '311':self.e_311,
             '220':self.e_220}
        
        y_tensor = {'111':self.e_111_tensor,
                    '200':self.e_200_tensor,
                    '311':self.e_311_tensor,
                    '220':self.e_220_tensor}
        
        x = {'111':self.phi_111,
             '200':self.phi_200,
             '311':self.phi_311,
             '220':self.phi_220}
        
        
        plt.plot(x[lattice][:, frame, :].flatten(), y[lattice][:, frame, :].flatten(), '.', alpha=alpha)
        plt.plot(np.linspace(0, np.pi, 1001), strain_transformation(np.linspace(0, np.pi, 1001), *y_tensor[lattice][:, frame]), 'r')
        
    def plot_lattice_select(self, lattice='200', y='stress', idx=1, 
                            alpha=0.1, window=10, plot_select=True):
        
        
        y_ = {'strain':self.e,
             'stress':self.s}
        
        phi_ = {'111':self.phi_111,
             '200':self.phi_200,
             '311':self.phi_311,
             '220':self.phi_220} 
        
        phi = phi_[lattice]
        
        if idx == 1:
            va = np.argwhere(np.abs(phi[:, 94]*180/np.pi - 90)> 90 - window/2)
        else:
            va = np.argwhere(np.abs(phi[:, 94]*180/np.pi - 90) < window/2)
            
        
        s = np.nanmean(y_[y], axis=1)[1]
        if plot_select:
            plt.plot(self.lat[idx].T[:, va[:, 0]], s, 'k', alpha=alpha)
        plt.plot(np.nanmean(self.lat[idx].T[:, va[:, 0]], axis=1), s, label=lattice)
    
    def plot_lattice_all(self, y='stress', idx=1, alpha=0.1, window=10):
        for lattice in ['111', '311', '220', '200']:
            self.plot_lattice_select(lattice=lattice, y=y, idx=idx, 
                            alpha=alpha, window=window, plot_select=False)
        


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
    
if __name__ == '__main__':
    folder = os.path.join(os.path.dirname(__file__), r'data') # should be sub [0]
    fpath = os.path.join(folder, r'test_cpex.npz')
    data = Load(fpath)
    data.plot(x='strain', y='stress')
    # data.plot(x='lattice', y='stress')
