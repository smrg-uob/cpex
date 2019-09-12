# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 19:32:23 2019

@author: cs17809
"""
import numpy as np

def nvec_extract(h=1, k=1, l=1):
    h, k, l = np.abs(h), np.abs(k), np.abs(l)
    hkl = np.array([h, k, l])
    if h == k == l:
        nvec = np.array([[h, k, l],
                         [-h, -k, -l],
                         [-h, k, l], 
                         [h, -k, -l],
                         [h, -k, l],
                         [-h, k, -l],
                         [-h, -k, l],
                         [h, k, -l]])
        nvec = (1 / np.sqrt(h**2 + k**2 + l**2)) * nvec[::2]
        
    elif h + k == 0 or h + l == 0 or k + l == 0:  
        
        ind = np.nonzero(hkl)
        a = hkl[ind][0]
        nvec = np.array([[a, 0, 0],
                         [-a, 0, 0],
                         [0, a, 0],
                         [0, -a, 0],
                         [0, 0, a],
                         [0, 0, -a]])
        nvec = (1 / np.sqrt(h**2 + k**2 + l**2)) * nvec[::2]
        
    elif (h == k or h == l or k == l) and (h==0 or k==0 or l==0):
        ind = np.nonzero(hkl)
        a = hkl[ind][0]
        
        nvec = np.array([[a, a, 0],
                          [a, 0, a],
                          [0, a, a],
                          [a, -a, 0],
                          [a, 0, -a],
                          [0, -a, a]])
            
        nvec = (1 / np.sqrt(h**2 + k**2 + l**2)) * nvec
        
        
    elif (h == k or h == l or k == l) and h>0 and k>0 and l>0:
        ind = np.nonzero(hkl)
        b = h if h ==l else k if h == k else l
        a = h if h != b else k
        
        nvec = np.array([[a, b, b],
                         [-a, -b, -b], #same
                         [-a, b, b],
                         [a, -b, -b],
                         [a, -b, b],
                         [-a, b, -b],
                         [a, b, -b], # same?
                         [-a, -b, b]])
        nvec = np.vstack([nvec, # hkl
                          nvec[:,[1, 2, 0]], # klh
                          nvec[:,[2, 0, 1]]]) # lhk
          
        nvec = (1 / np.sqrt(h**2 + k**2 + l**2)) * nvec
        nvec = nvec[::2]
   
    elif h != k != l:  
        
        ind = np.nonzero(hkl)
        b = h if h ==l else k if h == k else l
        a = h if h != b else k
        
        nvec = np.array([[a, b, b],
                         [-a, -b, -b], #same
                         [-a, b, b],
                         [a, -b, -b],
                         [a, -b, b],
                         [-a, b, -b],
                         [a, b, -b], # same?
                         [-a, -b, b]])
        nvec = np.vstack([nvec, # hkl
                          nvec[:,[1, 2, 0]], # klh
                          nvec[:,[2, 0, 1]], # lhk
                          nvec[:,[1, 0, 2]], # khl
                          nvec[:,[2, 1, 0]], # lkh
                          nvec[:,[1, 0, 2]]]) # hlk
          
        nvec = (1 / np.sqrt(h**2 + k**2 + l**2)) * nvec[:24:2] # second half is all identical
        #nvec = nvec[::2]
        
    return nvec
