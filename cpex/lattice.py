# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:48:44 2019

@author: cs17809
"""
from cpex.data_io import cpex_to_hdf5


class Lattice():
    def extract_lattice(self):
        pass
    
    def calc_lattice_tensor(self):
        pass
    
    def save_to_hdf5(self, fpath, overwrite=False):
        """ Save data back to hdf5 file format.

        Saves analyzed information and the raw, scraped data.

        Args:
            fpath (str): Abs. path for new file
            overwrite (bool): Overwrite file if it already exists
        """
        cpex_to_hdf5(fpath, self, overwrite)    