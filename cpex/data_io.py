# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 17:40:07 2015

@author: casimp
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from six import string_types, binary_type
import h5py
import numpy as np


def cpex_to_hdf5(fname, cpex, overwrite=False):
    """ Saves pyxe data object - specifically the analyzed data and basic
    detector/experimental setup - to a hdf5 file. Stores progression point
    so analysis can be continued where it was left off.

    Args:
        fname (str): File path.
        cpex: pyXe data object
        overwrite (bool): Option to overwrite if same filename is specified.
    """
    data_ids = ['sxx', 'syy', 'sxy', 'exx', 'eyy', 'exy', 'eSDV149',
                'ro1', 'ro2', 'ro3', 'x', 'y', 'z', 't', 'grain']

    write = 'w' if overwrite else 'w-'
    with h5py.File(fname, write) as f:

        for name in data_ids:
            
            try:
                d_path = 'cpex/%s' % name
                d = getattr(cpex, name)
                d = d.encode() if isinstance(d, string_types) else d
                if d is not None:
                    f.create_dataset(d_path, data=d)
            except AttributeError:
                pass


def data_extract(cpex, variable_id):
    """ Takes cpex hdf5 file and variable type and extract/returns data."""
    data_ids = {'dims': ['x', 'y', 'z', 'v', 't'],
                'rot': ['ro1', 'ro2', 'ro3'],
                'raw': ['sxx', 'syy', 'sxy', 'exx', 'eyy', 'exy', 'eSDV149']}

    extract = data_ids[variable_id]
    data = []
    for ext in extract:
        try:
            d = cpex['cpex/{}'.format(ext)]
            d = d[()].decode() if isinstance(d[()], binary_type) else d[()]
            data.append(d)
        except KeyError:
            data.append(None)
    return data


