# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 10:55:48 2019

@author: cs17809
"""

from cpex.scrape import ScrapeODB
import numpy as np

if __name__ == '__main__':
    data = ScrapeODB('/newhome/mi19356/chris_odb/chris_odb.odb')
    print(data.e.shape)
    np.savetxt('test_strain.txt', data.e[0], delimiter=',')
    data.save_cpex('test_cpex.npz')