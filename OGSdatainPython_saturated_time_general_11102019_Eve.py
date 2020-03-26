# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:20:27 2019

@author: khurana
"""

import numpy  as np
import csv
import Pythonfunctions_SK as sk
filename = 'model_domain_quad.tec'
fwithd = filename
print("Reading tec file....")
size, steps, Headers, D = sk.readTecfile(fwithd)
print("Converting to array....")
df = sk.Converttomarr(D)
print("Saving numpy array....")
np.save('numparray_test_df',df)
