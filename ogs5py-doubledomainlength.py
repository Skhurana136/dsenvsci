# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:20:04 2020

@author: khurana
"""
#Generating double domain lengths for testing simulations

from ogs5py import OGS

model = OGS(task_root = "double", task_id = "doublemodel")
model.msh.generate("rectangular", dim = 2, mesh_origin = (0.0, 0.0, 0.0), element_no = (30, 100, 0), element_size = (0.01, 0.01, 0))
model.write_input()

#Generating half domain lengths for testing simulations

from ogs5py import OGS

model = OGS(task_root = "half", task_id = "halfmodel")
model.msh.generate("rectangular", dim = 2, mesh_origin = (0.0, 0.0, 0.0), element_no = (30, 25, 0), element_size = (0.01, 0.01, 0))
model.write_input()

#Generating 1.25m domain lengths for testing simulations

from ogs5py import OGS

model = OGS(task_root = "half", task_id = "doublehalfmodel")
model.msh.generate("rectangular", dim = 2, mesh_origin = (0.0, 0.0, 0.0), element_no = (30, 125, 0), element_size = (0.01, 0.01, 0))
model.write_input()