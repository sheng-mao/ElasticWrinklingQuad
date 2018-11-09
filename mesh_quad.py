#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created Thu Nov 8 20:52:40 EST 2018
@author: ShengMao
"""

import os
from geometry import create_geometry_quad

def gen_mesh(L, H, Hf, N=4):
    dHf = Hf/N;
    # Read mesh and refine once
    fileDir = "./mesh/"
    create_geometry_quad(fileDir, Hf, H, L, dHf, dHf*10)
    # convert geofile to gmsh and then xml
    os.system('gmsh -2 -optimize -order 1 %sgeometry_quad.geo' % (fileDir))
    os.system('meshio-convert -z %sgeometry_quad.msh %sgeometry_quad.xdmf' % (fileDir, fileDir))


