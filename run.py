#!/usr/bin/env python
# coding: utf-8
"""
Created Thu Nov 8 20:52:40 EST 2018
@author: ShengMao

This file runs the hyperelastic
"""

from HyperElasticity import main

"""
Parameters for the main function
    dt, tEnd,
    length=10.0, height=5.0, numElementsFilm=1.5, reGen=True,
    ratLame=5.0, ratFilmSub=100.0
"""
main(0.01, 7.0， numElementsFilm=2)


