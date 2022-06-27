#!/usr/bin/env python
#  -*- coding: utf-8 -*-

""" Protein-Ligand Binding functions
"""

from numpy import sqrt

"""binding_fraction
Returns the binding fraction of protein
kd is the constant of dissociation
pt is the total concentration of protein
lt is the total concentration of ligand
"""
def binding_fraction(kd, pt, lt):
    a = pt
    b = (kd + lt + pt)
    c = lt
    return((b-sqrt(b*b-4*a*c))/(2*a))

"""koff
Returns the koff for given Kd and kon
"""
def koff(kd, kon):
    return(kd*kon)

