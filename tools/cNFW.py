#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 18:05:32 2021

@author: aknebe
"""
# this tiny code is just mimicking the cNFW calculation from AHF
#===============================================================
import matplotlib.pyplot as plt
import numpy as np

def cNFWroot(c, V2_ratio):
    return (0.216*c/(np.log(1+c)-c/(1+c)) - V2_ratio)


M200 = 1687300000000000.0
r200 = 2192.83/1000.0 * 0.819200 # conversion to Mpc/h by 1/1000., and the last factor is aexp!
vmax = 2165.92
G    = 4.3011e-9

V2_max = vmax**2
V2_vir = G*M200/r200

V2_ratio = V2_max/V2_vir

a=2.2
b=100

while (b-a > 1e-3):
    c = (a+b)/2;
    if (cNFWroot(a, V2_ratio)*cNFWroot(c, V2_ratio) > 0):
        a = c
    else:
        b = c
 
c = (a+b)/2.0


    
