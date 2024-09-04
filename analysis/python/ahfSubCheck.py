#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 10:37:22 2023

@author: aknebe
"""

import glob
import numpy             as np
import ahf               as ahf
import matplotlib.pyplot as plt

# AHF_halosFile = "/Users/aknebe/Office/DATA/Tests/AHF/SubhaloesGoingNotts/test-v111.snap_C02_400_1023.z0.000.AHF_halos"
# hostHaloID    = np.int64(1023000000000001)
# BoxSize       = 100000 # kpc/h

AHF_halosFile = "/Users/aknebe/Office/DATA/Tests/AHF/The300/test-v111.GadgetX-NewMDCLUSTER_0001.z0.000.AHF_halos"
hostHaloID    = np.int64(128000000000001)
BoxSize       = 1000000 # kpc/h

# read in the AHF halo catalogue
halos, Nhalos = ahf.Read_halos_DMonly(AHF_halosFile)

# get access to the hostHalo
pos      = np.where(halos['haloid']==hostHaloID)[0]
hostHalo = halos[pos]

# find all subhaloes based upon hostHalo pointer
pos                 = np.where(halos['hosthaloid']==hostHaloID)[0]
subHalos_hostHalo   = halos[pos]

# find all subhaloes based upon 'D<Rhost' criterion
DistToHost          = np.sqrt( (halos['Xhalo']-hostHalo['Xhalo'])**2 + (halos['Yhalo']-hostHalo['Yhalo'])**2 + (halos['Zhalo']-hostHalo['Zhalo'])**2  )
pos                 = np.where(DistToHost < hostHalo['Rhalo'])[0]
subHalos_DistHost   = halos[pos]
DistToHost          = DistToHost[pos]
subHalos_DistHost   = np.delete(subHalos_DistHost, 0)  # remove the host halo
DistToHost          = np.delete(DistToHost, 0)

# those D<Rhost' subhalos that do not point to hostHalo
pos                         = np.where(subHalos_DistHost['hosthaloid']==0)[0]
# pos                         = np.where(subHalos_DistHost['hosthaloid']!=hostHaloID)[0]
subHalos_DistHost_NoPointer = subHalos_DistHost[pos]
DistToHost_NoPointer        = DistToHost[pos]

print('Distance of NoPointer subhalos to host halo (in units of Rhost):')
print(DistToHost_NoPointer/hostHalo['Rhalo'])

# those D<Rhost' subhalos that do point to another hostHalo
pos                               = np.where((subHalos_DistHost['hosthaloid']!=hostHaloID) & (subHalos_DistHost['hosthaloid']>0))[0]
subHalos_DistHost_MultiplePointer = subHalos_DistHost[pos]
DistToHost_MultiplePointer        = DistToHost[pos]

# the host halo sphere
theta = np.linspace(0,  np.pi,100)
phi   = np.linspace(0,2*np.pi,100)
x     = hostHalo['Rhalo']*np.outer(np.cos(phi),          np.sin(theta)) + hostHalo['Xhalo']
y     = hostHalo['Rhalo']*np.outer(np.sin(phi),          np.sin(theta)) + hostHalo['Yhalo']
z     = hostHalo['Rhalo']*np.outer(np.ones(np.size(phi)),np.cos(theta)) + hostHalo['Zhalo']



# plot host halo sphere and the troublesome subhalos
fig = plt.figure(figsize=(15,15))
ax = plt.axes(projection='3d')

ax.view_init(azim=0)

ax.plot_surface(x, y, z, alpha=0.2)
ax.scatter(subHalos_DistHost_NoPointer['Xhalo'],
           subHalos_DistHost_NoPointer['Yhalo'],
           subHalos_DistHost_NoPointer['Zhalo'],
           s=15*(subHalos_DistHost_NoPointer['Rhalo']))#/np.max(np.log10(subHalos_DistHost_NoPointer['Mhalo'])))

print('Number of subHalos_hostHalo: ',np.size(subHalos_hostHalo))
print('Number of subHalos_DistHost: ',np.size(subHalos_DistHost))
print('Number of subHalos_DistHost w/ embedded hosthaloid: ',np.size(subHalos_DistHost_MultiplePointer))
print('Number of DistHost subhalos with hosthaloid=0:      ',np.size(DistToHost_NoPointer))
