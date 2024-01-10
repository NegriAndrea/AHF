#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import h5py
from pathlib import PurePath
from astropy.table import Table
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('halofile', type=str,
        help='File of haloes in hdf5 format')
parser.add_argument('-v',help='verbose output',
        action='store_true')
parser.add_argument('-o',help='overwrite',
        action='store_true')

args = parser.parse_args()


path = PurePath(args.halofile).parent
stem = PurePath(args.halofile).name.removesuffix('halos.h5')
subFilename = path/(stem+'substructure.h5')
particleFilename = path/(stem+'particles.h5')

if args.v:
    print('Reading ', args.halofile)
    print('Reading ', subFilename)
    print('Reading ', particleFilename)

with h5py.File(args.halofile, 'r') as h:
    # IDs of subhaloes that contain another subhalo
    IDFull = h['ID'][()]
    NsubStr = h['/numSubStruct'][()]
    hostHaloesID = h['/hostHalo'][()]
    NpartInHalo = h['npart'][()]
    NpartGASInHalo = h['n_gas'][()]
    NpartStarsInHalo = h['n_star'][()]

IDSub_withSubSub = IDFull[(NsubStr > 0) & (hostHaloesID > -1)]

with h5py.File(subFilename, 'r') as s:

    indSubstruc = np.flatnonzero(np.isin(s['IDhaloes'][()],
        IDSub_withSubSub))

    # get the number of subsubhalos
    subID_tot = s['subhaloesID'][()]
    dim = s['/Nsub_perHalo'][indSubstruc]
    off = s['/off'][indSubstruc]

N = dim.sum()
subIDs = []

for i in range(dim.size):
    subIDs.append(subID_tot[off[i]:off[i]+dim[i]])

subIDs = np.concatenate(subIDs)

with h5py.File(particleFilename, 'r') as p:
    ID_partFile = p['/total/descriptor/FoFID'][()]
    dim = p['/total/descriptor/dim'][()]

    DM_ID_partFile = p['DM/descriptor/FoFID'][()]
    DM_dim = p['/DM/descriptor/dim'][()]

    try:
        BH_ID_partFile = p['BH/descriptor/FoFID'][()]
        BH_dim = p['/BH/descriptor/dim'][()]
        BH_present = True
    except KeyError:
        BH_present = False

indParticles = np.flatnonzero(np.isin(ID_partFile, subIDs))

nPartsDM = np.zeros(subIDs.size, dtype=DM_dim.dtype)
if BH_present:
    nPartsBH = np.zeros(subIDs.size, dtype=BH_dim.dtype)

for i in range(subIDs.size):
    mask2 = np.isin(DM_ID_partFile, subIDs[i])
    if np.count_nonzero(mask2) == 1:
        nPartsDM[i] = DM_dim[mask2]
    elif np.count_nonzero(mask2) > 1:
        raise ValueError

    if BH_present:
        mask3 = np.isin(BH_ID_partFile, subIDs[i])
        if np.count_nonzero(mask3) == 1:
            nPartsBH[i] = BH_dim[mask3]
        elif np.count_nonzero(mask3) > 1:
            raise ValueError


ind = np.isin(IDFull, subIDs)

# id of the parent halo with the right sorting, it can contain multiple values
IDparent = hostHaloesID[ind]

# I'm forced to use a loop to take into account multiple values, maybe I
# can make this more efficient, but up to now I have only ~100 elements
NpartsParent = np.zeros_like(IDparent)
for i in range(IDparent.size):
    NpartsParent[i] = NpartInHalo[np.isin(IDFull, IDparent[i])]


np.testing.assert_equal(np.unique(IDparent), IDSub_withSubSub)

# ind2 = np.flatnonzero(ind)
# for i,j in zip(indParticles, ind2):
    # print(ID_partFile[i], dim[i], NpartInHalo[j], DM_dim[i])
    # assert dim[i] ==  NpartInHalo[j]

# for j in ind:
    # print(IDFull[j], NpartInHalo[j])


# quick check, do we have sub-sub-sub haloes?
if NsubStr[ind].sum() > 0:
    print('WARNING: there',
            np.count_nonzero(NsubStr[ind]>0), 'sub-sub-sub haloes!')
    print(IDFull[ind & (NsubStr>0)])

print('There are ', np.count_nonzero(nPartsDM >= 100),
        ' sub sub haloes with >= 100 DM particles')

diff = NpartInHalo[ind]- nPartsDM
print(IDFull[ind][diff <0])
# if any(diff < 0):
    # raise ValueError

ind_sort = np.argsort(subIDs)
np.testing.assert_equal(IDFull[ind], subIDs[ind_sort])
nPartsDM = nPartsDM[ind_sort]
if BH_present:
    nPartsBH = nPartsBH[ind_sort]
else:
    nPartsBH = np.zeros_like(nPartsDM)

t = Table([IDFull[ind], NpartInHalo[ind], nPartsDM,
    NpartGASInHalo[ind], NpartStarsInHalo[ind], nPartsBH,
    IDparent, NpartsParent], 
    names=['ID','nParts','nDMparts', 'nGasParts', 'nStarsParts',
        'nBHparts', 'IDparent', 'NpartsParent'])

np.testing.assert_equal(t['nParts'], t['nDMparts'] + t['nGasParts']+
        t['nStarsParts'] + t['nBHparts'])
print('Writing ', str(path/'SubSubHaloes.txt'))
t.write(path/'SubSubHaloes.txt', format='ascii.ecsv', overwrite=args.o)
