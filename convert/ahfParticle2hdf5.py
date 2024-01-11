#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import h5py
from numba import njit

def unique_ordered(x):
    """
    Same output of np.unique when return_counts=True,
    return_index=True, but 1000 faster since x is assumed sorted.
    It works on 1d arrays, untested for other shapes.
    """
    import numpy as np

    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError('x must be a 1D array-like')

    # compare the elements and find where the elements change. This
    # is better than using np.diff, since the comparison always
    # work, while for some arrays (like recarrays) the difference
    # may be not defined
    flt =  np.flatnonzero(x[1:] != x[:-1])

    # by construction, we to add 1 to flt
    off = np.concatenate([[0], flt+1])
    dim = np.concatenate([np.diff(off), [x.size-off[-1]]])
    unique = x[off]

    return unique, off, dim

# use njit only when the number of structures is really high
@njit
def extractParts(Nhaloes, ID, pID, pT, pFoFID, col1, col2, off, dim, off_new):

    for i in range(Nhaloes):
        pID[off_new[i]:off_new[i]+dim[i]] = col1[off[i]:off[i]+dim[i]]
        pT [off_new[i]:off_new[i]+dim[i]] = col2[off[i]:off[i]+dim[i]]

        # note that pFoFID is always in "blocks", but with the MPI version of ahf
        # is not sorted, I do need only to be in blocks
        pFoFID[off_new[i]:off_new[i]+dim[i]] = ID[i]


def extractSinglePartType(pType, pT, pID, pFoFID):
    ind = pT == pType
    N = np.count_nonzero(ind)
    if N > 0:
        newpT = pT[ind]
        newpID = pID[ind]
        un, off, dim = unique_ordered(pFoFID[ind])
    else:
        newpT = []
        newpID = []
        un = []
        off = []
        dim = []

    out = {'N':N, 'pT': newpT, 'pID': newpID, 'unFoFID':un, 'off':off,
            'dim':dim}

    return out

def writeParticleType(h5f, name, dic):
    g = h5f.create_group(name)
    g.attrs.create('Nparticles', dic['N'])

    if dic['N']>0:
        # g['ParticleID'] = dic['pID']
        g.create_dataset('ParticleID',
                data = dic['pID'],
                compression="gzip", compression_opts=9)

        assert dic['N'] == dic['pID'].size

        gD = g.create_group('descriptor')

        gD.create_dataset('FoFID',
                data = dic['unFoFID'],
                compression="gzip", compression_opts=9)
        gD.create_dataset('dim',
                data = dic['dim'],
                compression="gzip", compression_opts=9)
        gD.create_dataset('off',
                data = dic['off'],
                compression="gzip", compression_opts=9)

import time

import argparse

parser = argparse.ArgumentParser(description='Convert a particle file in hdf5')
parser.add_argument('filename', type=str, help='filename')
parser.add_argument('-v','--verbose', help='verbose output',
        action='store_true')
parser.add_argument('-o', action='store_true',
        help='overwrite output files')
parser.add_argument('--no-cast', action='store_true',
        help='Do not try to cast the IDs to a smaller datatype')

args = parser.parse_args()

fname = args.filename
t=time.time()

# read the first line, with the number of "haloes"
Nhaloes = int(np.loadtxt(fname, dtype=np.uint64, max_rows=1, unpack=True))

if Nhaloes == 0:
    print(f'No data in {fname}, skipping this file')
    exit()

# these two readings work, but loadtxt is faster
col1, col2 = np.loadtxt(fname, dtype=np.uint64, skiprows=1, unpack=True)
# tab = Table.read(fname, format='ascii.no_header',
        # guess=False, delimiter='\s', data_start=1,
        # fast_reader=True)
print('time=', time.time()-t)
# np.testing.assert_array_equal(col1, tab['col1'])
# np.testing.assert_array_equal(col2, tab['col2'])
# print(Table.read(fname, format='ascii',
        # delimiter='\s', data_end=1))


dim = np.zeros(Nhaloes, dtype=np.uint64)
ID = np.zeros(Nhaloes, dtype=np.uint64)

dim[0] = col1[0]
ID[0] = col2[0]

tmp = 0


for i in range(1,Nhaloes):
    tmp += int(dim[i-1])+1
    dim[i] = col1[tmp]
    ID[i] = col2[tmp]

off = np.zeros_like(dim)
off[1:] = np.cumsum(dim[:-1])
# this addition accounts for the fact that there is a header for each fof
off+= np.arange(1,dim.size+1, dtype=np.uint64)

# extract the particles
Nparts = dim.sum()
pID = np.zeros(Nparts, dtype=np.uint64)
pT = np.zeros(Nparts, dtype=np.uint64)
off_new = np.zeros_like(dim)
off_new[1:] = np.cumsum(dim[:-1])

# need this to separate the particles
pFoFID = np.zeros(Nparts, dtype=ID.dtype)


t = time.time()

extractParts(Nhaloes, ID, pID, pT, pFoFID, col1, col2, off, dim, off_new)

# # this is a check, recreate the ASCII table
# with open('particle_table.txt', 'w') as ptable:
    # ptable.write(str(dim.sum()) + '\n')
    # for i in range(pID.size):
        # # jj = np.isin(off_new, i)
        # # if np.count_nonzero(jj) == 1:
            # # ptable.write(str(dim[jj])+' '+str(ID[jj]) +' ' + '\n')
        # if i in off_new:
            # jj = np.searchsorted(off_new, i)
            # # ptable.write('              '+'\n')
            # ptable.write(str(dim[jj])+' '+str(ID[jj]) +' ' + '\n')
        # ptable.write(str(pID[i])+' '+str(pT[i])+' '+str(pFoFID[i])+'\n')


# # check ends here

del off
print(time.time()-t)

if not args.no_cast:
    # try to cast it to a lesser precision, if possible

    maxpT = pT.max()
    if np.can_cast(maxpT, np.uint8):
        pT = pT.astype(np.uint8, casting='same_kind')
    elif np.can_cast(maxpT, np.uint16):
        pT = pT.astype(np.uint16, casting='same_kind')
    elif np.can_cast(maxpT, np.uint32):
        pT = pT.astype(np.uint32, casting='same_kind')

    maxpID = pID.max()
    if np.can_cast(maxpID, np.uint8):
        pID = pID.astype(np.uint8, casting='same_kind')
    elif np.can_cast(maxpID, np.uint16):
        pID = pID.astype(np.uint16, casting='same_kind')
    elif np.can_cast(maxpID, np.uint32):
        pID = pID.astype(np.uint32, casting='same_kind')

gas   = extractSinglePartType(0, pT, pID, pFoFID)
DM    = extractSinglePartType(1, pT, pID, pFoFID)
BH    = extractSinglePartType(3, pT, pID, pFoFID)
stars = extractSinglePartType(4, pT, pID, pFoFID)

outname = args.filename+'.h5'
if args.o:
    writeMode = 'w'
else:
    writeMode = 'w-'

with h5py.File(outname, writeMode) as h5f:
    T = h5f.create_group('total')
    # T['pID']=pID
    # T['pT']=pT
    T.create_dataset('pID',
            data = pID,
            compression="gzip", compression_opts=9)
    T.create_dataset('pT',
            data = pT,
            compression="gzip", compression_opts=9)
    Td = T.create_group('descriptor')
    # Td['off']=off_new
    # Td['dim']=dim
    # Td['FoFID'] = ID

    Td.create_dataset('off',
            data = off_new,
            compression="gzip", compression_opts=9)
    Td.create_dataset('dim',
            data = dim,
            compression="gzip", compression_opts=9)
    Td.create_dataset('FoFID',
            data = ID,
            compression="gzip", compression_opts=9)

    writeParticleType(h5f, 'DM', DM)
    writeParticleType(h5f, 'gas', gas)
    writeParticleType(h5f, 'stars', stars)
    writeParticleType(h5f, 'BH', BH)
