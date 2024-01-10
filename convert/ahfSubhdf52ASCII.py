#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import h5py
from pathlib import PurePath
import argparse

def writeAHF_ASCII(filename, IDhaloes,
        Nsubhaloes_perHalo, subhaloesID, off):
    """
    Write the AHF ASCII subhalo file from the numpy array

    """
    with open(filename, 'w') as sub:
        for i in range(IDhaloes.size):

            sub.write(str(IDhaloes[i])+'  '+ str(Nsubhaloes_perHalo[i])+'\n')
            line = ''
            for j in subhaloesID[off[i]:off[i]+Nsubhaloes_perHalo[i]]:
                line+=str(j)+'   '
            line+='\n'
            sub.write(line)

parser = argparse.ArgumentParser(description='Convert the AHF hdf5,'
        ' substructure file obtained with ahfSub2hdf5 in a AHF ASCII file')

parser.add_argument('filename', type=str,
        help='Filename')

parser.add_argument('-v', help='verbose output',
        action='store_true')

args = parser.parse_args()

if args.v:
    print('writing ', outname)

with h5py.File(args.filename, 'r') as h5f:
    IDhaloes = h5f['IDhaloes'][()]
    Nsubhaloes_perHalo = h5f['Nsub_perHalo'][()]
    subhaloesID = h5f['subhaloesID'][()]
    off = h5f['off'][()]

# we are assuming the file has an extension
filename = PurePath(args.filename).with_suffix('.txt')

writeAHF_ASCII(filename, IDhaloes,
        Nsubhaloes_perHalo, subhaloesID, off)
