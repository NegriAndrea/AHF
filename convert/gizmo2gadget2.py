###############################################################################
# convert GIZMO runs to Gadget2 file
###############################################################################

import h5py
import numpy as np
import sys
import os


def write_head(fp, header, size):
    fp.write(np.int32(8))
    fp.write(header.encode())
    fp.write(np.int32(size+8))
    fp.write(np.int32(8))


def write_header(fl, attr, pattr):
    # Adam wrote two types of particles into the snapshot: non-existing gas[0] and DM[1],
    # but we need to fill space for all 6 species
    fl.write(np.int32(256))
    fl.write(np.asarray([0, attr['NumPart_ThisFile'][1], 0, 0, 0, 0], dtype=np.int32))
    fl.write(np.asarray([0, attr['MassTable'][1], 0, 0, 0, 0], dtype=np.float64))
    fl.write(np.asarray([attr['Time'], attr['Redshift']], dtype=np.float64))
    fl.write(np.array([1, 1], dtype=np.int32))
    fl.write(np.asarray([0, attr['NumPart_Total'][1], 0, 0, 0, 0], dtype=np.int32))
    fl.write(np.array([1, attr['NumFilesPerSnapshot']], dtype=np.int32))
    fl.write(np.asarray([pattr['BoxSize'], pattr['Omega0'], pattr['OmegaLambda'], pattr['HubbleParam']], dtype=np.float64))
    fl.write(np.zeros(np.int32(256/4)-6-12-4-2-6-2-8, dtype=np.int32))
    fl.write(np.int32(256))

###############################################################################
#                               MAIN
###############################################################################

# take folder name from command line
if len(sys.argv) != 2:
    print('I require only one argument -- the folder name!')
    print('I exit, because you providing is not fit!', str(sys.argv))
    exit(0)
else:
    h5files = [sys.argv[1]]

# h5files = ["/Users/aknebe/Office/DATA/Projects/WDM2023/simulations/wdm/snapshot_006.hdf5"]

# loop to get all required information
for filename in h5files:
    if filename[-4:] == 'hdf5' or filename[-4:] == 'HDF5':
        
        # open HDF5 file for reading
        print('opening file',filename)
        f = h5py.File(filename, "r")

        # open GADGET2 file for writing
        G2filename = filename[:-5]
        of = open(G2filename, 'wb')
        print('writing to file',G2filename)
        
        attrs  = f['/Header'].attrs
        pattrs = f['Parameters'].attrs     # BoxSize etc. are found here in Adam's WDM runs
        TNp = attrs['NumPart_ThisFile'][:6]
        TNs = np.sum(TNp)

        # HEAD
        write_head(of, "HEAD", 256)
        write_header(of, attrs, pattrs)

        # POS
        write_head(of, "POS ", np.uint32(TNs*4*3))
        of.write(np.uint32(TNs*4*3))
        if 'PartType0' in f.keys():
            of.write(np.float32(f['/PartType0/Coordinates'][:]))
        if 'PartType1' in f.keys():
            of.write(np.float32(f['/PartType1/Coordinates'][:]))
        if 'PartType2' in f.keys():
            of.write(np.float32(f['/PartType2/Coordinates'][:]))
        if 'PartType3' in f.keys():
            of.write(np.float32(f['/PartType3/Coordinates'][:]))
        if 'PartType4' in f.keys():
            of.write(np.float32(f['/PartType4/Coordinates'][:]))
        if 'PartType5' in f.keys():
            of.write(np.float32(f['/PartType5/Coordinates'][:]))
        of.write(np.uint32(TNs*4*3))

        # VEL
        write_head(of, "VEL ", np.uint32(TNs*4*3))
        of.write(np.uint32(TNs*4*3))
        if 'PartType0' in f.keys():
            of.write(np.float32(f['/PartType0/Velocities'][:]))
        if 'PartType1' in f.keys():
            of.write(np.float32(f['/PartType1/Velocities'][:]))
        if 'PartType2' in f.keys():
            of.write(np.float32(f['/PartType2/Velocities'][:]))
        if 'PartType3' in f.keys():
            of.write(np.float32(f['/PartType3/Velocities'][:]))
        if 'PartType4' in f.keys():
            of.write(np.float32(f['/PartType4/Velocities'][:]))
        if 'PartType5' in f.keys():
            of.write(np.float32(f['/PartType5/Velocities'][:]))
        of.write(np.uint32(TNs*4*3))

        # ID
        write_head(of, "ID  ", np.uint32(TNs*4))
        of.write(np.uint32(TNs*4))
        if 'PartType0' in f.keys():
            of.write(np.uint32(f['/PartType0/ParticleIDs'][:]))
        if 'PartType1' in f.keys():
            of.write(np.uint32(f['/PartType1/ParticleIDs'][:]))
        if 'PartType2' in f.keys():
            of.write(np.uint32(f['/PartType2/ParticleIDs'][:]))
        if 'PartType3' in f.keys():
            of.write(np.uint32(f['/PartType3/ParticleIDs'][:]))
        if 'PartType4' in f.keys():
            of.write(np.uint32(f['/PartType4/ParticleIDs'][:]))
        if 'PartType5' in f.keys():
            of.write(np.uint32(f['/PartType5/ParticleIDs'][:]))
        of.write(np.uint32(TNs*4))

        # MASS  -> there is no mass block needed as all particles have the same mass
        # write_head(of, "MASS", np.uint32(TNs*4))
        # of.write(np.uint32(TNs*4))
        # if 'PartType0' in f.keys():
        #     of.write(np.float32(f['/PartType0/Masses'][:]))
        # if 'PartType1' in f.keys():
        #     of.write(np.float32(f['/PartType1/Masses'][:]))
        # if 'PartType2' in f.keys():
        #     of.write(np.float32(f['/PartType2/Masses'][:]))
        # if 'PartType3' in f.keys():
        #     of.write(np.float32(f['/PartType3/Masses'][:]))
        # if 'PartType4' in f.keys():
        #     of.write(np.float32(f['/PartType4/Masses'][:]))
        # if 'PartType5' in f.keys():
        #     of.write(np.float32(f['/PartType5/Masses'][:]))
        # of.write(np.uint32(TNs*4))
