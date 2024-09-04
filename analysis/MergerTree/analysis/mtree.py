#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Just a collection of simple python routines that might or might not be useful
   ( Note, this is a trimmed down version of ahf/analysis/python/ahf.py )


Created on Thu Nov 29 10:39:16 2018

@author: Alexander Knebe

"""
# import all relevant modules
#=================================================
import numpy                as np
import matplotlib.pyplot    as plt
from   mpl_toolkits.mplot3d import axes3d
import gzip

#==============================================================================
# read a Sussing Merger Tree file: returns mtree[][] to be used like this..
#
# mtree[0][:] = haloid
# mtree[1][:] = progid
#
#==============================================================================
def Read_mtree(mtreefile):

#    print('  Reading merger tree file',mtreefile)
    
    # read the merger tree file
    f = open(mtreefile,'r')
    
    f.readline()                                                                # ignore the first two lines
    f.readline()
    
    buf    = f.readline().split()                                               # get the number of halos for which a tree exists
    nhalos = np.uint64(buf[0])
#    print('      ->',nhalos,' connections will be read...')
    
    haloid = np.uint64(np.zeros(nhalos))
    progid = np.uint64(np.zeros(nhalos))
    k      = 0
    for i in range(nhalos):
        buf     = f.readline().split()
        ihaloid = np.uint64(buf[0])
        nprog   = np.int32(buf[1])

        # are there any progenitors? they will be ignored...
        if (nprog>0):
            iprogid = np.uint64(f.readline())
            for j in range(nprog-1):                                            # we are only interested in the main progenitor
                f.readline()
        else:
            iprogid = np.uint64(0)
    
        # store information
        haloid[k] = ihaloid
        progid[k] = iprogid
        k         = k +1
    
    # the MergerTree file might contain entries for halos without a progenitor: remove them manually now...
    ihaveprog = np.where(progid > 0)[0]
    haloid    = haloid[ihaveprog]
    progid    = progid[ihaveprog]
    
    mtree = [haloid,progid]
    
    return mtree
    
#==============================================================================
# create array of same length as np.size(haloids) filling it with 
#          0: no progenitor found
#     progid: haloid of main progenitor
#
# snapid is a single number that indicates that we want progids for this snapid,
# i.e. halos are traced backwards until their Calc_snapids(haloids) matches snapid
# Note, if the snapid of the progenitor is smaller than our desired snapid, the current haloid is returned! 
# Note, snapid is optional: if it is not given, the next progenitor will be returned irrespective of its snapid
#==============================================================================
def Find_progids(mtree, haloids, snapid=[]):
            
    # number of halos to work with
    nhalos = np.size(haloids)
    
    # storage
    progids = np.zeros(nhalos, dtype='uint64')

    # check for usage of snapid=[]
    if(np.size(snapid)==1):
        if(snapid==0):
            return progids


    # only a single snapid is allowed -> return empty progids, if called wrongly
    if(np.size(snapid)>1):
        print('Find_progids() cannot be called with multiple snapids!')
        return progids
            
    # let the index magic begin...
    isin_mtree0        = np.isin(haloids,mtree[0])            # boolean array telling us whether haloid can be found in mtree[0]
    hpos_isinmtree0    = np.where(isin_mtree0 == True)[0]     # positions of the ones with progenitor in haloids[]  -> to be used with mpos_progids[]

    haloids_isinmtree0 = haloids[hpos_isinmtree0]

    
    pos_mtree0               = Find_arraypositions(mtree[0],haloids_isinmtree0)
    progids[hpos_isinmtree0] = mtree[1][pos_mtree0]

    # we synchronize to a preselected snapid:
    if(np.size(snapid)==1):
        # extract snapids from progids and...
        snapids = Calc_snapids(progids)

        # ...leave the haloids whose snapids<snapid untouched
        ibad = np.where((snapids<snapid) & (snapids>0))[0]   # 'bad' are the ones that should be copied untouched
        if(np.size(ibad)>0):
#            print('          ',np.size(ibad),'haloids are left untouched as their progenitors reach too far back...')
            progids[ibad] = haloids[ibad]
            
        # ...continue to follow those haloids whose snapids>snapid
        ibad = np.where(snapids>snapid)[0]   # 'bad' are the ones that require recursive following in mtree[][]
        if(np.size(ibad)>0):
#            print('          ',np.size(ibad),'haloids require consecutive searching...')
            badhaloids    = progids[ibad]
            badprogids    = Find_progids(mtree, badhaloids, snapid=snapid)
            progids[ibad] = badprogids
        
    return progids

#==============================================================================
# find the last possible progenitor id for all ids given in haloids[]
#==============================================================================
def Find_lastprogids(mtree, haloids):
    
    # number of halos to work with
	nhalos = np.size(haloids)
	
	# the core of this routines assumes to be working with numpy-arrays...
#	if(nhalos == 1):
#		haloids = np.array([haloids])

	print('  Finding the very last progids for',nhalos,' provided haloids')

    # storage
	progids = np.copy(haloids)
    
    # flag indicating whether we need to continue
	cont = True
    
	while (cont==True):
        # let the index magic begin...
		isin_mtree0 = np.isin(haloids,mtree[0])            # boolean array telling us whether haloid can be found in mtree[0]
                
		if(np.size(np.where(isin_mtree0 == True))>0):
			hpos_isinmtree0    = np.where(isin_mtree0 == True)[0]     # positions of the ones with progenitor in haloids[]  -> to be used with mpos_progids[]
			haloids_isinmtree0 = haloids[hpos_isinmtree0]
            
			pos_mtree0               = Find_arraypositions(mtree[0],haloids_isinmtree0)
			progids[hpos_isinmtree0] = mtree[1][pos_mtree0]
    
			# now we have the information:
			#  haloids[isnap] and progids[isnap+1]
	
			haloids = np.copy(progids)
		else:
			cont=False
        
	# the core of this routines assumed to be working with numpy-arrays...
#	if(nhalos == 1):
#		progids = progids[0]
		
	return progids

#==============================================================================
# find the last possible progenitor id for all ids given in haloids[]
# snapid_max is the snapid corresponding to the final snapshot of the simulation (i.e. redshift z=0)
#==============================================================================
def Construct_mtreeMatrix(mtree, haloids, snapid_max):

    # number of halos to work with
    nhalos = np.size(haloids)

    # the core of this routines assumes to be working with numpy-arrays...
    # if(nhalos == 1):
    #     haloids = np.array([haloids])

    print('  Constructiong a full merger tree matrix for',nhalos,' provided haloids')

    # make sure that all haloids start with the same snapid
    snapids = Calc_snapids(haloids)
    if(min(snapids) != max(snapids)):
        print('     -> ERROR: all haloids should be from the same snapid')
        return (0)
    else:
        snapid_halos = min(snapids)

    # storage
    progids = np.copy(haloids)
    matrix  = np.zeros((nhalos,snapid_max),dtype=np.uint64)

    # the first entry in the super-matrix is haloids itself
    isnap           = snapid_max-snapid_halos
    matrix[:,isnap] = haloids

    # flag indicating whether we need to continue
    cont = True
    

    while (cont==True):
        isnap += 1

        # let the index magic begin...
        isin_mtree0 = np.isin(haloids,mtree[0])            # boolean array telling us whether haloid can be found in mtree[0]

        if(np.size(np.where(isin_mtree0 == True))>0):
            hpos_isinmtree0    = np.where(isin_mtree0 == True)[0]     # positions of the ones with progenitor in haloids[]  -> to be used with mpos_progids[]
            haloids_isinmtree0 = haloids[hpos_isinmtree0]

            pos_mtree0               = Find_arraypositions(mtree[0],haloids_isinmtree0)
            progids[hpos_isinmtree0] = mtree[1][pos_mtree0]

            # we now need to get rid of progids[] that will not correspond to isnap
            snapid_matrix                   = snapid_max-isnap
            snapids_progids                 = Calc_snapids(progids)
            pos_badsnapids_progids          = np.where(snapids_progids != snapid_matrix)[0]
#           progids[pos_badsnapids_progids] = haloids[pos_badsnapids_progids]  # this drags along the haloid in the matrixwhenever the snapid does not match the extracted snapid from progid
            progids[pos_badsnapids_progids] = 0                                # this puts a 0 in the matrix whenever the snapid does not match the extracted snapid from progid

            # now we have the information haloids[isnap] and progids[isnap+1] and are going to store it in the super-matrix
            matrix[:,isnap] = progids
            
            # but we need to put the 'bad' progids back as we want to continue tracking them
            progids[pos_badsnapids_progids] = haloids[pos_badsnapids_progids]  # we only need this line, if we set progids[]=0 for the non-matching snapid ones

            haloids = np.copy(progids)
        else:
            cont=False

    # the core of this routines assumed to be working with numpy-arrays...
    # if(nhalos == 1):
    #     progids = progids[0]

    return matrix

#==============================================================================
# create array of same length as np.size(haloids) filling it with 
#          0: no descendant found
#     descid: haloid of main descendant
#==============================================================================
def Find_descids(mtree, haloids):
    
    # number of halos to work with
	nhalos = np.size(haloids)
	
	# the core of this routines assumes to be working with numpy-arrays...
#	if(nhalos == 1):
#		haloids = np.array([haloids])
    
    # storage
	descids = np.zeros(nhalos, dtype='uint64')
    
    # let the index magic begin...
	isin_mtree1        = np.isin(haloids,mtree[1])            # boolean array telling us whether haloid can be found in mtree[0]
	hpos_isinmtree1    = np.where(isin_mtree1 == True)[0]     # positions of the ones with progenitor in haloids[]  -> to be used with mpos_progids[]
	haloids_isinmtree1 = haloids[hpos_isinmtree1]
    
	pos_mtree1               = Find_arraypositions(mtree[1],haloids_isinmtree1)
	descids[hpos_isinmtree1] = mtree[0][pos_mtree1]

	# the core of this routines assumed to be working with numpy-arrays...
#	if(nhalos == 1):
#		descids = descids[0]
		
	return descids

#==============================================================================
# extract the snapid from haloid and convert it into isnap 
# (to be used with halos[isnap])
#==============================================================================
def Calc_isnaps(haloids, snapid_max):
    
    snapids = haloids/1e12
        
    if(np.size(haloids)>1):
        snapids = snapids.astype(int)
    else:
        if(isinstance(haloids,int) == False):
            snapids = int(snapids)
        else:
            snapids = np.asscalar(snapids)
    
    isnaps = snapid_max-snapids
    
    return isnaps

#==============================================================================
# extract the snapid from haloid 
#==============================================================================
def Calc_snapids(haloids):
    
    snapids = haloids/1e12
        
    if(np.size(haloids)>1):
        snapids = snapids.astype(int)
    else:
        if(isinstance(snapids,int) == False):
            snapids = int(snapids)
        else:
            snapids = np.asscalar(snapids)
        
    return snapids

#==============================================================================
# convert a Sussing2013 MergerTree file to binary Python file
#==============================================================================
def mtree2npy(file):
    mtree = Read_mtree(file)
    np.save('mtree.npy',mtree)

    return

#==============================================================================
# Returns the positions of the values provided in array1[] within the array2[]
#==============================================================================
def Find_arraypositions(array2, array1):
    
  if(np.size(array1)==1):
    positions = np.where(array2 == array1)[0]
  else:
    sorted_index2    = np.argsort(array2)
    sorted_array2    = array2[sorted_index2]
    sorted_positions = np.searchsorted(sorted_array2, array1)
    positions        = np.take(sorted_index2, sorted_positions, mode='clip')
        
  # tag the positions of array1[] whose values where *not* found in array2[]
  inotfound            = (array2[positions] != array1)
  if(sum(inotfound)>0):
    positions[inotfound] = -1

  return positions
