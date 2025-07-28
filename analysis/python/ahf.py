#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 10:39:16 2018

@author: Alexander Knebe
"""
# import all relevant modules
#=================================================
import MyLibrary            as MyL
import numpy                as np
import matplotlib.pyplot    as plt
from   mpl_toolkits.mplot3d import axes3d
import gzip

#==============================================================================
#   halo dtype
#==============================================================================
def def_halostruct():
    Hdesc = [
        ('haloid'           ,  np.uint64),
        ('hosthaloid'       ,  np.uint64),
        ('numSubStruct'     ,  np.int64),
        ('Mhalo'            ,  np.float32),
        ('Npart'            ,  np.int64),
        ('Xhalo'            ,  np.float32),
        ('Yhalo'            ,  np.float32),
        ('Zhalo'            ,  np.float32),
        ('VXhalo'           ,  np.float32),
        ('VYhalo'           ,  np.float32),
        ('VZhalo'           ,  np.float32),
        ('Rhalo'            ,  np.float32),
        ('Rmax'             ,  np.float32),
        ('r2'               ,  np.float32),
        ('mbp_offset'       ,  np.float32),
        ('com_offset'       ,  np.float32),
        ('Vmax'             ,  np.float32),
        ('v_esc'            ,  np.float32),
        ('sigV'             ,  np.float32),
        ('lambda'           ,  np.float32),
        ('lambdaE'          ,  np.float32),
        ('Lx'               ,  np.float32),
        ('Ly'               ,  np.float32),
        ('Lz'               ,  np.float32),
        ('b'                ,  np.float32),
        ('c'                ,  np.float32),
        ('Eax'              ,  np.float32),
        ('Eay'              ,  np.float32),
        ('Eaz'              ,  np.float32),
        ('Ebx'              ,  np.float32),
        ('Eby'              ,  np.float32),
        ('Ebz'              ,  np.float32),
        ('Ecx'              ,  np.float32),
        ('Ecy'              ,  np.float32),
        ('Ecz'              ,  np.float32),
        ('ovdens'           ,  np.float32),
        ('nbins'            ,  np.int32),
        ('fMhires'          ,  np.float32),
        ('Ekin'             ,  np.float32),
        ('Epot'             ,  np.float32),
        ('SurfP'            ,  np.float32),
        ('Phi0'             ,  np.float32),
        ('cNFW'             ,  np.float32),
        ('n_gas'            ,  np.int32),
        ('M_gas'            ,  np.float32),
        ('lambda_gas'       ,  np.float32),
        ('lambdaE_gas'      ,  np.float32),
        ('Lx_gas'           ,  np.float32),
        ('Ly_gas'           ,  np.float32),
        ('Lz_gas'           ,  np.float32),
        ('b_gas'            ,  np.float32),
        ('c_gas'            ,  np.float32),
        ('Eax_gas'          ,  np.float32),
        ('Eay_gas'          ,  np.float32),
        ('Eaz_gas'          ,  np.float32),
        ('Ebx_gas'          ,  np.float32),
        ('Eby_gas'          ,  np.float32),
        ('Ebz_gas'          ,  np.float32),
        ('Ecx_gas'          ,  np.float32),
        ('Ecy_gas'          ,  np.float32),
        ('Ecz_gas'          ,  np.float32),
        ('Ekin_gas'         ,  np.float32),
        ('Epot_gas'         ,  np.float32),
        ('n_stars'          ,  np.int32),
        ('M_stars'          ,  np.float32),
        ('lambda_stars'     ,  np.float32),
        ('lambdaE_stars'    ,  np.float32),
        ('Lx_stars'         ,  np.float32),
        ('Ly_stars'         ,  np.float32),
        ('Lz_stars'         ,  np.float32),
        ('b_stars'          ,  np.float32),
        ('c_stars'          ,  np.float32),
        ('Eax_stars'        ,  np.float32),
        ('Eay_stars'        ,  np.float32),
        ('Eaz_stars'        ,  np.float32),
        ('Ebx_stars'        ,  np.float32),
        ('Eby_stars'        ,  np.float32),
        ('Ebz_stars'        ,  np.float32),
        ('Ecx_stars'        ,  np.float32),
        ('Ecy_stars'        ,  np.float32),
        ('Ecz_stars'        ,  np.float32),
        ('Ekin_stars'       ,  np.float32),
        ('Epot_stars'       ,  np.float32),
        ('mean_z_gas'       ,  np.float32),
        ('mean_z_stars'     ,  np.float32),
        ('n_stars_excised'  ,  np.int32),
        ('M_stars_excised'  ,  np.float32),
        ('mean_z_stars_excised',np.float32)
        ]
    names   = [Hdesc[i][0] for i in range(len(Hdesc))]
    formats = [Hdesc[i][1] for i in range(len(Hdesc))]
    H = np.dtype({'names':names, 'formats':formats}, align=True)
    
    return H

def def_halostruct_DMonly():
    Hdesc = [
        ('haloid'           ,  np.uint64),
        ('hosthaloid'       ,  np.uint64),
        ('numSubStruct'     ,  np.int64),
        ('Mhalo'            ,  np.float32),
        ('Npart'            ,  np.int64),
        ('Xhalo'            ,  np.float32),
        ('Yhalo'            ,  np.float32),
        ('Zhalo'            ,  np.float32),
        ('VXhalo'           ,  np.float32),
        ('VYhalo'           ,  np.float32),
        ('VZhalo'           ,  np.float32),
        ('Rhalo'            ,  np.float32),
        ('Rmax'             ,  np.float32),
        ('r2'               ,  np.float32),
        ('mbp_offset'       ,  np.float32),
        ('com_offset'       ,  np.float32),
        ('Vmax'             ,  np.float32),
        ('v_esc'            ,  np.float32),
        ('sigV'             ,  np.float32),
        ('lambda'           ,  np.float32),
        ('lambdaE'          ,  np.float32),
        ('Lx'               ,  np.float32),
        ('Ly'               ,  np.float32),
        ('Lz'               ,  np.float32),
        ('b'                ,  np.float32),
        ('c'                ,  np.float32),
        ('Eax'              ,  np.float32),
        ('Eay'              ,  np.float32),
        ('Eaz'              ,  np.float32),
        ('Ebx'              ,  np.float32),
        ('Eby'              ,  np.float32),
        ('Ebz'              ,  np.float32),
        ('Ecx'              ,  np.float32),
        ('Ecy'              ,  np.float32),
        ('Ecz'              ,  np.float32),
        ('ovdens'           ,  np.float32),
        ('nbins'            ,  np.int32),
        ('fMhires'          ,  np.float32),
        ('Ekin'             ,  np.float32),
        ('Epot'             ,  np.float32),
        ('SurfP'            ,  np.float32),
        ('Phi0'             ,  np.float32),
        ('cNFW'             ,  np.float32)
        ]
    names   = [Hdesc[i][0] for i in range(len(Hdesc))]
    formats = [Hdesc[i][1] for i in range(len(Hdesc))]
    H = np.dtype({'names':names, 'formats':formats}, align=True)
    
    return H

#==============================================================================
#   halo dtype
#==============================================================================
def def_diskstruct():
    Hdesc = [
        ('haloid'         ,  np.uint64),
        ('r'                ,  np.float32),
        ('Mtot'             ,  np.float64),
        ('M_gas'            ,  np.float32),
        ('Ekin_gas'         ,  np.float32),
        ('k_gas'            ,  np.flot32),
        ('Lx_gas'           ,  np.float32),
        ('Ly_gas'           ,  np.float32),
        ('Lz_gas'           ,  np.float32),
        ('b_gas'            ,  np.float32),
        ('c_gas'            ,  np.float32),
        ('Eax_gas'          ,  np.float32),
        ('Eay_gas'          ,  np.float32),
        ('Eaz_gas'          ,  np.float32),
        ('Ebx_gas'          ,  np.float32),
        ('Eby_gas'          ,  np.float32),
        ('Ebz_gas'          ,  np.float32),
        ('Ecx_gas'          ,  np.float32),
        ('Ecy_gas'          ,  np.float32),
        ('Ecz_gas'          ,  np.float32),
        ('M_stars'            ,  np.float32),
        ('Ekin_stars'         ,  np.float32),
        ('k_stars'            ,  np.flot32),
        ('Lx_stars'           ,  np.float32),
        ('Ly_stars'           ,  np.float32),
        ('Lz_stars'           ,  np.float32),
        ('b_stars'            ,  np.float32),
        ('c_stars'            ,  np.float32),
        ('Eax_stars'          ,  np.float32),
        ('Eay_stars'          ,  np.float32),
        ('Eaz_stars'          ,  np.float32),
        ('Ebx_stars'          ,  np.float32),
        ('Eby_stars'          ,  np.float32),
        ('Ebz_stars'          ,  np.float32),
        ('Ecx_stars'          ,  np.float32),
        ('Ecy_stars'          ,  np.float32),
        ('Ecz_stars'          ,  np.float32)
        ]
    names   = [Hdesc[i][0] for i in range(len(Hdesc))]
    formats = [Hdesc[i][1] for i in range(len(Hdesc))]
    H = np.dtype({'names':names, 'formats':formats}, align=True)
    
    return H

#==============================================================================
# create a structure that holds the snapid-zred mapping
#==============================================================================
def def_snapidmapping():
    Cdesc = [
        ('snapid'            ,  np.uint64),
        ('zred'              ,  np.float32),
        ('aexp'              ,  np.float32),
        ('cosmictime'        ,  np.float32)
        ]
    names   = [Cdesc[i][0] for i in range(len(Cdesc))]
    formats = [Cdesc[i][1] for i in range(len(Cdesc))]
    C = np.dtype({'names':names, 'formats':formats}, align=True)
    
    return C

#==============================================================================
#   reading the file that contains the snapid mapping to z, a, and t
#==============================================================================
def Read_snapidmap(file):
    
    snapidmapping = def_snapidmapping()

    fin = open(file, 'r')                        # open file
    fin.readline()                               # ignore header line
    map = np.loadtxt(file, dtype=snapidmapping)  # read the full halo matrix    
    
    return map
    
#==============================================================================
# read AHF_halos file:
#   either ASCII or .npy, also possible saving the data as .npy file
#==============================================================================
def Read_AHF(halosfile, read_ascii_halos, save_as_npy, fields=[]):

    # read single _halos file into halos()
    if(read_ascii_halos==1):
        halos,nhalos = Read_halos(halosfile,fields)
        if(save_as_npy==1):
            np.save(halosfile+'.npy', halos)
    else:
        halos = np.load(halosfile)
		
    return halos

# read the DMonly version
def Read_AHF_DMonly(halosfile, read_ascii_halos, save_as_npy, fields=[]):

    # read single _halos file into halos()
    if(read_ascii_halos==1):
        halos,nhalos = Read_halos_DMonly(halosfile,fields)
        if(save_as_npy==1):
            np.save(halosfile+'.npy', halos)
    else:
        halos = np.load(halosfile)
        halos = halos[fields]
		
    return halos

#==============================================================================
#   reading _halos file using def_halostruct() as dtype
#
#       careful: it does not work for empty files!
#==============================================================================
def Read_halos(file, fields=[]):
    
    halostruct = def_halostruct()
    if len(fields)==0:
        fields=list(halostruct.names)

    fin = open(file, 'r')                     # open file
    fin.readline()                            # ignore header line
    H   = np.loadtxt(file, dtype=halostruct)  # read the full halo matrix    
    H   = H[fields]                           # restrict to desired fields
    
    return H, np.size(H)
    

# read the DMonly version
def Read_halos_DMonly(file, fields=[]):
    
    halostruct = def_halostruct_DMonly()
    if len(fields)==0:
        fields=list(halostruct.names)

    fin = open(file, 'r')                     # open file
    fin.readline()                            # ignore header line
    H   = np.loadtxt(file, dtype=halostruct)  # read the full halo matrix    
    H   = H[fields]                           # restrict to desired fields
    
    return H, np.size(H)
    
#==============================================================================
#   reading _disks file using def_diskstruct() and def_halostruct() as dtype
#
#       careful: it does not work for empty files!
#==============================================================================
def Read_disks(file):
	
	halostruct = def_halostruct()
	diskstruct = def_diskstruct()

	halos = np.array(1,dtype=halostruct)
	disks = np.array(1,dtype=diskstruct)

	f = open(file, 'r')                     # open file
	f.readline()                            # ignore header line

	# read the halos line
	iline  = 0
	buf    = f.readline().split()
	haloid = np.uint64(buf[0])   # haloid(1)
	nbins  = np.int32(buf[36])   # nbins(37)
	np.append(halos, buf)
	for i in range(nbins):
		buf   = f.readline().split()

	f.close()
	
	return disks
    
#==============================================================================
# read a _particles file
#
# the returndata[][][] shall be used like this:
#      returndata['haloid'][ihalo]    -> haloid of the ihalo-th entry
#      returndata['npart'][ihalo]     -> number of particles in the ihalo-th entry
#      returndata['pids'][ihalo][:]   -> all the (stellar) particle ids in ihalo-th entry
#      returndata['ptype'][ihalo][:]   -> all the (stellar) particle ids in ihalo-th entry
#==============================================================================
def Read_particles(file):
    
    # read the file
    if(file.endswith('.gz')):
        f = gzip.open(file,'r')
    else:
        f = open(file,'r')
        
    buf    = f.readline().split()
    nhalos = np.uint64(buf[0])
    # print('      ->',nhalos,' halos will be read...')
    
    npart    = []
    haloid   = []
    partid   = []
    parttype = []
    
    for i in range(nhalos):
        buf     = f.readline().split()
        inpart  = np.int32(buf[0])
        ihaloid = np.uint64(buf[1])
        
        npart.append(inpart)
        haloid.append(ihaloid)
                        
        ipartid   = []
        iparttype = []
        for j in range(inpart):
            buf = f.readline().split()
            ipartid.append(np.uint64(buf[0]))
            iparttype.append(np.uint32(buf[1]))
        
        partid.append(np.array(ipartid))
        parttype.append(np.array(iparttype))
    
    f.close()

    returndata = {'haloid': np.array(haloid),
                  'npart':  np.array(npart),
                  'pids':   partid,
                  'ptype':  parttype }
        
    return returndata
    
#==============================================================================
# read a _particlesSTARDUST file
#
# the returndata[][][] shall be used like this:
#      returndata['haloid'][ihalo]    -> haloid of the ihalo-th entry
#      returndata['npart'][ihalo]     -> number of particles in the ihalo-th entry
#      returndata['pids'][ihalo][:]   -> all the (stellar) particle ids in ihalo-th entry
#==============================================================================
def Read_particlesSTARDUST(file):
    
    # read the ile
    if(file.endswith('.gz')):
        f = gzip.open(file,'r')
    else:
        f = open(file,'r')
        
    buf    = f.readline().split()
    nhalos = np.uint64(buf[0])
    # print('      ->',nhalos,' halos will be read...')
    
    npart    = []
    haloid   = []
    partid   = []
    partmass = []
    partage  = []
    partmetal= []
    
    for i in range(nhalos):
        buf     = f.readline().split()
        inpart  = np.int32(buf[0])
        ihaloid = np.uint64(buf[1])
        
        npart.append(inpart)
        haloid.append(ihaloid)
                        
        ipartid   = []
        ipartmass = []
        ipartage  = []
        ipartmetal= []
        for j in range(inpart):
            buf = f.readline().split()
            ipartid.append(np.uint64(buf[0]))
            ipartmass.append(np.double(buf[1]))
            ipartage.append(np.double(buf[2]))
            ipartmetal.append(np.double(buf[3]))
        
        partid.append(np.array(ipartid))
        partmass.append(np.array(ipartmass))
        partage.append(np.array(ipartage))
        partmetal.append(np.array(ipartmetal))
    
    f.close()

    returndata = {'haloid': np.array(haloid),
                  'npart':  np.array(npart),
                  'pids':   partid,
                  'mass':   partmass,
                  'age':    partage,
                  'metal':  partmetal}
        
    return returndata
    
#==============================================================================
# read a _starparticlesBCG file (as generated by Ana Contreras)
#
# the returndata[][][] shall be used like this:
#      returndata['BCGsize'][iBCG]   -> size of the BCG
#      returndata['npart'][iBCG]     -> number of particles in the BCG
#      returndata['pids'][iBCG][:]   -> all the (stellar) particle ids
#==============================================================================
def Read_starparticlesBCG(file):
    
    # read the file
    f = open(file,'r')
        
    buf    = f.readline().split()
    nBCGs  = np.uint64(buf[0])
    
    BCGsize = []
    npart   = []
    partid  = []
    masses  = []
    ages    = []
    metals  = []
    
    for i in range(nBCGs):
        buf      = f.readline().split()
        inpart   = np.int32(buf[0])
        iBCGsize = np.uint64(buf[1])
        
        npart.append(inpart)
        BCGsize.append(iBCGsize)
                        
        ipartid = []
        imass   = []
        iage    = []
        imetal  = []
        for j in range(inpart):
            buf = f.readline().split()
            ipartid.append(np.uint64(buf[0]))
            imass.append(np.double(buf[1]))
            iage.append(np.double(buf[2]))
            imetal.append(np.double(buf[3]))
        
        partid.append(np.array(ipartid))
        masses.append(np.array(imass))
        ages.append(np.array(iage))
        metals.append(np.array(imetal))
    f.close()

    returndata = {'BCGsize': np.array(BCGsize),
                  'npart':   np.array(npart),
                  'pids':    partid,
                  'masses':  masses,
                  'ages':    ages,
                  'metals':  metals}
        
    return returndata

#==============================================================================
# read a MergerTree AHF_croco file...
#
# halos[isnap][ihalo]['Prop']
# mtree[]
#
#==============================================================================
def Read_croco(crocofile):
    
    # read the file
    if(crocofile.endswith('.gz')):
        f = gzip.open(crocofile,'r')
    else:
        f = open(crocofile,'r')
        

    # store (at least) the number of particles
    halostruct = def_halostruct_DMonly()
    halos      = []

    # store the main branch info (later to be converted to mtree[][])
    haloid = []
    progid = []

    # ignore the header lines
    buf    = f.readline().split()
    buf    = f.readline().split()
    buf    = f.readline().split()


    # we need to read the file until we hit EOF (= len(buf)==0)
    while True:
        
        # get next data from file (haloid, npart, nprog)
        buf = f.readline()
        if len(buf) == 0:
            break
        
        # store data
        buf = buf.split()
        haloid_i = np.uint64(buf[0])
        npart_i  = np.uint64(buf[1])
        nprog_i  = np.uint32(buf[2])
            
        haloid.append(haloid_i)
        
        halos_i = np.zeros(1,dtype=halostruct)
        halos_i['Npart']  = npart_i
        halos_i['haloid'] = haloid_i
        halos.append(halos_i)
    
        for j in range(nprog_i):
            buf     = f.readline().split()
            nshared_i_j      = np.uint32(buf[0])
            progid_i_j       = np.uint64(buf[1])
            npart_prog_i_j   = np.uint64(buf[2])
            merit_i_j        = np.float64(buf[3])
            
            # we are only interested in the main branch
            if(j==0):
                progid.append(progid_i_j)
            
    f.close()

    # construct the usual mtree[] matrix    
    mtree = []
    for i in range(len(haloid)):
        mtree.append([haloid[i],progid[i]])
        
    return halos, mtree
    

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
    
    pos_mtree0               = MyL.Find_arraypositions(mtree[0],haloids_isinmtree0)
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
            
			pos_mtree0               = MyL.Find_arraypositions(mtree[0],haloids_isinmtree0)
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

            pos_mtree0               = MyL.Find_arraypositions(mtree[0],haloids_isinmtree0)
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
    
	pos_mtree1               = MyL.Find_arraypositions(mtree[1],haloids_isinmtree1)
	descids[hpos_isinmtree1] = mtree[0][pos_mtree1]

	# the core of this routines assumed to be working with numpy-arrays...
#	if(nhalos == 1):
#		descids = descids[0]
		
	return descids

#==============================================================================
# find all the backsplash galaxies
#
#  logic for naming convention:
#
#     hpos_ : index of position in halos[] array
#     mpos_ : index of position in mtree[] array
#    _id(s) : id of halo(s)
#
# Note: assumes that we start at z=0 going backwards
#==============================================================================
def Find_BacksplashGalaxies(halos, mtree, Rcut0=3.0, hosthaloid=128000000000001, show_plot=0):

    print('  Finding backsplash galaxies starting with objects in a region of radius',Rcut0,'x Rhost about host halo with id=',hosthaloid)
    
    # extract host halo informtion at z=0
    hpos_hosthalo = np.where(halos[0]['haloid'] == hosthaloid)[0]            # gives the position of the host halo in the array
    Xhost0        = (halos[0]['Xhalo'])[hpos_hosthalo]
    Yhost0        = (halos[0]['Yhalo'])[hpos_hosthalo]
    Zhost0        = (halos[0]['Zhalo'])[hpos_hosthalo]
    Rhost0        = (halos[0]['Rhalo'])[hpos_hosthalo]

    # extract all possible backsplash candidates from halos[0] (avoiding the hosthalo and restricting to reasonable objects, too)
    Dist0           = MyL.norm2(halos[0]['Xhalo'],halos[0]['Yhalo'],halos[0]['Zhalo'],Xhost0,Yhost0,Zhost0)/Rhost0
    hpos_backsplash = np.where((Dist0 < Rcut0) & (halos[0]['fMhires'] > 0.95) & (halos[0]['haloid'] != hosthaloid))[0]

    Dist0           = Dist0[hpos_backsplash]
    backsplash0     = halos[0][hpos_backsplash]  
    backsplash_ids0 = backsplash0['haloid']

    # keep backsplash_ids0 un-modified
    backsplash_ids = np.copy(backsplash_ids0)
            
    # calculate distance to host for all backsplash candidates
    D = Calc_Dist2HostTree(halos, mtree, backsplash_ids, hosthaloid=hosthaloid)

    # remove zero distances for minimum calculation
    iZero    = np.where(D<1e-16)
    D[iZero] = 1e40

    # find minimum value along nsnaps axis of D[,]
    DistMin = np.min(D,axis=1)

    # eventually identify the ids of the actual backsplash objects
    bpos_backsplash  = np.where((Dist0>1) & (DistMin<1))[0]
    backsplash_ids   = backsplash_ids0[bpos_backsplash]
    backsplash_d0    = Dist0[bpos_backsplash]
    backsplash_dmin  = DistMin[bpos_backsplash]
    backsplash_d     = D[bpos_backsplash]

    # non-found progenitors had distance set to 1e40: reset to 0.0 now...
    pos_reset               = np.where(backsplash_d > 1e39)
    backsplash_d[pos_reset] = 0.0
    
    if(show_plot):
        plt.figure()
        plt.scatter(Dist0,DistMin,s=0.5)
        plt.plot([1,1],[0,1],'y--')
        plt.plot([1,Rcut0],[1,1],'y--')
        plt.xlabel('$D(z=0)/R_{host}$')
        plt.ylabel('$D_{min}/R_{host}$')
        plt.savefig('DminDnow.pdf')
    
    print(' ')
    
    backsplashdata = {'backsplash_ids'  : backsplash_ids,
                      'backsplash_d0'   : backsplash_d0,
                      'backsplash_dmin' : backsplash_dmin,
                      'backsplash_d'    : backsplash_d,
                      'd0'              : Dist0,
                      'dmin'            : DistMin}
    
    return backsplashdata

#==============================================================================
# find all progenitor ids up to end point in tree of a single(!) haloid
#==============================================================================
def Find_allprogids(mtree, haloid):
    
    # number of snapshots read in and number of halos to be traced backwards
    snapid_min = Calc_snapids(mtree[0][-1])
    snapid_max = Calc_snapids(mtree[0][0])
    nsnaps     = snapid_max-snapid_min

    print('  Tracing ',haloid,' using',nsnaps,'snapshots in total')

    snapid  = []
    progids = [haloid]
    isnap   = snapid_max-Calc_snapids(haloid)
    while isnap < nsnaps:
        
        # extract host halo
        print('      tracing halo',isnap,haloid)

        snapid.append(Calc_snapids(haloid))
            
        # move to next halo
        haloid_next = Find_progids(mtree, np.array([haloid]))
        haloid_next = haloid_next[0]
        
        if(haloid_next>0):
            snapid_now  = Calc_snapids(haloid)
            snapid_next = Calc_snapids(haloid_next)
            dsnapid     = snapid_now-snapid_next
            progids.append(haloid_next)
        else:
            print('      -> cannot trace halo anymore')
            break
            
        # next snapshot
        isnap  = isnap + dsnapid
        haloid = haloid_next
                      
    
    print(' ')    
    return progids

#==============================================================================
# trace (single!) haloid
#==============================================================================
def Trace_haloid(halos, mtree, haloid):
    
    # number of snapshots read in and number of halos to be traced backwards
    nsnaps = np.size(halos)
    snapid_max = Calc_snapids(halos[0]['haloid'][0])

    print('  Tracing ',haloid,' using',nsnaps,'snapshots in total')

    X       = []
    Y       = []
    Z       = []
    MAH     = []
    snapid  = []
    isnap   = snapid_max-Calc_snapids(haloid)
    while isnap < nsnaps:
        
        # extract host halo
        pos_hosthalo = np.where(halos[isnap]['haloid'] == haloid)[0] # gives the position of the host halo in the array
        Xhalo        = (halos[isnap]['Xhalo'])[pos_hosthalo]
        Yhalo        = (halos[isnap]['Yhalo'])[pos_hosthalo]
        Zhalo        = (halos[isnap]['Zhalo'])[pos_hosthalo]
        Rhalo        = (halos[isnap]['Rhalo'])[pos_hosthalo]
        Mhalo        = (halos[isnap]['Mhalo'])[pos_hosthalo]
        print('      tracing halo',isnap,haloid,Xhalo,Yhalo,Zhalo,Rhalo,Mhalo)

        snapid.append(Calc_snapids(haloid))
        MAH.append(Mhalo)
        X.append(Xhalo)
        Y.append(Yhalo)
        Z.append(Zhalo)
            
        # move to next halo
        haloid_next = Find_progids(mtree, np.array([haloid]))
        haloid_next = haloid_next[0]
        
        if(haloid_next>0):
            snapid_now  = Calc_snapids(haloid)
            snapid_next = Calc_snapids(haloid_next)
            dsnapid     = snapid_now-snapid_next
        else:
            print('      -> cannot trace halo anymore')
            break
            
        # next snapshot
        isnap  = isnap + dsnapid
        haloid = haloid_next
                      
    returndata = {'MAH'    : MAH,
                  'snapid' : snapid,
                  'X'      : X,
                  'Y'      : Y,
                  'Z'      : Z}
    
    print(' ')    
    return returndata


#==============================================================================
# calculate properties at infall of halos with provided haloids
# Note: assumes that we start at z=0 going backwards
# Note: it also provides some information about the hosthalo...
#==============================================================================
def Find_InfallData(halos, mtree, haloids, hosthaloid=[128000000000001]):
                
    # number of snapshots read in and number of halos to be traced backwards
    nhalos     = np.size(haloids)
    nsnaps     = len(halos)
    snapid_max = Calc_snapids(halos[0]['haloid'][0])

    print('  Finding infall snapids for',nhalos,' halos using',nsnaps,'snapshots in total (assuming snapid_max=',snapid_max,')')

    # create memory for holding information about subhalo status
    #      is_subhalo[ihalo,isnap] =  haloid, if haloid is a subhalo  at isnap (haloid is the id of the progenitor at isnap)
    #      is_subhalo[ihalo,isnap] = -haloid, if it is not a subhalo  at isnap
    #      is_subhalo[ihalo,isnap] = -1,      if it had no progenitor at isnap
    is_subhalo  = np.ones((nhalos,nsnaps), dtype='int64')
    is_subhalo *= -1
    
    # create memory for array holding the infall snapshot and halo ids
    infallhaloids = np.zeros(nhalos, dtype='uint64')
    infallisnaps  = np.ones(nhalos, dtype='int32')
    infallisnaps  = infallisnaps*(nsnaps-1) # set infall time to the maximum possible value, post-correcting it below...

    # create memory for array holding all distances of all halos at all times (-1 for isnaps without progenitor)
    D  = np.ones((nhalos,nsnaps), dtype='float')
    D *= -1

    # create memory for keeping track of mass of host halo
    Mhosthalo = []
    Vmaxhosthalo = []
    isnaphalo = []
    Xhosthalo = []
    Yhosthalo = []
    Zhosthalo = []
    
    # trace hosthalo backwards
    while (hosthaloid > 0):
        isnap      = Calc_isnaps(hosthaloid,snapid_max)
        
        pos_hosthalo = np.where(halos[isnap]['haloid'] == hosthaloid)[0] # gives the position of the host halo in the array
        Xhost        = (halos[isnap]['Xhalo'])[pos_hosthalo]
        Yhost        = (halos[isnap]['Yhalo'])[pos_hosthalo]
        Zhost        = (halos[isnap]['Zhalo'])[pos_hosthalo]
        Rhost        = (halos[isnap]['Rhalo'])[pos_hosthalo]
        Mhost        = (halos[isnap]['Mhalo'])[pos_hosthalo]
        Vmaxhost     = (halos[isnap]['Vmax'])[pos_hosthalo]
        Mhosthalo.append(Mhost)
        Vmaxhosthalo.append(Vmaxhost)
        isnaphalo.append([isnap])
        Xhosthalo.append(Xhost)
        Yhosthalo.append(Yhost)
        Zhosthalo.append(Zhost)
        print('      tracing host halo',hosthaloid,Xhost,Yhost,Zhost,Rhost,Mhost,'in isnap',isnap)

        # mask out bad haloids (i.e. the ones for which we could not find a progenitor at this isnap)
        pos_good_haloids = np.where(Calc_snapids(haloids) == snapid_max-isnap)[0]

        # check subhalo status for the good ones
        if(len(pos_good_haloids) > 0):
            good_haloids = haloids[pos_good_haloids]
            
            # positions in halos[i] where our haloids can be found
            pos_good_halos    = MyL.Find_arraypositions(halos[isnap]['haloid'],good_haloids)                
            
            # store distance, too
            X                         = halos[isnap]['Xhalo'][pos_good_halos]
            Y                         = halos[isnap]['Yhalo'][pos_good_halos]
            Z                         = halos[isnap]['Zhalo'][pos_good_halos]
            D[pos_good_haloids,isnap] = MyL.norm2(X,Y,Z, Xhost,Yhost,Zhost)/Rhost

            # classify haloids as subhalo or not subhalo
#            pos_subhalos    = np.where(halos[isnap]['hosthaloid'][pos_good_halos] == hosthaloid)[0]    # AHF criterion
#            pos_notsubhalos = np.where(halos[isnap]['hosthaloid'][pos_good_halos] != hosthaloid)[0]    # AHF criterion
            pos_subhalos    = np.where((D[pos_good_haloids,isnap] <= 1) & (D[pos_good_haloids,isnap] > 0))[0]    # distance criterion
            pos_notsubhalos = np.where( D[pos_good_haloids,isnap] >  1)[0]                                       # distance criterion
            
#            subhaloids = good_haloids[pos_subhalos] 
            
            is_subhalo[pos_good_haloids[pos_subhalos],   isnap] =  good_haloids[pos_subhalos] 
            is_subhalo[pos_good_haloids[pos_notsubhalos],isnap] = -good_haloids[pos_notsubhalos]

        # move to next snapshot (allowing for skips in the MergerTree!)
        hosthaloid = Find_progids(mtree, hosthaloid)
        haloids    = Find_progids(mtree, haloids, snapid=Calc_snapids(hosthaloid)) 
    
    # use the matrix is_subhalo[][] to identify infall times (now we loop over each halo...)
    for ihalo in range(nhalos):
        
        # first infall
        pos_isnap_subhalo = np.where(is_subhalo[ihalo] > 0)[0]
        if(len(pos_isnap_subhalo)>0):
            
            # this is the first time the halo entered its host
            infallisnap = pos_isnap_subhalo.max()
            
            # but we want its properties while it was still outside (if possible)
            if(is_subhalo[ihalo][min(infallisnap+1,nsnaps-1)] < -1):
                infallisnap = min(infallisnap+1,nsnaps-1)
                
            infallisnaps[ihalo]  = infallisnap           
            infallhaloids[ihalo] = abs(is_subhalo[ihalo][infallisnap])
        else:
            infallisnaps[ihalo]  = nsnaps-1  
            infallhaloids[ihalo] = 0
            
           
    # convert to snapid
    infallsnapids = snapid_max-infallisnaps
    
    # prepare dictionary with return data
    infalldata = {'infallsnapids' : infallsnapids,
                  'infallhaloids' : infallhaloids,
                  'distances'     : D,
                  'MAHhost'       : np.array(Mhosthalo),
                  'Vmaxhost'      : np.array(Vmaxhosthalo),
                  'isnap'         : np.array(isnaphalo),
                  'Xhost'         : np.array(Xhosthalo),
                  'Yhost'         : np.array(Yhosthalo),
                  'Zhost'         : np.array(Zhosthalo)
                  }
    
    print(' ')    
    return infalldata  # we also return the distance to the host: we were able to calculate that and hence it might be of use...

#==============================================================================
# calculate the distance to the host for provided haloids
# Note, distances are normalized to the host radius
#==============================================================================
def Calc_Dist2Host(halos, haloids, hosthaloid=128000000000001):
    
    # number of snapshots read in and number of halos to be traced backwards
    nhalos     = np.size(haloids)
    
    # create memory for array holding all distances of all halos at all times
    D = np.zeros((nhalos), dtype='float')        
            
    # extract host halo
    hpos_hosthalo = np.where(halos['haloid'] == hosthaloid)[0] # gives the position of the host halo in the array
    Xhost         = (halos['Xhalo'])[hpos_hosthalo]
    Yhost         = (halos['Yhalo'])[hpos_hosthalo]
    Zhost         = (halos['Zhalo'])[hpos_hosthalo]
    Rhost         = (halos['Rhalo'])[hpos_hosthalo]
    Mhost         = (halos['Mhalo'])[hpos_hosthalo]
    
    # mask out bad haloids (i.e. the ones for which we could not find a progenitor)
    pos_good_haloids = np.where((haloids > 0) & (haloids != hosthaloid))[0]

    if(len(pos_good_haloids) > 0):
        good_haloids = haloids[pos_good_haloids]
        
        # positions in halos[isnap] where our haloids are found
        pos_haloids  = MyL.Find_arraypositions(halos['haloid'],good_haloids)
        X            = halos['Xhalo'][pos_haloids]
        Y            = halos['Yhalo'][pos_haloids]
        Z            = halos['Zhalo'][pos_haloids]
        
        # store distance
        D[pos_good_haloids] = MyL.norm2(X,Y,Z, Xhost,Yhost,Zhost)/Rhost
        
    
#    print(' ')    
    return D

#==============================================================================
# follow all halos as selected by haloids and track them backwards in time
#   -> here we calculate the distance to the host with hosthaloid
# Note: - assumes that we start at z=0 going backwards
#       - distances are normalized to the host radius
#==============================================================================
def Calc_Dist2HostTree(halos, mtree, haloids, hosthaloid=128000000000001):
    
    # number of snapshots read in and number of halos to be traced backwards
    nhalos     = np.size(haloids)
    nsnaps     = len(halos)
    snapid_max = Calc_snapids(halos[0]['haloid'][0])
    
	# the core of this routines assumes to be working with numpy-arrays...
    if(nhalos == 1):
        haloids = np.array([haloids])

    print('  Calculating distance to host',hosthaloid,' for',nhalos,'halos across',nsnaps,'snapshots')

    # create memory for array holding all distances of all halos at all times
    D = np.zeros((nhalos,nsnaps), dtype='float')
          
    # fill in D[,]
    isnap = 0
    while isnap < nsnaps:
        
            
        # extract host halo
        hpos_hosthalo = np.where(halos[isnap]['haloid'] == hosthaloid)[0] # gives the position of the host halo in the array
        Xhost         = (halos[isnap]['Xhalo'])[hpos_hosthalo]
        Yhost         = (halos[isnap]['Yhalo'])[hpos_hosthalo]
        Zhost         = (halos[isnap]['Zhalo'])[hpos_hosthalo]
        Rhost         = (halos[isnap]['Rhalo'])[hpos_hosthalo]
        Mhost         = (halos[isnap]['Mhalo'])[hpos_hosthalo]
        print('      tracing host halo',hosthaloid,Xhost,Yhost,Zhost,Rhost,Mhost)
        
        # mask out bad haloids (i.e. the ones for which we could not find a progenitor)
        pos_good_haloids = np.where((haloids > 0) & (haloids != hosthaloid) & (Calc_isnaps(haloids,snapid_max) == isnap))[0]

        if(len(pos_good_haloids) > 0):
            good_haloids = haloids[pos_good_haloids]
            
            # positions in halos[isnap] where our haloids are found
            pos_haloids  = MyL.Find_arraypositions(halos[isnap]['haloid'],good_haloids)
            X            = halos[isnap]['Xhalo'][pos_haloids]
            Y            = halos[isnap]['Yhalo'][pos_haloids]
            Z            = halos[isnap]['Zhalo'][pos_haloids]
            
            # store distance
            D[pos_good_haloids,isnap] = MyL.norm2(X,Y,Z, Xhost,Yhost,Zhost)/Rhost
        
        # move to next host halo
        hosthaloid_next = Find_progids(mtree, np.array([hosthaloid]))
        hosthaloid_next = hosthaloid_next[0]
                
        if(hosthaloid_next>0):
            hostsnapid      = Calc_snapids(hosthaloid)
            hostsnapid_next = Calc_snapids(hosthaloid_next)
            dsnapid         = hostsnapid-hostsnapid_next
            
            # in case the merger tree was not done using consecutive snapids:
            # dsnapid = 1
            
            # next snapshot
            isnap      = isnap + dsnapid
            hosthaloid = hosthaloid_next
    
            # locate progenitor haloids @ hostsnapid_next
            haloids = Find_progids(mtree, haloids, snapid=hostsnapid_next)
        else:
            print('      -> cannot trace host halo anymore')
            break

    
    print(' ')    
    return D

#==============================================================================
# find the isnap and actual value of the minimum distance to the host halo
#   for a given list of haloids
#==============================================================================
def Find_minD(halos, mtree, haloids, hosthaloid=128000000000001):
    
    # get the matrix that contains the distance at all times, i.e. D[ihalo,isnap]
    D = Calc_Dist2HostTree(halos, mtree, haloids, hosthaloid)
    
    # scan that matrix for the minimum in each line
    pos_bad     = D<=0
    D[pos_bad]  = 1e40  # this puts the non-tracked values out of the min() range
    minD        = np.amin(D,axis=1)

    # find the position of the minimum for each halo
    nhalos     = np.size(haloids)
    minD_isnap = np.zeros(nhalos)
    for i in range(nhalos):
        pos           = MyL.Find_arraypositions(D[i,:], minD[i])
        minD_isnap[i] = pos[0]
        
    returndata = {'minD'       : minD,
                  'minD_isnap' : minD_isnap}
    
    return returndata

#==============================================================================
# extract halos with haloids[] from halos[][][]
#==============================================================================
def Extract_halos(halos, haloids):
    
    # number of halos
    nhalos     = np.size(haloids)  # this should be identical to np.size(snapids)!
    snapid_max = Calc_snapids(halos[0]['haloid'][0])
    
    # create space to extracted halos
    halostruct      = def_halostruct()
    halos_extracted = np.zeros(nhalos,dtype=halostruct)
    
    for ihalo in range(nhalos):
        haloid                 = haloids[ihalo]
        if(haloid>0):
            isnap                  = Calc_isnaps(haloid, snapid_max)
            pos_halos              = np.where(halos[isnap]['haloid'] == haloid)[0]
            if(len(pos_halos > 0)):
                halos_extracted[ihalo] = halos[isnap][pos_halos]
    
    return halos_extracted

#==============================================================================
# plot a Dnow-Dmin plot (typical for backsplash galaxies)
#==============================================================================
def Plot_DminDnow(halos_d0, halos_dmin):
    
    plt.figure()
    plt.scatter(halos_d0,halos_dmin,s=0.1)
    plt.xlabel('$D(z=0)/R_{host}$')
    plt.ylabel('$D_{min}/R_{host}$')
    plt.savefig('DminDnow.pdf')

    return


#==============================================================================
# plot the orbits of all selected halos
#==============================================================================
def Plot_Orbits(halos, mtree, haloids, hosthaloid=128000000000001, endsnapids=[], distances=[], snapidmap=[]):
    
    nhalos     = np.size(haloids)
    nsnaps     = np.size(halos)
    snapid_max = Calc_snapids(halos[0]['haloid'][0])
    
    if(len(endsnapids) == 0):
        endisnaps = np.ones(nhalos,dtype='int32')*nsnaps
    else:
        endisnaps = np.array(snapid_max-endsnapids)
            
    if(len(distances) == 0):
        distances = Calc_Dist2HostTree(halos, mtree, haloids, hosthaloid=hosthaloid)
           
    plt.figure()
    # extract lines out of D[,] plotting the orbits as a function of time
    for ihalo in range(nhalos):
        snapidrange = snapid_max-np.linspace(0,endisnaps[ihalo]-1,endisnaps[ihalo])
        if(len(snapidmap)==0):
            plt.plot(snapidrange, distances[ihalo,0:endisnaps[ihalo]])
        else:
            zreds = Find_zreds(snapidmap,snapidrange)
            plt.plot(zreds, distances[ihalo,0:endisnaps[ihalo]])
            
    if(len(snapidmap)==0): 
        plt.plot([endisnaps.min(),snapid_max],[1,1])
        plt.xlabel('snapid')
    else:
        plt.plot([0,Find_zreds(snapidmap,endisnaps.min())],[1,1])
        plt.xlabel('redshift $z$')
        
    plt.ylabel('$D/R_{host}$')
    plt.savefig('orbits.pdf')

    return distances
    

#==============================================================================
# follow all the halos given by haloids[] extracting full information and
# placing it into a matrix halos_traced[nhalos,nsnaps]
# Note:
#        - assumes that we start at z=0 going backwards
#==============================================================================
def Trace_halos(halos, mtree, haloids_original):
    
    haloids = np.copy(haloids_original)
    
    # number of snapshots read in and number of halos to be traced backwards
#    nsnaps     = np.size(halos)
    nsnaps     = len(halos)
    nhalos     = np.size(haloids)
    
    print('  Tracing',nhalos,' halos across',nsnaps,'snapshots using their progenitor ids to extract full halo information')
    
    # create space for extracted halos
    halostruct   = def_halostruct()
    halos_traced = np.zeros((nhalos,nsnaps),dtype=halostruct)
    
    # fill in starting values
    pos_now           = MyL.Find_arraypositions(halos[0]['haloid'], haloids)
    ifound_now        = np.where(pos_now >= 0)  # 08/11/2023: changed '>' to '>='
    halos_traced[:,0] = halos[0][pos_now[ifound_now]]
    
    # loop over all snapshots extracting progenitors
    for isnap in range(1,nsnaps):
        
        # fill 'next' values
        progids               = Find_progids(mtree, haloids)
        
        # the ones that have been found
        pos                   = MyL.Find_arraypositions(halos[isnap]['haloid'], progids)
        ifound                = np.where(pos >= 0)   # 08/11/2023: changed '>' to '>='
        if (np.size(ifound) > 0):
            halos_traced[ifound,isnap] = halos[isnap][pos[ifound]]
    
            # update current haloids[] to those values for which a progenitor had been found (and used...)
            haloids[ifound] = progids[ifound]

        # the ones that have not been found
        # pos                   = MyL.Find_arraypositions(halos[isnap]['haloid'], progids)
        # inotfound             = np.where(pos < 0)
        # if (np.size(inotfound) > 0):
        #     halos_traced[inotfound,isnap] = -1
        
    print(' ')    
    return halos_traced

#==============================================================================
# follow all the halos given by haloids[] extracting full information and
# placing it into a matrix halos_traced[nhalos,nsnaps]
# Note:
#        - assumes that we start at z=0 going backwards
#==============================================================================
def Trace_halos_DMonly(halos, mtree, haloids_original):
    
    haloids = np.copy(haloids_original)
    
    # number of snapshots read in and number of halos to be traced backwards
#    nsnaps     = np.size(halos)
    nsnaps     = len(halos)
    nhalos     = np.size(haloids)
    
    print('  Tracing',nhalos,' halos across',nsnaps,'snapshots using their progenitor ids to extract full halo information')
    
    # create space for extracted halos
    halostruct   = def_halostruct_DMonly()
    halos_traced = np.zeros((nhalos,nsnaps),dtype=halostruct)
    
    # fill in starting values
    pos_now           = MyL.Find_arraypositions(halos[0]['haloid'], haloids)
    ifound_now        = np.where(pos_now >= 0)  # 08/11/2023: changed '>' to '>='
    halos_traced[:,0] = halos[0][pos_now[ifound_now]]
    
    # loop over all snapshots extracting progenitors
    for isnap in range(1,nsnaps):
        
        # print(isnap)
        # print(len(halos[isnap]))
        
        # fill 'next' values
        progids               = Find_progids(mtree, haloids)
        pos                   = MyL.Find_arraypositions(halos[isnap]['haloid'], progids)
        ifound                = np.where(pos >= 0)   # 08/11/2023: changed '>' to '>='
        if (np.size(ifound) > 0):
            halos_traced[ifound,isnap] = halos[isnap][pos[ifound]]
    
            # update current haloids[] to those values for which a progenitor had been found (and used...)
            haloids[ifound] = progids[ifound]

    print(' ')    
    return halos_traced

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
# return zreds[] for provided snapids[] using snapidmap
#==============================================================================
def Find_zreds(snapidmap, snapids):

    pos_snapids = MyL.Find_arraypositions(snapidmap['snapid'],snapids)
    zreds       = snapidmap['zred'][pos_snapids]
    
    return zreds

#==============================================================================
# return times[] for provided snapids[] using snapidmap
#==============================================================================
def Find_cosmictimes(snapidmap, snapids):

    pos_snapids = MyL.Find_arraypositions(snapidmap['snapid'],snapids)
    cosmictimes = snapidmap['cosmictime'][pos_snapids]
    
    return cosmictimes

#==============================================================================
# return subhaloids for a preselected host,
#  either using AHF subhalo assignment or a radial cut Rfrac
#==============================================================================
def Find_subhaloids(halos, hosthaloid, Rfrac_min=0.0, Rfrac=-1.0, isnap=0):
    	
    if(Rfrac>0):
        pos_hosthalo = np.where(halos[isnap]['haloid']==hosthaloid)[0]
        Xhost = halos[isnap]['Xhalo'][pos_hosthalo]
        Yhost = halos[isnap]['Yhalo'][pos_hosthalo]
        Zhost = halos[isnap]['Zhalo'][pos_hosthalo]
        Rhost = halos[isnap]['Rhalo'][pos_hosthalo]
        
        # subhalos: radial distance selection
        pos_nonhost = np.where(halos[isnap]['haloid']!=hosthaloid)[0]
        dX   = halos[isnap]['Xhalo'][pos_nonhost]-Xhost
        dY   = halos[isnap]['Yhalo'][pos_nonhost]-Yhost
        dZ   = halos[isnap]['Zhalo'][pos_nonhost]-Zhost
        Dist = np.sqrt(dX**2+dY**2+dZ**2)
        
        pos_subhalos = np.where((Rfrac_min*Rhost <= Dist) & (Dist <= Rfrac*Rhost))[0]			
        subhaloids   = halos[isnap]['haloid'][pos_nonhost][pos_subhalos]
    else:
        # subhalos: AHF selection criterion
        subhaloids = halos[isnap]['haloid'][np.where(halos[isnap]['hosthaloid'] == hosthaloid)[0]]
    
    return subhaloids

#==============================================================================
# convert a Sussing2013 MergerTree file to binary Python file
#==============================================================================
def mtree2npy(file):
    mtree = Read_mtree(file)
    np.save('mtree.npy',mtree)

    return

#==============================================================================
# convert ASCII AHF_halos file to binary Python file
#==============================================================================
def halos2npy(file):
    halos = Read_halos(file)
    np.save('halos.npy',halos)
    
    return

#==============================================================================
# visualize the 3D positions of supplied haloids
#==============================================================================
def Plot_Positions3D(halos, haloids):
	
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	ax = fig.gca(projection='3d')
	
	ihalos = MyL.Find_arraypositions(halos['haloid'], haloids)
	
	Xhalo = halos['Xhalo'][ihalos]
	Yhalo = halos['Yhalo'][ihalos]
	Zhalo = halos['Zhalo'][ihalos]
	
	# choose which coordinates to plot
	x = Xhalo
	y = Yhalo
	z = Zhalo
	
	# get maximum range
	xmin = x.min()
	ymin = y.min()
	xmax = x.max()
	ymax = y.max()
	min  = np.min(np.array([xmin, ymin]))
	max  = np.max(np.array([xmax, ymax]))
	
	ax.scatter(x,y,z,'.',s=0.75)
	
#	plt.xlim((min,max))
#	plt.ylim((min,max))
#	plt.zlim((min,max))
	plt.axis('equal')    
    
#	plt.xlabel('x [kpc/h]')
#	plt.ylabel('y [kpc/h]')
#	plt.zlabel('z [kpc/h]')
	
	plt.savefig('positions.pdf')
	
	return

#==============================================================================
# visualize the 2D positions of supplied haloids
#==============================================================================
def Plot_Positions2D(halos, haloids):
	
	plt.figure()
	plt.axis('equal')    

	ihalos = MyL.Find_arraypositions(halos['haloid'], haloids)
	
	Xhalo = halos['Xhalo'][ihalos]
	Yhalo = halos['Yhalo'][ihalos]
	
	# choose which coordinates to plot
	x = Xhalo
	y = Yhalo
	
	# get maximum range
	xmin = x.min()
	ymin = y.min()
	xmax = x.max()
	ymax = y.max()
	min  = np.min(np.array([xmin, ymin]))
	max  = np.max(np.array([xmax, ymax]))
	
	plt.plot(x,y,'.',markersize=0.75)
	
	plt.xlim((min,max))
	plt.ylim((min,max))
    
	plt.xlabel('x [kpc/h]')
	plt.ylabel('y [kpc/h]')
	
	plt.savefig('positions.pdf')
	
	return

#==============================================================================
# remove all halos with haloids
#==============================================================================
def Remove_halos(halos, haloids):
    
    mask  = ~np.isin(halos['haloid'],haloids)
    halos = halos[mask]
                    
    return halos

#==============================================================================
# remove all halos with haloids
#==============================================================================
def Keep_halos(halos, haloids):
    
    mask  = np.isin(halos['haloid'],haloids)
    halos = halos[mask]
                    
    return halos

#==============================================================================
# extract only the unique pids from two _particlesSTARDUST files
#  -> this will give you all the ids of stars found in subhaloes
#==============================================================================
def Extract_uniqueSTARDUSTids(file1, file2):
    
    pSTARDUST1 = Read_particlesSTARDUST(file1)
    pSTARDUST2 = Read_particlesSTARDUST(file2)
    
    nhalos1 = len(pSTARDUST1['haloid'])
    nhalos2 = len(pSTARDUST2['haloid'])
    
    if(nhalos1 != nhalos2):
        print(' halo numbers do not match:',nhalos1,' vs ',nhalos2)
        
    else:
        uniqueids = []
    
        for ihalo in range(nhalos1):
        
            pids       = pSTARDUST1['partid'][ihalo]
            pidsremove = pSTARDUST2['partid'][ihalo]
        
            uniqueids.append(np.array(MyL.Remove_ids(pids,pidsremove)))
    
    returndata = {'nhalos':  nhalos1,
                  'pids':    uniqueids}
    
    return returndata

#==============================================================================
# convert all lengths to physical units
#==============================================================================
def physical_units(AHF, a, h):
    
    AHF['Xhalo']      = AHF['Xhalo']      * a/h
    AHF['Yhalo']      = AHF['Yhalo']      * a/h
    AHF['Zhalo']      = AHF['Zhalo']      * a/h
    AHF['Rhalo']      = AHF['Rhalo']      * a/h
    AHF['Rmax']       = AHF['Rmax']       * a/h
    AHF['r2']         = AHF['r2']         * a/h
    AHF['mbp_offset'] = AHF['mbp_offset'] * a/h
    AHF['com_offset'] = AHF['com_offset'] * a/h

    AHF['Mhalo']           = AHF['Mhalo']      /h
    AHF['M_gas']           = AHF['M_gas']      /h
    AHF['M_stars']         = AHF['M_stars']    /h
    AHF['M_stars_excised'] = AHF['M_stars_excised']    /h

    return AHF


#==============================================================================
# the prev111/112 sub-subhalo bug is annoying :-(
#==============================================================================
def Clean_particlesSTARDUSTexcised(STARDUST, clusterid):
    
    # extract all pids belonging to clusterid
    icluster_pos = np.where(STARDUST['haloid'] == clusterid)[0]
    clusterpids  = np.uint64(STARDUST['pids'][icluster_pos[0]])

    # find all duplicate pids
    N       = len(STARDUST['pids'])
    allpids = np.array([],dtype='uint64')
    for i in range(N):
        if(len(STARDUST['pids'][i])>0):
            allpids = np.append(allpids,(STARDUST['pids'][i]))
    duppids = MyL.Find_Duplicates(allpids)
        
    # if there are duplicates, remove them from 'clusterid'
    if(len(duppids)>0):
        tag_bad  = MyL.Find_arraypositions(duppids, clusterpids)
        pos_good = np.where(tag_bad < 0)[0] # if the pid was not found in duppids there will be a '-1' in tag_bad
    
        # only keep the unique particles
        STARDUST['pids'][icluster_pos[0]]  = STARDUST['pids'][icluster_pos[0]][pos_good]
        STARDUST['mass'][icluster_pos[0]]  = STARDUST['mass'][icluster_pos[0]][pos_good]
        STARDUST['age'][icluster_pos[0]]   = STARDUST['age'][icluster_pos[0]][pos_good]
        STARDUST['metal'][icluster_pos[0]] = STARDUST['metal'][icluster_pos[0]][pos_good]
        STARDIST['npart'][icluster_pos[0]] = len(STARDUST['pids'][icluster_pos[0]])
    
    return STARDUST

#==============================================================================
# remove duplicate pids from clusterid
#==============================================================================
def Excise_particles(PARTICLES, clusterid):
    
    # extract all pids belonging to clusterid
    icluster_pos = np.where(PARTICLES['haloid'] == clusterid)[0]
    clusterpids  = np.uint64(PARTICLES['pids'][icluster_pos[0]])

    # create array allpids[] that contains all pids, but not the one of clusterid
    N       = len(PARTICLES['pids'])
    allpids = np.array([],dtype='uint64')
    for i in range(N):
        if((len(PARTICLES['pids'][i])>0) & (PARTICLES['haloid'][i]!=clusterid)):
            allpids = np.append(allpids,(PARTICLES['pids'][i]))
            
    # remove allpids from the clusterpids, METHOD A
    mask = np.isin(clusterpids, allpids, invert=True)

    PARTICLES['pids'][icluster_pos[0]]  = PARTICLES['pids'][icluster_pos[0]][mask]
    PARTICLES['ptype'][icluster_pos[0]] = PARTICLES['ptype'][icluster_pos[0]][mask]
    PARTICLES['npart'][icluster_pos[0]] = len(PARTICLES['pids'][icluster_pos[0]])
          
    # remove allpids from the clusterpids, METHOD B
    # duppids = MyL.Find_Duplicates(allpids)
        
    # # if there are duplicates, remove them from 'clusterid'
    # if(len(duppids)>0):
    #     tag_bad  = MyL.Find_arraypositions(duppids, clusterpids)
    #     pos_good = np.where(tag_bad < 0)[0] # if the pid was not found in duppids there will be a '-1' in tag_bad
    
    #     # only keep the unique particles
    #     PARTICLES['pids'][icluster_pos[0]]  = PARTICLES['pids'][icluster_pos[0]][pos_good]
    #     PARTICLES['ptype'][icluster_pos[0]] = PARTICLES['ptype'][icluster_pos[0]][pos_good]
    #     PARTICLES['npart'][icluster_pos[0]] = len(PARTICLES['pids'][icluster_pos[0]])
        
    return PARTICLES, icluster_pos[0]
