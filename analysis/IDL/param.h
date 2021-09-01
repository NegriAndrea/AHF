
;***********************************************************************
; choice of units:         Msun/h, kpc/h, km/sec
;***********************************************************************

;=======================================================================
; limitations to be applied when plotting
;=======================================================================

RLIMIT          = 250   ; kpc/h
NPLIMIT         = 20
VMAXLIMIT       = 0.1

;=======================================================================
; the AHF files
;=======================================================================
SubhaloesGoingNotts = 1
LGR2Mpc             = 2
HaloesGoingMAD      = 3
B20a                = 4
Box160              = 5
Box20b              = 6
d2100               = 7

;-----------------------;
; pick data for analyse ;
;-----------------------;
iDATA = SubhaloesGoingNotts


; Subhalos Going Notts
;----------------------
if iDATA EQ SubhaloesGoingNotts then begin
  SUBHALOS = 1
  IHOST    = 0
  BoxSize  = 100000  ; kpc/h
  Omega0   = 0.25
  hubble   = 0.73
  DATApath = "/Users/aknebe/Office/DATA/AHFtests/SubhaloesGoingNotts/"
  halofile = ["test-v003.snap_C02_400_1023.z0.000.AHF_halos",$
              "test-v004.snap_C02_400_1023.z0.000.AHF_halos"]
  ncolumns = [43, 43]
endif

; LGR2Mpc
;---------
if iDATA EQ LGR2Mpc then begin
  SUBHALOS = 1
  IHOST    = 0
  BoxSize  = 64000  ; kpc/h
  Omega0   = 0.24
  hubble   = 0.73
  DATApath = "/Users/aknebe/Office/DATA/AHFtests/LGR2Mpc/"
  halofile = ["test-v003.LGR2Mpc_sfrh_497.z0.000.AHF_halos",$
              "test-v004.LGR2Mpc_sfrh_497.z0.000.AHF_halos"]
  ncolumns = [85, 85]
endif

; Haloes Going MAD
;------------------
if iDATA EQ HaloesGoingMAD then begin
  SUBHALOS = 0
  BoxSize  = 500000  ; kpc/h
  Omega0   = 0.3
  hubble   = 0.7
  DATApath = "/Users/aknebe/Office/DATA/AHFtests/HaloesGoingMAD/"
  halofile = ["test-v003-MPI.FullBox_z0.0_0256.AHF_halos",$
              "test-v004-MPI.FullBox_z0.0_0256.AHF_halos"]
  ncolumns = [43, 43]
endif

; B20a
;------
if iDATA EQ B20a then begin
  SUBHALOS = 0
  BoxSize  = 20000  ; kpc/h
  Omega0   = 0.3
  hubble   = 0.7
  DATApath = "/Users/aknebe/Office/DATA/AHFtests/B20/"
  halofile = ["test-v003.B20a.z0.000.AHF_halos",$
              "test-v004.B20a.z0.000.AHF_halos"]
  ncolumns = [43, 43]
endif

; Box160
;--------
if iDATA EQ Box160 then begin
  SUBHALOS = 0
  BoxSize  = 160000  ; kpc/h
  Omega0   = 0.27
  hubble   = 0.7
  DATApath = "/Users/aknebe/Office/DATA/AHFtests/Box160/"
  halofile = ["test-v003.snap_067.z0.768.AHF_halos",$
              "test-v004.snap_067.z0.768.AHF_halos"]
  ncolumns = [83, 83]
endif

; Box20b
;--------
if iDATA EQ Box20b then begin
  SUBHALOS = 1
  IHOST    = 0
  BoxSize  = 20000  ; kpc/h
  Omega0   = 0.3
  hubble   = 0.7
  DATApath = "/Users/aknebe/Office/DATA/AHFtests/Box20b/"
  halofile = ["test-v003.Box20b.z0.000.AHF_halos",$
              "test-v004.Box20b.z0.000.AHF_halos"]
  ncolumns = [43, 43]
endif

; d2100
;-------
if iDATA EQ d2100 then begin
  SUBHALOS = 0
  BoxSize  = 55520  ; kpc/h
  Omega0   = 0.295
  hubble   = 0.7
  DATApath = "/Users/aknebe/Office/DATA/AHFtests/d2100/"
  halofile = ["test-v003.d2100.z0.000.AHF_halos",$
              "test-v004.d2100.AHF_halos"]
  ncolumns = [82, 82]
endif

NFILES = n_elements(halofile)

;=======================================================================
; misc parameters
;=======================================================================
VERBOSE_READING = 0

colTable      = 39   ; Rainbow+white
c             = { black  :      0,  $
                  blue   :      .25*!d.n_colors, $
                  ltblue :      .40*!d.n_colors, $
                  green  :      .65*!d.n_colors, $
                  yellow :      .75*!d.n_colors, $
                  orange :      .80*!d.n_colors, $
                  red    :      .90*!d.n_colors, $
                  white  :      !d.n_colors-1   }
col           = [c.black,c.ltblue,c.red]
sym           = [1, 3, 5]

;=======================================================================
; the columns within halodata[] array
;=======================================================================
NPROPERTIES = 43   ; this is the number of properties to be used nad stored

Hid      = 0
Hihost   = 1
Hnsub    = 2
HMvir    = 3
Hnpart   = 4
HXc      = 5
HYc      = 6
HZc      = 7
HVXc     = 8
HVYc     = 9
HVZc     = 10
HRvir    = 11
HRmax    = 12
Hr2      = 13
Hmbp_off = 14
Hcom_off = 15
HVmax    = 16
Hvesc    = 17
HsigV    = 18
Hlambda  = 19
HlambdaE = 20
HLx      = 21
HLy      = 22
HLz      = 23
Hb       = 24
Hc       = 25
HEax     = 26
HEay     = 27
HEaz     = 28
HEbx     = 29
HEby     = 30
HEbz     = 31
HEcx     = 32
HEcy     = 33
HEcz     = 34
Hovdens  = 35
Hnbins   = 36
HfMhires = 37
HEkin    = 38
HEpot    = 39
HSurfP   = 40
HPhi0    = 41
HcNFW    = 42
;
; here we can place additional quantities
; so that they are at the same position as
; for the DM_ONLY run
;
nvarEXTRA    = 0
HnvarExtra1  = HcNFW+1
HnvarExtra2  = HcNFW+2
; etc.
;
;
Hngas          = HcNFW+nvarEXTRA+1
HMgas          = HcNFW+nvarEXTRA+2
Hlambda_gas    = HcNFW+nvarEXTRA+3
HlambdaE_gas   = HcNFW+nvarEXTRA+4
HLx_gas        = HcNFW+nvarEXTRA+5
HLy_gas        = HcNFW+nvarEXTRA+6
HLz_gas        = HcNFW+nvarEXTRA+7
Hb_gas         = HcNFW+nvarEXTRA+8
Hc_gas         = HcNFW+nvarEXTRA+9
HEax_gas       = HcNFW+nvarEXTRA+10
HEay_gas       = HcNFW+nvarEXTRA+11
HEaz_gas       = HcNFW+nvarEXTRA+12
HEbx_gas       = HcNFW+nvarEXTRA+13
HEby_gas       = HcNFW+nvarEXTRA+14
HEbz_gas       = HcNFW+nvarEXTRA+15
HEcx_gas       = HcNFW+nvarEXTRA+16
HEcy_gas       = HcNFW+nvarEXTRA+17
HEcz_gas       = HcNFW+nvarEXTRA+18
HEkin_gas      = HcNFW+nvarEXTRA+19
HEpot_gas      = HcNFW+nvarEXTRA+20
Hnstar         = HcNFW+nvarEXTRA+21
HMstar         = HcNFW+nvarEXTRA+22
Hlambda_star   = HcNFW+nvarEXTRA+23
HlambdaE_star  = HcNFW+nvarEXTRA+24
HLx_star       = HcNFW+nvarEXTRA+25
HLy_star       = HcNFW+nvarEXTRA+26
HLz_star       = HcNFW+nvarEXTRA+27
Hb_star        = HcNFW+nvarEXTRA+28
Hc_star        = HcNFW+nvarEXTRA+29
HEax_star      = HcNFW+nvarEXTRA+30
HEay_star      = HcNFW+nvarEXTRA+31
HEaz_star      = HcNFW+nvarEXTRA+32
HEbx_star      = HcNFW+nvarEXTRA+33
HEby_star      = HcNFW+nvarEXTRA+34
HEbz_star      = HcNFW+nvarEXTRA+35
HEcx_star      = HcNFW+nvarEXTRA+36
HEcy_star      = HcNFW+nvarEXTRA+37
HEcz_star      = HcNFW+nvarEXTRA+38
HEkin_star     = HcNFW+nvarEXTRA+39
HEpot_star     = HcNFW+nvarEXTRA+40
Hmean_z_gas    = HcNFW+nvarEXTRA+41
Hmean_z_star   = HcNFW+nvarEXTRA+42


