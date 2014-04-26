;===============================================================================
; read AHF halos from binary file
;===============================================================================
pro ReadAHFbinaryfile, halodata, nhalos, ifile, AHFfile, ncol
close,1

if VERBOSE_READING EQ 1 then begin
   print,'   reading ',nhalos[ifile],' halos in file ',AHFfile, format='(A,i12,A,A)'
endif


openr,1,DATApath+AHFfile
one        = read_binary(1,DATA_DIMS=0,DATA_TYPE=2)  ; int32_t
numHalos   = read_binary(1,DATA_DIMS=0,DATA_TYPE=15)
numColumns = read_binary(1,DATA_DIMS=0,DATA_TYPE=15)

print,'numHalos   = ',numHalos
print,'numColumns = ',numColumns

if nhalos[ifile] NE numHalos then begin
  print,' you are trying to read a file with ',numHalos,' halos but only allowed for ',nhalos[ifile],' in halodata[]'
  STOP
endif
if ncol NE numColumns then begin
  print,' you are trying to read a file with ',numColumns,' columns but only allowed for ',ncol,' in halodata[]'
  STOP
endif


;for ihalo=0l,numHalos-1l do begin
for ihalo=0l,1l do begin
   halodata[ifile,ihalo,0]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=3)  ; ID
   halodata[ifile,ihalo,1]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=3)  ; hostHalo
   halodata[ifile,ihalo,2]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=3)  ; numSubStruct
   halodata[ifile,ihalo,3]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Mvir
   halodata[ifile,ihalo,4]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=3)  ; npart
   halodata[ifile,ihalo,5]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Xc
   halodata[ifile,ihalo,6]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Yc
   halodata[ifile,ihalo,7]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Zc
   halodata[ifile,ihalo,8]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; VXc
   halodata[ifile,ihalo,9]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; VYc
   halodata[ifile,ihalo,10] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; VZc
   halodata[ifile,ihalo,11] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Rvir
   halodata[ifile,ihalo,12] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Rmax
   halodata[ifile,ihalo,13] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; r2
   halodata[ifile,ihalo,14] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; mbp_offset
   halodata[ifile,ihalo,15] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; com_offset
   halodata[ifile,ihalo,16] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Vmax
   halodata[ifile,ihalo,17] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Vesc
   halodata[ifile,ihalo,18] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; sigV
   halodata[ifile,ihalo,19] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; lambda
   halodata[ifile,ihalo,20] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; lambdaE
   halodata[ifile,ihalo,21] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Lx
   halodata[ifile,ihalo,22] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Ly
   halodata[ifile,ihalo,23] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Lz
   halodata[ifile,ihalo,24] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; b
   halodata[ifile,ihalo,25] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; c
   halodata[ifile,ihalo,26] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Eax
   halodata[ifile,ihalo,27] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Eay
   halodata[ifile,ihalo,28] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Eaz
   halodata[ifile,ihalo,29] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Ebx
   halodata[ifile,ihalo,30] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Eby
   halodata[ifile,ihalo,31] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Ebz
   halodata[ifile,ihalo,32] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Ecx
   halodata[ifile,ihalo,33] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Ecy
   halodata[ifile,ihalo,34] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Ecz
   halodata[ifile,ihalo,35] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; ovdens
   halodata[ifile,ihalo,36] = read_binary(1,DATA_DIMS=0,DATA_TYPE=3)  ; nbins
   halodata[ifile,ihalo,37] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; fMhires
   halodata[ifile,ihalo,38] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Ekin
   halodata[ifile,ihalo,39] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Epot
   halodata[ifile,ihalo,40] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; SurfP
   halodata[ifile,ihalo,41] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; Phi0
   halodata[ifile,ihalo,42] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4)  ; cNFW
   if numColumns GT 43 then begin
      halodata[ifile,ihalo,43] = read_binary(1,DATA_DIMS=0,DATA_TYPE=3) ; n_gas
      halodata[ifile,ihalo,44] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; M_gas
      halodata[ifile,ihalo,45] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; lambda_gas
      halodata[ifile,ihalo,46] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; lambdaE_gas
      halodata[ifile,ihalo,47] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Lx_gas
      halodata[ifile,ihalo,48] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ly_gas
      halodata[ifile,ihalo,49] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Lz_gas
      halodata[ifile,ihalo,50] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; b_gas
      halodata[ifile,ihalo,51] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; c_gas
      halodata[ifile,ihalo,52] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eax_gas
      halodata[ifile,ihalo,53] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eay_gas
      halodata[ifile,ihalo,54] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eaz_gas
      halodata[ifile,ihalo,55] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ebx_gas
      halodata[ifile,ihalo,56] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eby_gas
      halodata[ifile,ihalo,57] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ebz_gas
      halodata[ifile,ihalo,58] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ecx_gas
      halodata[ifile,ihalo,59] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ecy_gas
      halodata[ifile,ihalo,60] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ecz_gas
      halodata[ifile,ihalo,61] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ekin_gas
      halodata[ifile,ihalo,62] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Epot_gas
      halodata[ifile,ihalo,63] = read_binary(1,DATA_DIMS=0,DATA_TYPE=3) ; n_star
      halodata[ifile,ihalo,64] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; M_star
      halodata[ifile,ihalo,65] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; lambda_star
      halodata[ifile,ihalo,66] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; lambdaE_star
      halodata[ifile,ihalo,67] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Lx_star
      halodata[ifile,ihalo,68] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ly_star
      halodata[ifile,ihalo,69] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Lz_star
      halodata[ifile,ihalo,70] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; b_star
      halodata[ifile,ihalo,71] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; c_star
      halodata[ifile,ihalo,72] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eax_star
      halodata[ifile,ihalo,73] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eay_star
      halodata[ifile,ihalo,74] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eaz_star
      halodata[ifile,ihalo,75] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ebx_star
      halodata[ifile,ihalo,76] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eby_star
      halodata[ifile,ihalo,77] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ebz_star
      halodata[ifile,ihalo,78] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ecx_star
      halodata[ifile,ihalo,79] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ecy_star
      halodata[ifile,ihalo,80] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ecz_star
      halodata[ifile,ihalo,81] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ekin_star
      halodata[ifile,ihalo,82] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Epot_star
      if numColumns GT 83 then begin
         halodata[ifile,ihalo,83] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; mean_z_gas
         halodata[ifile,ihalo,84] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; mean_z_star        
      endif
   endif
endfor

close,1
end

;===============================================================================
; read AHF profiles from binary file
;===============================================================================
pro ReadAHFbinaryprofile, profiledata, nhalos, ifile, AHFfile
;pro ReadAHFbinaryprofile,filename,profiledata
close,1

openr,1,DATApath+AHFfile
one        = read_binary(1,DATA_DIMS=0,DATA_TYPE=2)  ; int32_t
numHalos   = read_binary(1,DATA_DIMS=0,DATA_TYPE=15)
numColumns = read_binary(1,DATA_DIMS=0,DATA_TYPE=15)

print,'numHalos   = ',numHalos
print,'numColumns = ',numColumns

if nhalos[ifile] NE numHalos then begin
  print,' you are trying to read a file with ',numHalos,' halos but only allowed for ',nhalos[ifile],' in halodata[]'
  STOP
endif

for ihalo=0l,numHalos-1l do begin
   nbins = read_binary(1,DATA_DIMS=0,DATA_TYPE=15)
   for ibin=0l,nbins-1l do begin
      profiledata[ifile,ihalo,0,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; r
      profiledata[ifile,ihalo,1,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=3) ; npart
      profiledata[ifile,ihalo,2,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; M_in_r
      profiledata[ifile,ihalo,3,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; ovdens
      profiledata[ifile,ihalo,4,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; dens
      profiledata[ifile,ihalo,5,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; vcirc
      profiledata[ifile,ihalo,6,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; vesc
      profiledata[ifile,ihalo,7,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; sigv
      profiledata[ifile,ihalo,8,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Lx
      profiledata[ifile,ihalo,9,ibin]  = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ly
      profiledata[ifile,ihalo,10,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Lz
      profiledata[ifile,ihalo,11,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; b
      profiledata[ifile,ihalo,12,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; c
      profiledata[ifile,ihalo,13,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eax
      profiledata[ifile,ihalo,14,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eay
      profiledata[ifile,ihalo,15,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eaz
      profiledata[ifile,ihalo,16,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ebx
      profiledata[ifile,ihalo,17,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Eby
      profiledata[ifile,ihalo,18,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ebz
      profiledata[ifile,ihalo,19,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ecx
      profiledata[ifile,ihalo,20,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ecy
      profiledata[ifile,ihalo,21,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ecz
      profiledata[ifile,ihalo,22,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Ekin
      profiledata[ifile,ihalo,23,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Epot
      if numColumns GT 24 then begin
         profiledata[ifile,ihalo,24,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; M_gas
         profiledata[ifile,ihalo,25,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; M_star
         profiledata[ifile,ihalo,26,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; u_gas
         if numColumns GT 27 then begin
            profiledata[ifile,ihalo,27,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Z_gas_sh
            profiledata[ifile,ihalo,28,ibin] = read_binary(1,DATA_DIMS=0,DATA_TYPE=4) ; Z_star_sh
         endif
      endif
   endfor
endfor

close,1
end

;===============================================================================
; read AHF halos from ASCII file
;===============================================================================
pro ReadAHFfile, halodata, nhalos, ifile, AHFfile, ncol
@param.h


if VERBOSE_READING EQ 1 then begin
   print,'   reading ',nhalos[ifile],' halos in file ',AHFfile, format='(A,i12,A,A)'
endif

line   = 'c'
idummy = dblarr(ncol)
ihalo  = 0l

openr,1,DATApath+AHFfile
while not EOF(1) do begin

    readf,1,line

    if strcmp(line,'#',1) NE 1 then begin
        reads,line,idummy

        ; move over to halodata[] (assuming the same order!)
        ntransfer = min([ncol,NPROPERTIES])
        for iprop=0,ntransfer-1 do begin
           halodata[ifile,ihalo,iprop] = idummy[iprop]
        endfor

        ihalo = ihalo + 1l
    endif else begin
        ;print,line  ; for debugging purposes
    endelse
    
endwhile
close,1

if VERBOSE_READING EQ 1 then begin
   if ihalo NE nhalos[ifile] then print,'   ** something odd: ihalo=',ihalo,' nhalos[ifile]=',nhalos[ifile]

   print,'       x-range=',min(halodata[ifile,0l:nhalos[ifile]-1l,HXc]),max(halodata[ifile,0l:nhalos[ifile]-1l,HXc])
   print,'       y-range=',min(halodata[ifile,0l:nhalos[ifile]-1l,HYc]),max(halodata[ifile,0l:nhalos[ifile]-1l,HYc])
   print,'       z-range=',min(halodata[ifile,0l:nhalos[ifile]-1l,HZc]),max(halodata[ifile,0l:nhalos[ifile]-1l,HZc])
   print
endif

end

;===============================================================================
; skim all files to get nmaxhalos
;===============================================================================
pro SkimAHFfile, ifile, AHFfile, nhalos
@param.h


nhalos[ifile] = 0l
line = 'c'
openr,1,DATApath+AHFfile
while not EOF(1) do begin

    readf,1,line

    if strcmp(line,'#',1) NE 1 then begin
        nhalos[ifile] = nhalos[ifile] + 1l
    endif else begin
        ;print,line  ; for debugging purposes
    endelse

endwhile
close,1

if VERBOSE_READING EQ 1 then begin
   print,'   found ',nhalos[ifile],' halos in file ',AHFfile, format='(A,i12,A,A)'
endif

end


;===============================================================================
; read all halo catalogue
;===============================================================================
pro initialize, halodata,nhalos, iplot
@param.h

if VERBOSE_READING EQ 1 then begin
   print,'-------------------------------------------------------------'
   print,' reading data for ',NFILES,' codes'
   print,'-------------------------------------------------------------'


   print
   print,' skimming files to get nhalos[] and nmaxhalos...'
   ;--------------------------------------------------------
endif

nhalos = lonarr(NFILES)

for ifile=0,NFILES-1 do begin
   SkimAHFfile, ifile, halofile[ifile], nhalos
endfor

nmaxhalos = max(nhalos)


if VERBOSE_READING EQ 1 then begin
   print
   print,' allocating double-sized halodata[] array:',NFILES,nmaxhalos,NPROPERTIES
   ;-----------------------------------------------------------------------
endif

halodata = dblarr(NFILES,nmaxhalos,NPROPERTIES)


if VERBOSE_READING EQ 1 then begin
   print
   print,' filling halodata[] with meaning...'
   ;------------------------------------------
endif

for ifile=0,NFILES-1 do begin
   ReadAHFfile, halodata, nhalos, ifile, halofile[ifile], ncolumns[ifile]
endfor

print
print

end

;===============================================================================
; close eps output file
;===============================================================================
pro CloseEPS

device,/close
set_plot,'X'
!p.multi = [0,1,1,0,0]

end

;===============================================================================
; open eps output file (or a window)
;===============================================================================
pro OpenEPS,outfile,xxx,yyy
@param.h

SET_PLOT,'ps'
device,filename=outfile
DEVICE, xsize=xxx, ysize=yyy
device,/color,/encapsul

LoadCT, colTable	; Rainbow+white

end

;===============================================================================
; smooth data array
;===============================================================================
pro smooth,data,NBINS,NSMOOTH

for ismooth=1,NSMOOTH do begin

   tmp = data

   data(1) = (tmp(0)+tmp(1)+tmp(2))/3.
   data(2) = (tmp(1)+tmp(2)+tmp(3))/3.

   for i=3,NBINS-4 do begin
      data(i) = (tmp(i-1)+tmp(i)+tmp(i+1))/3.
   endfor

   data(NBINS-3) = (tmp(NBINS-4)+tmp(NBINS-2)+tmp(NBINS-2))/3.
   data(NBINS-2) = (tmp(NBINS-3)+tmp(NBINS-2)+tmp(NBINS-1))/3.

endfor

end

;===============================================================================
; function for overplotting a mass function (not really used!?)
;===============================================================================
pro massfunc,halodata,nhalos,itime,Hmass,lsty,lcol
@param.h

idx    = where(halodata[itime,0:nhalos[itime]-1,Hmass] GT MCUT AND halodata[itime,0:nhalos[itime]-1,Hcom_off]/halodata[itime,0:nhalos[itime]-1,HRvir] LT com_off)
nidx   = N_ELEMENTS(idx)
isort  = sort(halodata[itime,idx,Hmass])

mass  = halodata[itime,idx[isort],Hmass]
Pmass = REVERSE((findgen(nidx)+1))

oplot,mass,Pmass,psym=10,linestyle=lsty,thick=10,color=lcol

end

;=======================================================================
; distance taking into account periodic boundaries
;=======================================================================
function box_distance,U,V,B
D = U-V

idx=where(D GT B/2.)
if idx[0] NE -1 then D[idx] = D[idx]-B

idx=where(D LT -B/2.)
if idx[0] NE -1 then D[idx] = D[idx]+B

;if D GT  B/2. then D = D-B
;if D LE -B/2. then D = D+B
return, D
end


;=======================================================================
; draw a random sample of N points out of an array
;=======================================================================
function draw_sample,array_full,N,ISEED

sample = [-1]
array  = array_full

for i=0l,N-1l do begin
nlength = N_ELEMENTS(array)
iuse    = fix(randomu(ISEED)*float(nlength-1) + 0.5)
sample  = [sample,array[iuse]]
array   = array[where(array NE array[iuse])]
endfor

return, sample[where(sample GE 0)]
end

;===============================================================================
; for fitting a Gaussian
;===============================================================================
pro Gaussian,x,AAA,Func,pder

x0    = AAA[0]
sigma = AAA[1]

Func  = 1./(sqrt(2*!pi)*sigma) * exp(-(x-x0)^2/(2.*sigma^2))

pder  = [ [Func * (x-x0)/sigma^2],$
          [2.*Func/sigma^2*((x-x0)^2/(2.*sigma)-1.)] ]

end

;===============================================================================
; for fitting a lognormal
;===============================================================================
pro lognormal,lambda,AAA,Func,pder

lambda0 = AAA(0)
sigma0  = AAA(1)

x       = (lambda/lambda0)

Func    = 1/sqrt(2*!pi) / sigma0 / lambda * exp(-alog(x)^2/(2*sigma0^2))

pder  = [ [Func * lambda/sigma0^2/lambda0^2 * alog(x)],$
          [Func * (alog(x)^2/sigma0^3 - 1/sigma0)] ]
end

;===============================================================================
; for fitting a power-law
;===============================================================================
pro powerlaw,x,AAA,Func,pder

a0 = AAA[0]
n  = AAA[1]

Func  = a0*x^n

pder  = [ [Func / a0],$
          [Func * alog(x)] ]

end

;===============================================================================
; read the fitparam files
;===============================================================================
pro ReadFitData, fitdata, fitfile
@param.h

nfiles = n_elements(fitfile)
line   = 'dummyline'
for ifile=0,nfiles-1 do begin

    openr,1,DATApath+fitpath+fitfile[ifile]

    while not EOF(1) do begin

        readf,1,line
        reads,line,ilevel,a0,siga0
        readf,1,line
        reads,line,ilevel,n,sign

        fitdata[ifile,ilevel-1,Hn]  = -n
        fitdata[ifile,ilevel-1,Ha0] = a0

    endwhile

    close,1
endfor

end


