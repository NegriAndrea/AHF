;===============================================================================
pro plot_vmaxfunc
close,1
close,11
close,12
@param.h

;===============================================================================
irestore    = 0
restorefile = 'AHF.DAT'
;===============================================================================

;===============================================================================
; plot specific parameter
;===============================================================================
outfile = "vmaxfunc.eps"
xlabel  = '!8V!6!dmax!n [km/s]'
ylabel  = '!8N!6(>!8V!6!dmax!n)'

;===============================================================================
; open eps file
;===============================================================================
OpenEPS, outfile, 30, 24

;===============================================================================
; read in all data files
;===============================================================================
if irestore EQ 1 then begin
    print,'---------------------------------------------------------------------------------------------------------------'
    print,' restoring data from ',DATApath+restorefile
    print,'---------------------------------------------------------------------------------------------------------------'
    restore,DATApath+restorefile
endif else begin
    initialize, halodata,nhalos, iplot
    ;save,halodata,nhalos,filename=restorefile
endelse

;===============================================================================
; plotting
;===============================================================================
xmin = +1e40
xmax = -1e40
for ifile=0,NFILES-1 do begin
  xmin = 0.75*min([min(halodata[ifile,0l:nhalos[ifile]-1l,HVmax]), xmin])
  xmax = 1.25*max([max(halodata[ifile,0l:nhalos[ifile]-1l,HVmax]), xmax])
endfor
ymin = 0.5
ymax = 1.5*max(nhalos[*])
if SUBHALOS EQ 1 then xmax = 100

print,' Vmax range = ',xmin,xmax
print

!x.range=[xmin,xmax]
!y.range=[ymin,ymax]

print,'--------------------------------------------------------------------'
print,' plotting data'
print,'--------------------------------------------------------------------'

; subhalo mass function ala Springel et al. (2008, Aquarius)
;------------------------------------------------------------
plot_OO,[xmin],[ymin],xstyle=1,ystyle=1,xtitle=xlabel,ytitle=ylabel,xrange=[xmin,xmax],yrange=[ymin,ymax],charsize=csize,charthick=cthick

; numerical functions
;---------------------
CodePlot,halodata,nhalos


;===============================================================================
; open eps file
;===============================================================================
CloseEPS

end

; ===========================================================================

pro CodePlot,halodata,nhalos
@param.h

for ifile=0,NFILES-1 do begin

  if SUBHALOS EQ 1 then begin
    ; set host position
    ;-------------------
    Xhost = halodata[ifile,IHOST,HXc]
    Yhost = halodata[ifile,IHOST,HYc]
    Zhost = halodata[ifile,IHOST,HZc]
    Mhost = halodata[ifile,IHOST,HMvir]

    ; find subhaloes
    ;----------------
    dx   = box_distance(halodata[ifile,0l:nhalos[ifile]-1l,HXc],Xhost,BoxSize)
    dy   = box_distance(halodata[ifile,0l:nhalos[ifile]-1l,HYc],Yhost,BoxSize)
    dz   = box_distance(halodata[ifile,0l:nhalos[ifile]-1l,HZc],Zhost,BoxSize)
    dist = sqrt(dx^2+dy^2+dz^2)
  endif else begin
    dist  = replicate(0.,nhalos[ifile])
    Mhost = 1e20
  endelse
   
  uhalo   = where(dist LE RLIMIT                                          AND $
                  halodata[ifile,0l:nhalos[ifile]-1l,Hnpart] GT NPLIMIT   AND $
                  halodata[ifile,0l:nhalos[ifile]-1l,HMvir]  LT Mhost          , nobj)
  print,' ifile = ',ifile,',  number of haloes for ',halofile[ifile],' =',nobj,nhalos[ifile]

  ; only plot if there are actually any subhaloes
  ;-----------------------------------------------
  if nobj GT 1 then begin

    ; the cumulative distribution
    ;-----------------------------
    Vmax  = halodata[ifile,uhalo,HVmax]
    Vmax  = Vmax[sort(Vmax)]
    PVmax = REVERSE((findgen(nobj)+1))

    oplot,Vmax,PVmax,psym=10,linestyle=sym[ifile]-1,color=col[ifile],thick=5        
  endif

endfor ; ifile

EndPlot:
end




