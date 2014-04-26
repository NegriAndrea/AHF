C============================================================================================
C
C
C        OBSOLETE CODE: HAS BEEN REPLACED BY THE MORE SOPHISTICATED ahfHaloHistory.c
C
C
C HaloHistory
C
C purpose:     follow a specified halo throughout *.AHF_halos files
C
C
C input:       - halo id
C              - NTIMES     redshifts
C              - NTIMES     AHF/BDM_halos files
C              - (NTIMES-1) MergerTree_idx files
C
C output:      - one file with the info:
C
C                z_1       Npart NVpart Xc YC Zc etc.
C                z_2       Npart NVpart Xc Yc Zc etc.
C                 ...
C                z_NTIMES  Npart NVpart Xc Yc Zc etc.   (same format as AHF/BDM_halos)
C
C============================================================================================
      program HaloHistory

      parameter(NMAXTIMES = 1000)

      real*4        zred(NMAXTIMES)
      character*512 halofile(NMAXTIMES), idxfile(NMAXTIMES), outfile
      character*512 dummyline
      integer       npart, nbins
      real*8        nvpart,Xc,Yc,Zc,VXc,VYv,VZc,Mvir,Rvir,Vmax
      real*8        Rmax,sigV,lambda,Lx,Ly,Lz,ovdens,Redge
      real*8        a,E1x,E1y,E1z,b,E2x,E2y,E2z,c,E3x,E3y,E3z
      real*8        Ekin,Epot,mbp_offset,com_offset,r2

      print*,'================'
      print*,'  Halo History'
      print*,'================'
 	

C  USER INPUT
C=============
      write(*,'(A,$)')' please give id of halo to trace       => '
      read*,ihalo
      print*,ihalo
      write(*,'(A,$)')' please give number of redshifts       => '
      read*,NTIMES
      print*,NTIMES

      print*
      do itime=1,NTIMES
         write(*,'(A,i3,A,$)')' please give        ',itime,
     &                   '. redshift       => '
         read*,zred(itime)
         print*,zred(itime)
      enddo

      print*
      do itime=1,NTIMES
         write(*,'(A,i3,A,$)')' please give name of',itime,
     &                   '. AHF_halos file => '
         read(*,'(A)')halofile(itime)
         CALL get_name(halofile(itime),id,jd)
         print*,halofile(itime)(id:jd)
      enddo

      print*
      do itime=1,NTIMES-1
         write(*,'(A,i3,A,$)')' please give name of',itime,
     &                   '. MergerTree_idx file => '
         read(*,'(A)')idxfile(itime)
         CALL get_name(idxfile(itime),id,jd)
         print*,idxfile(itime)(id:jd)
      enddo
      
      print*
      write(*,'(A,$)')' please give name for output file => '
      read(*,'(A)')outfile
      CALL get_name(outfile,id,jd)
      print*,outfile(id:jd)


C  START TRACKING
C==================
      open(11,file=outfile)
      iprev = ihalo
      do itime=1,NTIMES

         ! find ihalo in AHF/BDM_halos file
         if(ihalo .NE. iprev)then
            CALL get_name(halofile(itime),id,jd)
            print*,' o tracing halo',ihalo,
     &             ' in file ',halofile(itime)(id:jd)
         endif
         open(12,file=halofile(itime))
         read(12,'(A)')dummyline
         jhalo = -1
         do while(jhalo .NE. ihalo)
c            read(12,*,end=666)npart, nvpart, 
c     &                Xc, Yc, Zc, VXc, VYv, VZc, 
c     &                Mvir, Rvir, Vmax, Rmax, sigV, 
c     &                lambda, Lx, Ly, Lz, 
c     &                a, E1x, E1y, E1z,
c     &                b, E2x, E2y, E2z,
c     &                c, E3x, E3y, E3z,
c     &                ovdens, Redge, nbins,
c     &                Ekin, Epot, mbp_offset, com_offset, r2
            read(12,*,end=666)haloID, ihostHalo, numSubStruct, 
     &                Mvir, npart, Xc, Yc, Zc, VXc, VYv, VZc, 
     &                Rvir, Rmax, r2, mbp_off, com_off, 
     &                Vmax, v_esc, sigV, xlambda, xlambdaE,
     &                xLx,xLy,xLz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,
     &                Ecx,Ecy,Ecz,ovdens,nbins,fMhires,
     &                Ekin,Epot,SurfP,Phi0,cNFW
            jhalo = jhalo + 1
         enddo
         close(12)

         ! write info to outfile
c         write(11,999)zred(itime),npart, nvpart, 
c     &                Xc, Yc, Zc, VXc, VYv, VZc, 
c     &                Mvir, Rvir, Vmax, Rmax, sigV, 
c     &                lambda, Lx, Ly, Lz, 
c     &                a, E1x, E1y, E1z,
c     &                b, E2x, E2y, E2z,
c     &                c, E3x, E3y, E3z,
c     &                ovdens, Redge, nbins,
c     &                Ekin,Epot,mbp_offset,com_offset,r2
          write(11,888)haloID, ihostHalo, numSubStruct, 
     &                Mvir, npart, Xc, Yc, Zc, VXc, VYv, VZc, 
     &                Rvir, Rmax, r2, mbp_off, com_off, 
     &                Vmax, v_esc, sigV, xlambda, xlambdaE,
     &                xLx,xLy,xLz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,
     &                Ecx,Ecy,Ecz,ovdens,nbins,fMhires,
     &                Ekin,Epot,SurfP,Phi0,cNFW
 888     FORMAT(3i12,g14.6,i12,3f16.8,3f10.2,5f16.8,2f10.2,
     &          17f10.6,i10,f10.6,4g16.6,f12.6)
 999     FORMAT(f8.3,i10,f10.0,3f10.4,3f10.2,g14.6,f10.2,f8.2,f10.2,
     &       f8.2,f10.6,3g12.4,12f10.6,f8.2,f10.2,i8,2g16.6,3f14.5)

         ! find ihalo in idxfile
         if (itime .NE. NTIMES)then
            open(12,file=idxfile(itime))
            jhalo = -1
            do while(jhalo .NE. ihalo)
               read(12,*,end=666)jhalo,idx
            enddo
            close(12)
         endif

         iprev = ihalo
         ihalo = idx

      enddo
 666  close(11)


      STOP
      end
*----------------------------------------------------------------------------

*****************************************************************************

*----------------------------------------------------------------------------
      SUBROUTINE get_name(name,i,j)
      CHARACTER name*(*)
      DO i=1,256
         IF(name(i:i).NE.' ') GOTO 1
      ENDDO
      i=1
      j=1
      RETURN
 1    CONTINUE
      DO j=i,256
         IF(name(j:j).EQ.' ') GOTO 2
      ENDDO
      j=121
 2    j=j-1
      IF(name(i:i).EQ.'''') i=i+1
      IF(name(j:j).EQ.'''') j=j-1
      END
*-----------------------------------------------------------------------------

