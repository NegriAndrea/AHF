C==========================================================================================
C HaloHaloXi
C
C calculates halo-halo two-point correlation function based upon *.AHF_halos
C
C==========================================================================================
      program HaloHaloXi

c      include 'omp_lib.h'


*  propmeans-Files
*  ---------------
      character*120 datapath,data,filename
      character*3   char_z
      real z(30),zi,MMIN,MMAX
      integer fcount

      print*,'============================================'
      print*,'  halo-halo two-point correlation function'
      print*,'============================================'

      write(*,'(A,$)')'name of halo catalogues                   => '
      read(*,'(A)')data
      print*,data

      write(*,'(A,$)')'how many points for xi(r)                 => '
      read*,NC
      print*,NC
      write(*,'(A,$)')'mass cut (0) oder number density (1)      => '
      read*,icut
      print*,icut
      if(icut.EQ.0)then
         write(*,'(A,$)')'  Mmin   = '
         read*,MMIN
         print*,MMIN
         write(*,'(A,$)')'  Mmax   = '
         read*,MMAX
         print*,MMAX
      else
         write(*,'(A,$)')'  nhalos = '
         read*,NHALOS
         print*,NHALOS
      end if
      !write(*,'(A,$)')'logarithmic (0) oder linear (1) binning    => '
      !read*,ilogi
      !print*,ilogi
      ilogi = 0

      write(datapath,'(''./'')')
      call get_name(datapath,ip,jp)

      write(*,'(A,$)')'box [same units as in input file]          => '
      read*,rLbox
      print*,rLbox
      print*


      print*
      print*,'--------------------------------------------------'
      print*

      virial = 0
      call CorrelationFunction(datapath,data,rLbox,NC,NHALOS,MMIN,MMAX,
     .                 ilogi,icut)

      print*,'STOP'
      STOP
      end
*=============================================================================

******************************************************************************

*-----------------------------------------------------------------------------
      subroutine CorrelationFunction(datapath,data,rLbox,NC,NHALOS,MMIN,
     .                               MMAX,ilogi,icut)

      character*120 datapath,data,corrfile
      real rLbox,virial,MMIN,MMAX

      parameter (MAX_NCL = 5000000)    ! Maximalzahl aller Cluster
      parameter (RADFAKMIN  = 0.00004)
      parameter (RADFAKMAX  = 0.3)        ! Faktor, der angibt, bis zuwelchem Anteil
                                          ! der Boxlaenge die Korrelationsfunktion 
                                          ! berechnet werden soll
*  Koordianten der Cluster
*  -----------------------
      real r(3,MAX_NCL)

*  Fuer die Korrelationsfunktion
*  -----------------------------
      DIMENSION CR(NC), WP(NC), IC(NC), WC(NC)
      integer NSUB,isub(MAX_NCL),ISEED
      real    lRA,lRE

      ISEED = -1001               ! ISEED muss eine Variable sein !
      NSUB  = 2000000             ! NSUB  muss ebenfalls eine Variable sein !

      call get_name(data,id,jd)

      N=MAX_NCL                   ! Vorabschaetzung

      call read_halos(datapath,data,r,N,NHALOS,MMIN,MMAX,icut)

      print*,'total number of clusters used :',N

      if(N .LE. NSUB)then
         NSUB = N
         do i=1,NSUB
            isub(i) = i
         end do
      else
         call subsample(N,NSUB,ISEED,isub)
      end if

c      open(1,file='subsample.dat')
c      do i=1,NSUB
c         write(1,*) r(1,isub(i)),r(2,isub(i)),r(3,isub(i))
c      end do
c      close(1)

      if(icut .EQ. 0)then
         corrfile=data(id:jd)//'_xi.mass'
      else
         corrfile=data(id:jd)//'_xi.ncl'
      endif

      OPEN(18,FILE=corrfile)

C   Radiusbereich der  Korrelationsfunktion
*   --------------------------------------- 
      RA  = RADFAKMIN*rLbox
      RE  = RADFAKMAX*rLbox
      lRA = alog10(RA)
      lRE = alog10(RE)
      DO I = 1, NC
         if(ilogi.EQ.0)then
            CR(I) = 10.**(lRA + REAL(I-1)*real(lRE-lRA)/real(NC))
         else
            CR(I) = RA + REAL(I-1)*real(RE-RA)/real(NC)
         end if
      END DO

      RNORM = REAL(NSUB)*REAL(NSUB-1)/2.0

C   Analysen
*   --------
      CALL ANA(r,N,isub,RA,RE,NSUB,IC,NC,rLbox,ilogi)

C   Normierung
*   ----------
      CALL NORM(rLbox,RA,RE,WP,NC,ilogi)

      print*,'writing correlation function: ',corrfile
      print*
      print*,'-----------------------------------------------'

      DO I = 1, NC
         RNC   = REAL(IC(I))
         WC(I) = RNC/(WP(I)*RNORM)-1.0
         WRITE(18,3) CR(I), WC(I), IC(I)
      END DO

    3 FORMAT(2E14.3,I12)
      CLOSE(18)
      END
*---------------------------------------------------------------------------

****************************************************************************

*---------------------------------------------------------------------------
C     Analyseprogramm, RA - Anfangsradius, RD - Radienzuwachs,
C     LF - Datenfilelaenge
C     NA - Feld mit Anzahl der Paare in Analyse
C     NC - Laenge der Korrelationsfunktion       
*---------------------------------------------------------------------------
      SUBROUTINE ANA(r,N,isub,RA,RE,LF,IA,NC,rLbox,ilogi)

      real r(3,N)
      integer isub(LF)
      real    lRA,lRE,lDR

      DIMENSION IA(NC)

      print*
      write(*,*) 'Startradius, Endradius, Samplegroesse, Stuetzstellen'
      WRITE(*,*) RA,RE,LF,NC
      print*

      DO NR = 1, NC
         IA(NR) = 0
      END DO

      lRA = alog10(RA)
      lRE = alog10(RE)

c!$OMP PARALLEL DODEFAULT(SHARED)
c!$OMP+PRIVATE(I,J1,J,DR2,DR,lDR,IR)
      DO I = 1, LF-1
         !print*,I,LF-1
         J1 = I + 1
         DO J = J1, LF
            DR2 = ABST(rLbox,r(1,isub(I)),r(1,isub(J)),
     |                       r(2,isub(I)),r(2,isub(J)),
     |                       r(3,isub(I)),r(3,isub(J)))
            DR  = SQRT(DR2)
            lDR = alog10(DR)
            if(ilogi.EQ.0)then
               IR  = int((lDR-lRA)*real(NC)/real(lRE-lRA)) + 1
            else
               IR  = int((DR-RA)*real(NC)/real(RE-RA)) + 1
            end if

            IF ((IR.GE.1).AND.(IR.LE.NC)) IA(IR) = IA(IR) + 1 
 210     END DO
      END DO
      RETURN
      END
*---------------------------------------------------------------------------

****************************************************************************

*---------------------------------------------------------------------------
      SUBROUTINE NORM(rLbox,RA,RE,WP,NC,ilogi)

      DIMENSION WP(NC)

      real lRA,lRE

      lRA = alog10(RA)
      lRE = alog10(RE)

C   Dimension of the box
*   --------------------
      BOXL = rLbox
      VBOX = BOXL**3
      PI = 3.14157
      FAK = 4.0*PI/3.0 
      VA = FAK*RA**3
      DO I = 1, NC
         if(ilogi.EQ.0)then
            CR = 10.**(lRA + I * real(lRE-lRA)/real(NC))
         else
            CR = RA + I * real(RE-RA)/real(NC) 
         end if
         CV = FAK*CR**3
         WP(I) = (CV-VA)/VBOX
         VA=CV
      END DO
      END
*---------------------------------------------------------------------------

****************************************************************************

*---------------------------------------------------------------------------
      FUNCTION ABST(rLbox,X1,X2,Y1,Y2,Z1,Z2)

C     Euklidischer Abstand, Periodische Randbedingungen
C     Dimension of the box
      dx  = abs(x1 - x2)
      dy  = abs(y1 - y2)
      dz  = abs(z1 - z2)

      if(dx .GT. rLbox/2) dx=rLbox-dx
      if(dy .GT. rLbox/2) dy=rLbox-dy
      if(dz .GT. rLbox/2) dz=rLbox-dz

      abst = (dx**2 + dy**2 + dz**2)
      END
*------------------------------------------------------------------------------

******************************************************************************

*-----------------------------------------------------------------------------
      subroutine subsample(Nmax,NSUB,ISEED,isub)

      real ran3,random

      integer isub(NSUB)

      do i=1,NSUB
         random  = ran3(ISEED)
         isub(i) = int(random*(Nmax-1.)+1.)
      end do

      RETURN
      end
*-----------------------------------------------------------------------------

******************************************************************************
     
*-----------------------------------------------------------------------------
      FUNCTION RAN3(IDUM)
      SAVE
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END
*------------------------------------------------------------------------------

******************************************************************************

*-----------------------------------------------------------------------------
      SUBROUTINE get_name(name,i,j)
      CHARACTER*120 name
      DO i=1,120
         IF(name(i:i).NE.' ') GOTO 1
      ENDDO
      i=1
      j=1
      RETURN
 1    CONTINUE
      DO j=i,120
         IF(name(j:j).EQ.' ') GOTO 2
      ENDDO
      j=121
 2    j=j-1
      IF(name(i:i).EQ.'''') i=i+1
      IF(name(j:j).EQ.'''') j=j-1
      END
*------------------------------------------------------------------------------
      
*******************************************************************************

*------------------------------------------------------------------------------
      subroutine read_halos(datapath,data,r,N,NHALOS,MMIN,MMAX,icut)

      character*120 datapath,data,cdummy
      character*3   char_z
      real r(3,N),dummy,MMIN,MMAX
      
      parameter(MAXLINES = 5000000)

      call get_name(datapath,ip,jp)
      call get_name(data,i,j)

      print*,'reading file:  ',data(i:j)!//'_halos'
      print*

      !open(1,file=data(i:j)//'_halos')
      open(1,file=data(i:j))
      do i=1,1     ! AHF
        READ(1,*) cdummy
      enddo
      N=0
      do 200 iii=1,MAXLINES

c         READ(1,*,end=666) lkl,fMhires,rc1,rc2,rc3,v1,v2,v3,xMvir   ! AHF v0.0
         READ(1,*,end=666) ihalo,ihost,nsub,xMvir,lkl,rc1,rc2,rc3    ! AHF v1.0
c         READ(1,*,end=666) rc1,rc2,rc3                              ! legacy

         if(icut.EQ.0)then
            if(xMvir.GE.MMIN .AND. xMvir.LT.MMAX)then
               N=N+1
               r(1,N)=rc1
               r(2,N)=rc2       ! Koordinaten in Mpc/h
               r(3,N)=rc3
            end if
         else if(N.LT.NHALOS)then
            ! print*,N,NHALOS,rc1,rc2,rc3,lkl
            N=N+1
            r(1,N)=rc1
            r(2,N)=rc2          ! Koordinaten in Mpc/h
            r(3,N)=rc3
         end if
  200 continue
 666  CLOSE(1)

      RETURN
      end
*-----------------------------------------------------------------------------

******************************************************************************

*-----------------------------------------------------------------------------
      subroutine make_string(a,string)

      real a
      character*3 string
      write(string,123) a
 123  FORMAT(f3.1)
      do i=1,3
         if(string(i:i).EQ.' ') string(i:i)='0'
         if(string(i:i).EQ.'*') STOP 'Rotverschiebung nicht verwendbar'
      end do

      RETURN
      end
*-----------------------------------------------------------------------------
