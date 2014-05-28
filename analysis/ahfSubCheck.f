      program NsatInHost

      parameter(MAXHALOS = 1000000)
      parameter(Rfac     = 1.d0)
      parameter(xMmin    = 1.d+01)

      real*8          Xc(MAXHALOS),Yc(MAXHALOS),Zc(MAXHALOS)
      real*8          Mvir(MAXHALOS),Rvir(MAXHALOS)
      integer*8       npart(MAXHALOS), idhalo(MAXHALOS)
      integer*8       hostHalo(MAXHALOS)
      integer*8       numSubStruct(MAXHALOS)
      real*8          BOX, Msat, Mhost, rr, gg, bb, ihost
      character*512   cdummy, infile

      write(*,'(A)')
     &'======================================================'
      write(*,'(A)')
     &'  read *_halos and count satellites of specific host'
      write(*,'(A)')
     &'======================================================'

      write(*,'(A,$)')' Please give name of *.AHF_halos file     => '
      read(*,'(A)')infile
      print*,infile(1:80)
      write(*,'(A,$)')' Please give boxsize [kpc/h]              => '
      read*,BOX
      print*,BOX
      write(*,'(A,$)')' Please give id of host halo              => '
      read*,ihost
      print*,ihost

      !==========================================
      ! open all relevant files
      !==========================================
      open(11,file=infile)
      open(12,file="subhaloesRvir.geom")
      open(13,file="subhaloesID.geom")


      !==========================================
      ! read in the full halo catalogue
      !==========================================
      read(11,*)cdummy
      do i=1,MAXHALOS
         read(11,*,end=666)
     &                     idhalo(i),
     &                     hostHalo(i),
     &                     numSubStruct(i),
     &                     Mvir(i),
     &                     npart(i),
     &                     Xc(i),Yc(i),Zc(i),
     &                     VXc,VYc,VZc,
     &                     Rvir(i)
         ! print*,idhalo(i)
      enddo
 666  nhalos = i-1
      close(11)
      print*
      print*,'o nhalos in file = ',nhalos
      print*

      ! find host halo using idhalo()
      do i=1,nhalos
         if(idhalo(i) .EQ. ihost)then
            Xhost    = Xc(i)
            Yhost    = Yc(i)
            Zhost    = Zc(i)
            Rhost    = Rvir(i)
            Mhost    = Mvir(i)
            goto 777
         endif
      enddo
 777  continue
      print*,'o locating subhaloes of host ',
     &        Xhost,Yhost,Zhost,Rhost,Mhost
      write(12,222)'s ',Xhost,Yhost,Zhost,Rhost,0,0,1,Mhost
      write(13,222)'s ',Xhost,Yhost,Zhost,Rhost,0,0,1,Mhost
 222  FORMAT(A,4f16.8,3i5,g16.8)

      !==========================================
      ! find subhaloes based upon Rvir criterion
      !==========================================
      Nsat = 0
      Msat = 0.0
      do i=1,nhalos

         if (idhalo(i) .NE. ihost)then
            dX = Xc(i)-Xhost
            dY = Yc(i)-Yhost
            dZ = Zc(i)-Zhost
            
            if(dX .GT.  BOX/2.) dX = dX-BOX
            if(dY .GT.  BOX/2.) dY = dY-BOX
            if(dZ .GT.  BOX/2.) dZ = dZ-BOX
            if(dX .LT. -BOX/2.) dX = dX+BOX
            if(dY .LT. -BOX/2.) dY = dY+BOX
            if(dZ .LT. -BOX/2.) dZ = dZ+BOX
            
            Dist = sqrt(dX**2 + dY**2 + dZ**2) + 0.5*Rvir(i)
            
            if( Dist     .LE. Rfac*Rhost .AND.
     &          Mvir(i)  .GT. xMmin      .AND.
     &          npart(i) .GE. 20              ) then
               
               rr = 1.
               gg = 0.
               bb = 0.
               
               write(12,111)'s ',
     &              Xc(i),Yc(i),Zc(i),Rvir(i),
     &              rr,gg,bb,
     &              '   ',
     &              Mvir(i),
     &              idhalo(i)

 111           FORMAT(A,4f16.8,3f12.8,A,g16.8,i22)
               
               Nsat = Nsat + 1
               Msat = Msat + Mvir(i)
            endif
         endif
      enddo ! nhalos
      close(12)

      print*,' Subhalo information -- Rvir criterion:'
      print*,'========================================'
      print*,'  - no. of satellites in host #',ihost,'       = ',Nsat
      print*,'  - total mass in satellites            = ',Msat
      print*,'  - fraction of host mass in satellites = ',Msat/Mhost
      print*


      !==========================================
      ! find subhaloes based upon AHF's hostHalo
      !==========================================
      Nsat = 0
      Msat = 0.0
      do i=1,nhalos

         if( hostHalo(i) .EQ. ihost .AND. 
     &       Mvir(i)     .GT. xMmin .AND. 
     &       npart(i)    .GE. 20          )then

            rr = 0.
            gg = 1.
            bb = 0.
            
            write(13,333)'s ',
     &           Xc(i),Yc(i),Zc(i),Rvir(i),
     &           rr,gg,bb,
     &           '   ',
     &           Mvir(i),
     &           idhalo(i)
            
 333        FORMAT(A,4f16.8,3f12.8,A,g16.8,i22)
            
            Nsat = Nsat + 1
            Msat = Msat + Mvir(i)
         endif
      enddo ! nhalos
      close(13)

      print*,' Subhalo information -- hostHalo criterion:'
      print*,'============================================'
      print*,'  - no. of satellites in host #',ihost,'       = ',Nsat
      print*,'  - total mass in satellites            = ',Msat
      print*,'  - fraction of host mass in satellites = ',Msat/Mhost
      print*


      STOP
      end
         
