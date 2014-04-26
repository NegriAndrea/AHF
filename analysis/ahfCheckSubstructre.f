      program NsatInHost

      parameter(MAXHALOS = 1000000)
      parameter(Rfac     = 1.)
      parameter(xMmin    = 1.e+01)

      real          Xc(MAXHALOS),Yc(MAXHALOS),Zc(MAXHALOS)
      real          VXc(MAXHALOS),VYc(MAXHALOS),VZc(MAXHALOS)
      real          Mvir(MAXHALOS),Rvir(MAXHALOS)
      real          SigV(MAXHALOS)
      real          b(MAXHALOS),c(MAXHALOS)
      integer       npart(MAXHALOS)
      real          BOX, Msat, Mhost, rr, gg, bb
      character*512 cdummy, infile

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
      write(*,'(A,$)')' Please give id (C-notation) of host halo => '
      read*,ihost
      print*,ihost
      ! convert to FORTRAN index
      ihost = ihost + 1

      open(11,file=infile)
      open(12,file="satellites.geom")


      read(11,*)cdummy
      do i=1,MAXHALOS
         read(11,*,end=666)npart(i),vnpart,
     &                     Xc(i),Yc(i),Zc(i),
     &                     VXc(i),VYc(i),VZc(i),
     &                     Mvir(i),
     &                     Rvir(i),
     &        Vmax,Rmax,SigV(i),xlambda,xLx,xLy,xLz,
     &        a,Eax,Eay,Eaz,
     &        b(i),Ebx,Eby,Ebz,
     &        c(i),Ecx,Ecy,Ecz
      enddo
 666  nhalos = i-1
      close(11)
      print*
      print*,'o found nhalos     = ',nhalos

      Xhost    = Xc(ihost)
      Yhost    = Yc(ihost)
      Zhost    = Zc(ihost)
      VXhost   = VXc(ihost)
      VYhost   = VYc(ihost)
      VZhost   = VZc(ihost)
      Rhost    = Rvir(ihost)
      Mhost    = Mvir(ihost)
      SigVhost = SigV(ihost)
      
      !Rhost = 162./1000.

      ! write host halo
      write(12,222)'s ',Xhost,Yhost,Zhost,Rhost,0,0,1,Mhost,SigVhost
 222  FORMAT(A,4f16.8,3i5,g16.8,f16.8)
      Nsat = 0
      Msat = 0.0
      do i=1,nhalos

         if (i .NE. ihost)then
            dX = Xc(i)-Xhost
            dY = Yc(i)-Yhost
            dZ = Zc(i)-Zhost
            
            if(dX .GT.  BOX/2.) dX = BOX-dX
            if(dY .GT.  BOX/2.) dY = BOX-dY
            if(dZ .GT.  BOX/2.) dZ = BOX-dZ
            if(dX .LT. -BOX/2.) dX = BOX+dX
            if(dY .LT. -BOX/2.) dY = BOX+dY
            if(dZ .LT. -BOX/2.) dZ = BOX+dZ
            
            Dist = sqrt(dX**2 + dY**2 + dZ**2) + 0.5*Rvir(i)
            
            if(Dist .LE. Rfac*Rhost .AND. Mvir(i) .GT. xMmin 
     &           .AND. npart(i) .GE. 20) then
               
               rr = 1.
               gg = 0.
               bb = 0.
               
               Vsat = sqrt( (VXc(i)-VXhost)**2 + 
     &                      (VYc(i)-VYhost)**2 + 
     &                      (VZc(i)-VZhost)**2)
               
               
               write(12,111)'s ',
     &              Xc(i),Yc(i),Zc(i),Rvir(i),
     &              rr,gg,bb,
     &              '   ',
     &              Mvir(i),b(i),c(i),(1.-b(i)**2)/(1.-c(i)**2),
     &              Vsat,i-1

 111           FORMAT(A,4f16.8,3f12.8,A,g16.8,3f16.8,f16.8,i12)
               
               Nsat = Nsat + 1
               Msat = Msat + Mvir(i)
            endif

         endif

      enddo

      print*
      print*,'o no. of satellites in host #',ihost-1,'       = ',Nsat
      print*,'o total mass in satellites            = ',Msat
      print*,'o fraction of host mass in satellites = ',Msat/Mhost

      STOP
      end
         
