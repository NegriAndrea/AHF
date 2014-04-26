!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ramses2gadget version 1.0                                                       !
! =========================                                                       !
! (c) 2009 by Timur Doumler                                                       !
! + some corrections and bug fixes by Robert Feldmann (2010)                      !
!                                                                                 !
! This tool reads data from cosmological AMR simulations (in RAMSES format),      !
! converts the hydrodynamic grid quantities into a gas particle distribution      !
! and writes all particles (dm, hydro and stars) to output files in GADGET        !
! format. This tool can be run in parallel using MPI.                             !
! The original purpose is to create data files readable by the AMIGA Halo Finder  !
! (AHF); but this tool can also be used for other purposes where GADGET-format    !
! is handy.                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! README:
!=========
! you need to call the code giving two additional arguments on the command line, i.e.
!
!./ramses2gadget -g ./output_00001
!
! - the first specifies the mode of operation
! - the second is the directory containing the RAMSES files
!
!
!
!
! 1. possible modes of operation are:
!
! -d
! only the dark matter particles found in the RAMSES files will be written to the GADGET file
!
!
! -g
! the RAMSES files contain gas information on a grid: for every AMR cell one gas particle will be generated
!
!
! -i
! as there can be far more gas cells than dark matter particles, this option ensures that the final GADGET
! file will contain of order the same gas and dark matter particles by interpolating gas cells (taking into
! account the AMR hierarchy)
!
!
!
!
! 2. RAMSES directory:
!
! the converter expects the RAMSES snapshot to be located in a directory, e.g. output_00001/
! which in turn will contain files of the format
!   info_00001.txt
!   amr_00001.out00001
!   hydro_00001.out00001
!   part_00001.out00001
!   amr_00001.out00002
!   hydro_00001.out00002
!   part_00001.out00002
!   amr_00001.out00003
!   hydro_00001.out00003
!   part_00001.out00003
!   ...
!
!
!
!
! Some more remarks:
!
! - the code is presently limited to 1024^3 particles of each type (see "nteger (kind=4)" for npart, ngaspart, etc. below)
! - for the MPI version the code needs to be compiled with -DMPI; the code itself figures out the best distribution
! - as the code uses pre-processor statements such as #ifdef, #define, etc. , this needs to be switched on as a compiler flag, e.g. "gfortran -x f95-cpp-input"



program ramses2gadget  

#ifdef MPI
  use mpi
#endif

  implicit none

  ! counters
  integer :: i, j, ilevel, idim, idm, istar, ifile, icpu, ivar, ind,col
  integer :: nfiles, ncpus, ndim, twotondim, firstfile, lastfile

  ! flags
  logical :: gas, stars, interpolate, refined(8), allrefined

  ! cosmological and physical parameters
  real(kind=8)            :: omega_m, omega_l, H0, h, a, gamma, boxsize, dbldummy

  ! unit conversion factors and parameters needed for that
  real(kind=8)            :: ramses_length_unit, ramses_dens_unit, ramses_time_unit
  real(kind=8)            :: gadget_length_unit, gadget_mass_unit, gadget_velocity_unit, &
       & gadget_etherm_unit
  real(kind=8)            :: xfactor, vfactor, mfactor, ufactor
  real(kind=8), parameter :: kpc_in_cm = 3.08d21 ! the factor is from ramses: dont use more digits!
  real(kind=8), parameter :: Msun_in_cm = 1.9891d33

  ! variables needed for file management
  character(len=128) :: dir_name, part_filename, amr_filename, hydro_filename, info_filename
  character(len=128) :: output_filename
  character(len=5)   :: dir_number,suffix
  character(len=5)   :: output_file_suffix
  integer :: part_file, amr_file, hydro_file, info_file, output_file, scratch_file !i/o unit numbers

  ! arrays that will hold the particle data (change integer kind to 8 for more than 2^31 particles!)
  integer (kind=4) :: npart, ngaspart, nramsespart, ndmpart, nstarpart, nsink, intdummy
  integer (kind=4) :: ngaspart_total, ndmpart_total, nstarpart_total
  real(kind=8), dimension(:,:), allocatable :: gaspart_pos, ramsespart_pos, dmpart_pos, starpart_pos
  real(kind=8), dimension(:,:), allocatable :: gaspart_vel, ramsespart_vel, dmpart_vel, starpart_vel
  real(kind=8), dimension(:),   allocatable :: gaspart_m,   ramsespart_m,   dmpart_m,   starpart_m
  integer(kind=4),dimension(:), allocatable :: gaspart_id,  ramsespart_id,  dmpart_id,  starpart_id
  real(kind=8), dimension(:),   allocatable :: gaspart_u
  real(kind=8), dimension(:),   allocatable :: ramsespart_age

  ! temporary variables for the particle data
  real(kind=8) :: x, y, z, vx, vy, vz, m, u, m_group
  integer(kind=4) :: id

  ! variables needed to read AMR data  
  integer           :: nboundary,nx,ny,nz,nlevelmax,ngrida,nvarh,ix,iy,iz
  real(kind=8)      :: dx
  character(len=80) :: ordering
  integer,      dimension(:,:),   allocatable :: son,ngridfile,ngridlevel,ngridbound
  real(kind=8), dimension(1:8,1:3)            :: xc
  real(kind=8), dimension(1:3)                :: xbound
  real(kind=8), dimension(:,:),   allocatable :: xg
  real(kind=8), dimension(:,:,:), allocatable :: var  

  ! little helpers for the gadget header    (change integer kind to 8 ffor more than 2^31 particles!)
  integer(kind=4),parameter :: gadget_fillheader(15) = 0
  integer(kind=4)           :: flagsfr, ngaspart_dummy, ndmpart_dummy, nstarpart_dummy

#ifdef MPI
  ! variables only needed if the code is run with mpi
  integer         :: mpi_ierr
  integer(kind=4) :: mpi_reduce_buffer, files_per_cpu
#endif


!!! START UP AND GET PARAMETERS FROM FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! initialise mpi environment
#ifdef MPI
  call MPI_INIT(mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, icpu, mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpus, mpi_ierr)
#else
  icpu  = 0
  ncpus = 1
#endif

  ! display hello message
  if (icpu==0) call hello

  ! read arguments from program call; set directory name and flags for gas and interpolation mode
  call read_args(icpu, dir_name, gas, interpolate)

  ! display status message
  if (icpu == 0) then
     if (.not.gas) write(*,*) 'Running in dark matter mode.'
     if (gas.and.(.not.interpolate)) write (*,*) 'Running with hydro in one-particle-per-cell mode.'
     if (gas.and.interpolate) write (*,*) 'Running with hydro in interpolation mode.'
  endif

  ! we need to open the info_...txt file in order to gather some important parameters.

  ! construct filename, inquire existence and assign I/O unit number
  ! dir_number=dir_name( index(dir_name,'output_')+7 : index(dir_name,'output_')+13 ) 
  dir_number='00036' 
  info_filename   = trim(dir_name) // '/info_'   // trim(dir_number) // '.txt'
  call inquire_file (info_filename)
  info_file   = 988

  ! in case we have gas, we need to open some more files and gather some more parameters...
  if (gas) then 
     amr_filename    = trim(dir_name) // '/amr_'    // trim(dir_number) // '.out00001'
     hydro_filename  = trim(dir_name) // '/hydro_'  // trim(dir_number) // '.out00001'
     call inquire_file (amr_filename)
     call inquire_file (hydro_filename)  
     amr_file    = 989
     hydro_file  = 990  
  end if

  ! open info_...txt file and read basic parameters and unit conversion factors
  open(unit=info_file,file=info_filename,form='formatted',status='old',action='read', &
       & position = 'rewind')
  read(info_file,'(13X,I11)')    nfiles
  read(info_file,'(13X,I11)')    ndim
  read(info_file,*)              ! skip levelmin
  read(info_file,*)              ! skip levelmax
  read(info_file,*)              ! skip ngridmax
  read(info_file,*)              ! skip nstep_coarse
  read(info_file,*)              ! skip empty line
  read(info_file,'(13X,E23.15)') boxsize ! called 'boxlen' in the file
  read(info_file,*)              ! skip time
  read(info_file,'(13X,E23.15)') a ! called 'aexp' in the file
  read(info_file,'(13X,E23.15)') H0
  read(info_file,'(13X,E23.15)') omega_m
  read(info_file,'(13X,E23.15)') omega_l
  read(info_file,*)              ! skip omega_k
  read(info_file,*)              ! skip omega_b
  read(info_file,'(13X,E23.15)') ramses_length_unit
  read(info_file,'(13X,E23.15)') ramses_dens_unit
  read(info_file,'(13X,E23.15)') ramses_time_unit
  read(info_file,*)              ! skip empty line
  read(info_file,'(14X,A80)')    ordering
  close(info_file)

  if (gas) then
     ! open first amr file (.out00001) and read some required quantities from the header
     open(unit=amr_file,file=amr_filename,status='old',form='unformatted',action='read', &
          & position = 'rewind')
     read(amr_file)  ! skip nfiles
     read(amr_file)  ! skip ndim
     read(amr_file)  nx,ny,nz
     read(amr_file)  nlevelmax
     read(amr_file)  ! skip ngridmax
     read(amr_file)  nboundary
     close(amr_file)

     ! open first hydro file (.out00001) and read gamma, this is needed for unit conversion
     open(unit=hydro_file,file=hydro_filename,status='old',form='unformatted',action='read', &
          & position = 'rewind')
     read(hydro_file) ! skip nfiles
     read(hydro_file) ! skip nvar
     read(hydro_file) ! skip ndim
     read(hydro_file) ! skip nlevelmax
     read(hydro_file) ! skip nboundary
     read(hydro_file) gamma
     close(hydro_file)
  end if

  ! display status message
  if (icpu==0) write(*,'(A,I4,A,I4,A)') ' Working with', ncpus, ' cpu(s) on', nfiles, ' fileset(s).'


!!! SET UP SOME MORE PARAMETERS, DO SOME CHECKING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! set some parameters
  twotondim = 2**ndim
  h         = 0.01 * H0

  ! define units of destination format (= GADGET standard units)
  ! (the numeric value must be set to the desired unit expressed in cgs units)
  gadget_length_unit          = kpc_in_cm / h
  gadget_mass_unit            = 1d10 * Msun_in_cm / h
  gadget_velocity_unit        = 1d5 * sqrt(a) ! sqrt(a) factor: see GADGET users guide
  if (gas) gadget_etherm_unit = 1d10 ! units of v^2, but without the sqrt(a) factor

  ! calculate unit conversion factors from ramses to gadget
  xfactor          = ramses_length_unit / gadget_length_unit / a ! [ RF ]
  vfactor          = (ramses_length_unit/ramses_time_unit) / gadget_velocity_unit
  mfactor          = ramses_dens_unit * (ramses_length_unit**3) / gadget_mass_unit
  if (gas) ufactor = (ramses_length_unit/ramses_time_unit)**2 / gadget_etherm_unit

  ! we will need boxsize in comoving kpc/h for the GADGET header
  boxsize = boxsize * xfactor

  if (gas) then
     ! allocate arrays (and do some other stuff) needed for AMR grid processing
     allocate(ngridfile(1:nfiles+nboundary,1:nlevelmax))
     allocate(ngridlevel(1:nfiles,1:nlevelmax))
     if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))
     xbound = (/dble(nx/2),dble(ny/2),dble(nz/2)/)
  end if

  ! set counters for *total* particle counts (needed for the GADGET header)
  ngaspart_total  = 0
  ndmpart_total   = 0
  nstarpart_total = 0

  ! check dimensions
  if(ndim/=3)then
     if (icpu==0) write(*,*)'Error in input files: not a 3D snapshot'
     call terminate(1)
  endif

  if (icpu==0) write(*,*) 'Reading parameters successful. Processing data in snapshot files.'



!!! START THE LOOP OVER SNAPSHIOT FILES/CPUS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPI
  ! tell every cpu which files it will work on
  files_per_cpu = ( nfiles+(ncpus-1) ) / ncpus
  firstfile = (icpu*files_per_cpu) + 1
  lastfile  = min ( (icpu+1)*files_per_cpu, nfiles )
#else
  ! in the serial case we simply loop over all files
  firstfile = 1
  lastfile  = nfiles
#endif

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
#endif

  ! start the loop over files
  do ifile = firstfile, lastfile

     ! construct file names

     write(suffix,'(i5.5)') ifile
     write(output_file_suffix,'(i5.1)') ifile - 1 ! to comply with GADGET file ending convention

     part_filename =   trim(dir_name) // '/part_'  // trim(dir_number) // '.out' // trim(suffix)
     info_filename =   trim(dir_name) // '/info_'  // trim(dir_number) // '.txt'
     amr_filename =    trim(dir_name) // '/amr_'   // trim(dir_number) // '.out' // trim(suffix)
     hydro_filename =  trim(dir_name) // '/hydro_' // trim(dir_number) // '.out' // trim(suffix)
     output_filename = trim(dir_name) // '/ramses2gadget_' // dir_number(3:5) // &
          & '.' // trim(adjustl(output_file_suffix))

     ! check if all the input files are there
     call inquire_file(part_filename)
     call inquire_file(info_filename)
     if (gas) then
        call inquire_file(amr_filename)
        call inquire_file(hydro_filename)
     end if

     ! assign unique I/O unit numbers
     part_file    = 10*ifile + 1
     info_file    = 10*ifile + 2
     output_file  = 10*ifile + 3
     if (gas) then
        amr_file     = 10*ifile + 4
        hydro_file   = 10*ifile + 5
        scratch_file = 10*ifile + 6
     end if

     ! Display a message saying which files we are working on.
     write (*,'(A,I4,A,I4,A)') ' cpu#', icpu, ', working on fileset#', ifile, ':'
     write (*,*)  '   reading files ' // trim(dir_name) // '/..._'// &
          & trim(dir_number) // '.out' // trim(suffix) // ' ;'
     write (*,*)  '   writing file  ' // trim(output_filename) // ' .'


     if (gas) then

!!! PREPARE AMR GRID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Open AMR file (again...) and skip header
        open(unit=amr_file,file=amr_filename,status='old',form='unformatted',action='read', &
             & position = 'rewind')
        do i=1,21
           read(amr_file)
        end do

        ! Read grid numbers
        read(amr_file)ngridlevel
        ngridfile(1:nfiles,1:nlevelmax)=ngridlevel
        read(amr_file)
        if(nboundary>0)then
           do i=1,2
              read(amr_file)
           end do
           read(amr_file)ngridbound
           ngridfile(nfiles+1:nfiles+nboundary,1:nlevelmax)=ngridbound
        endif
        read(amr_file)

        ! ROM: comment the single follwing line for old stuff
        read(amr_file)
        if(trim(ordering).eq.'bisection')then
           do i=1,5
              read(amr_file)
           end do
        else
           read(amr_file)
        endif
        read(amr_file)
        read(amr_file)
        read(amr_file)

        ! Open hydro file (again...) and read header    
        open(unit=hydro_file,file=hydro_filename,status='old',form='unformatted',action='read', &
             & position = 'rewind')
        read(hydro_file) ! skip nfiles
        read(hydro_file) nvarh ! called 'nvar' in RAMSES output routine
        read(hydro_file) ! skip ndim
        read(hydro_file) ! skip nlevelmax
        read(hydro_file) ! skip nboundary
        read(hydro_file) ! skip gamma     


!!! CREATING GAS PARTICLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! open scratch file for temporary storage of gas particle data.
        ! note: this is usually stored in /tmp, make sure there is enough disk space. 
        ! the scratch file will be deleted automatically after program execution.
        open(unit=scratch_file,status='scratch',form='unformatted',action='readwrite')

        ! this counter is very important for getting the gas particle
        ! number and generating their ids
        ngaspart = 0

        ! Loop over levels
        do ilevel=1,nlevelmax

           ! Geometry
           dx=0.5**ilevel
           do ind=1,twotondim
              iz=(ind-1)/4
              iy=(ind-1-4*iz)/2
              ix=(ind-1-2*iy-4*iz)
              xc(ind,1)=(dble(ix)-0.5d0)*dx
              xc(ind,2)=(dble(iy)-0.5d0)*dx
              xc(ind,3)=(dble(iz)-0.5d0)*dx
           end do

           ! Allocate work arrays
           ngrida=ngridfile(ifile,ilevel)        
           if(ngrida>0)then
              allocate(xg(1:ngrida,1:ndim))
              allocate(son(1:ngrida,1:twotondim))
              allocate(var(1:ngrida,1:twotondim,1:nvarh))
           endif

           ! Loop over domains
           do j=1,nboundary+nfiles

              ! Read AMR data
              if(ngridfile(j,ilevel)>0)then
                 read(amr_file) ! Skip grid index
                 read(amr_file) ! Skip next index
                 read(amr_file) ! Skip prev index
                 ! Read grid center
                 do idim=1,ndim
                    if(j.eq.ifile)then
                       read(amr_file)xg(:,idim)
                    else
                       read(amr_file)
                    endif
                 end do
                 read(amr_file) ! Skip father index
                 do ind=1,2*ndim
                    read(amr_file) ! Skip neighbour index
                 end do
                 ! Read son index
                 do ind=1,twotondim
                    if(j.eq.ifile)then
                       read(amr_file)son(:,ind)
                    else
                       read(amr_file)
                    end if
                 end do
                 ! Skip cpu map
                 do ind=1,twotondim
                    read(amr_file)
                 end do
                 ! Skip refinement map
                 do ind=1,twotondim
                    read(amr_file)
                 end do
              endif

              ! Read hydro data
              read(hydro_file)
              read(hydro_file)
              if(ngridfile(j,ilevel)>0)then
                 ! Read hydro variables
                 do ind=1,twotondim
                    do ivar=1,nvarh
                       if(j.eq.ifile)then
                          read(hydro_file)var(:,ind,ivar)
                       else
                          read(hydro_file)
                       end if
                    end do
                 end do
              end if
           end do !j=1,nboundary+nfiles ! end loop over domains

           if (ngrida>0) then

              ! loop over octs in grid
              do i=1,ngrida

                 ! loop over cells in oct; check which cells are refined
                 do ind=1,twotondim
                    refined(ind) = son(i,ind)>0 .and. ilevel<nlevelmax
                 end do

                 if (interpolate) then

!!!!!! in this case we interpolate between hydro cells !!!!!!

                    ! essentially what we do is: gather an appropriate group of hydro cells
                    ! and then place a particle in the centre-of-mass with appropriately
                    ! interpolated values of all quantities.

                    ! reset data variables
                    x=0; y=0; z=0; vx=0; vy=0; vz=0; m=0; u=0; m_group=0; allrefined=.true. 

                    ! loop over cells in oct
                    do ind=1,twotondim
                       if ( .not.refined(ind) ) then ! "group" means: non-refined cells in this oct

                          ! set a flag that we will use at least one cell from this oct
                          allrefined = .false.

                          ! get the mass in this cell, add to mass of group                 
                          m    = var(i,ind,1) * dx**3 ! mass = density * cell volume
                          m_group = m_group + m

                          ! we will place the particle at the centre-of-mass of these cells.
                          ! for this we weight everything with the mass:                       
                          x  = x  + m * ( (xg(i,1)+xc(ind,1)-xbound(1)) )
                          y  = y  + m * ( (xg(i,2)+xc(ind,2)-xbound(2)) )
                          z  = z  + m * ( (xg(i,3)+xc(ind,3)-xbound(3)) )
                          vx = vx + m * var(i,ind,2)
                          vy = vy + m * var(i,ind,3)
                          vz = vz + m * var(i,ind,4)
                          u  = u  + m * var(i,ind,5) / ((gamma-1d0)*var(i,ind,1)) ! u=p/(gamma-1)rho

                       end if ! ( .not.refined(ind) )
                    end do ! ind=1,twotondim


                    if (.not.allrefined) then ! meaning: if we have considered any cells in this oct

                       ! divide by group mass to get correct interpolated values
                       x =  x  / m_group
                       y =  y  / m_group
                       z =  z  / m_group
                       vx = vx / m_group
                       vy = vy / m_group
                       vz = vz / m_group
                       u  = u  / m_group

                       ! the interpolated mass equals the mass of the group
                       m  = m_group

                       ! now we generate a new gas particle.
                       ngaspart = ngaspart + 1   ! it is very important to update this counter!                 
                       id = -( (ngaspart*(nfiles-1))+ifile ) ! this is a *unique* id

                       ! write particle data into scratch file
                       write(scratch_file) x, y, z, vx, vy, vz, m, id, u

                    end if ! (.not.allrefined)

                 else ! (interpolate)

!!!!!! in this case we create particles on a simple per-cell basis !!!!!!

                    ! loop over cells in oct
                    do ind=1,twotondim
                       if ( .not.refined(ind) ) then   

                          ! a gas particle will be generated. this counter is very important:
                          ngaspart = ngaspart + 1

                          ! compute particle position:
                          x = (xg(i,1)+xc(ind,1)-xbound(1))
                          y = (xg(i,2)+xc(ind,2)-xbound(2))
                          z = (xg(i,3)+xc(ind,3)-xbound(3))

                          ! get the other variables:
                          vx = var(i,ind,2)
                          vy = var(i,ind,3)
                          vz = var(i,ind,4)
                          m  = var(i,ind,1) * dx**3 ! mass = density * cell volume
                          id = -( (ngaspart*(nfiles-1))+ifile ) ! this is a *unique* id
                          u  = var(i,ind,5) / ( (gamma-1d0)*var(i,ind,1) ) ! u = p/(gamma-1)rho

                          ! write particle data into scratch file
                          write(scratch_file) x, y, z, vx, vy, vz, m, id, u 

                       end if ! ( .not.refined )   
                    end do ! ind=1,twotondim                 
                 end if ! (interpolate)
              end do ! i=1,ngrida

              deallocate(xg,son,var)
           end if ! (ngrida>0)

        end do ! End loop over levels

        ! close ramses gas data files
        close(amr_file)
        close(hydro_file) 

        ! allocate gas particle data arrays  
        allocate ( gaspart_pos(1:3,1:ngaspart) )
        allocate ( gaspart_vel(1:3,1:ngaspart) )
        allocate ( gaspart_m(1:ngaspart) )
        allocate ( gaspart_id(1:ngaspart) )
        allocate ( gaspart_u(1:ngaspart) )

        ! rewind the scratch file that now holds the data, copy into data arrays, and delete file.
        rewind(scratch_file)

        do i=1,ngaspart
           read (scratch_file) gaspart_pos(1,i), gaspart_pos(2,i), gaspart_pos(3,i), &
                & gaspart_vel(1,i), gaspart_vel(2,i), gaspart_vel(3,i), &
                & gaspart_m(i), gaspart_id(i), gaspart_u(i)
        end do

        close(scratch_file) ! this also deletes the scratch file from disk

     end if ! (gas)

!!! READING DM + STAR PARTICLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     open(unit=part_file,file=part_filename,status='old',form='unformatted',action='read')

     ! read particle file header
     read(part_file) ! skip nfiles
     read(part_file) ! skip ndim
     read(part_file) nramsespart
     read(part_file) ! skip localseed
     read(part_file) nstarpart ! [ total number of stars in sim - RF ]
     read(part_file) ! skip mstar_tot
     read(part_file) ! skip mstar_lost
     read(part_file) nsink

     ! check for sink particles
     if (nsink /= 0) then
        if (icpu==0) write(*,*) 'Error in input files: sink particles encountered.'
        call terminate(2)
     end if

     ! check for different particle types; set flags and particle counters accordingly
     if (gas .and. nstarpart > 0) then
        stars    = .true.
        ndmpart  = nramsespart
 ! [ RF ] 2 lines deleted
 ! [ RF ] note: in this branch neither ndmpart nor nstarpart are yet set correctly
     else if (gas .and. nstarpart == 0) then
        stars    = .false.
        ndmpart  = nramsespart
        ! [ RF ] 1 line deleted
     else if ( (.not.gas) .and. (nstarpart==0) ) then
        stars    = .false.
        ndmpart  = nramsespart
        ! [ RF ] 1 line deleted
        ngaspart = 0
     else
        if(icpu==0) write(*,*) 'Error in input files: numbers of different particle types are wrong.'
        call terminate(3)
     endif

     npart    = nramsespart + ngaspart ! [ RF ]

     ! allocate data arrays for dm and star (=ramses) particles
     allocate ( ramsespart_pos(1:3,1:nramsespart) )
     allocate ( ramsespart_vel(1:3,1:nramsespart) )
     allocate ( ramsespart_m(1:nramsespart) )
     allocate ( ramsespart_id(1:nramsespart) )
     if (stars) allocate ( ramsespart_age(1:nramsespart) )

     ! read data from ramses particle file into these data arrays
     ! careful: ramses files contain pos and vel in column major order and separate records
     ! for x,y,z. for gadget output we need row major order and everything in one record.

     do i=1,3 ! loop over dimensions
        read (part_file) (ramsespart_pos(i,j), j=1,nramsespart)
     end do
     do i=1,3 ! loop over dimensions
        read (part_file) (ramsespart_vel(i,j), j=1,nramsespart)
     end do
     read (part_file) ramsespart_m
     read (part_file) ramsespart_id
     if (stars) then
        read (part_file) ! skip level
        read (part_file) ramsespart_age ! will use that to distinguish dm and star particles
     end if

     close(part_file) 


!!! UNIT CONVERSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! note: the conversion factors are defined above.

     ramsespart_pos = ramsespart_pos * xfactor  
     ramsespart_vel = ramsespart_vel * vfactor  
     ramsespart_m   = ramsespart_m   * mfactor

     if (gas) then
        gaspart_pos = gaspart_pos    * xfactor
        gaspart_vel = gaspart_vel    * vfactor
        gaspart_m   = gaspart_m      * mfactor
        gaspart_u   = gaspart_u      * ufactor
     end if


!!! OPEN OUTPUT FILE, WRITE GADGET HEADER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! see GADGET users guide to understand the structure of this header. it is quite simple.

     open(unit=output_file,file=output_filename,status='replace',form='unformatted',action='write')

     ! set star formation flag
     if (stars) then
        flagsfr = 1
     else
        flagsfr = 0
     end if

     ! write block identifier
     write(output_file) 'HEAD', 264_4

     ! now we must use a dirty workaround. the problem is that the gadget header requires
     ! the total number of particles in *all* files to be written into nall[6]. but at this
     ! point there is no way of getting this number (at least if the code is run in serial). 
     ! so instead we write dummy values (ngaspart = 666, ndmpart = 777, nstarpart = 888) 
     ! to these positions and overwrite them later by reopeining the output files in direct 
     ! access mode.

     ! [ RF ] if gas&stars -> need to compute nstarpart & ndmpart before header can be written
     if ((gas .eqv. .false.).or.(nstarpart .eq. 0)) then ! [ RF ]

        ! write header block:
        write(output_file) &
             &  ngaspart, ndmpart, 0_4, 0_4, nstarpart, 0_4,          & ! npart[6]
             &  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,             & ! massarr[6]
             &  a, (1/a - 1d0), flagsfr, 0_4,                         & ! a, z, flagsfr, flagfeedback                 
             &  666_4, 777_4, 0_4, 0_4, 888_4, 0_4,                   & ! nall[6] - see comment above
             &  0_4, nfiles, boxsize, omega_m, omega_l, H0, 0_4, 0_4, & ! flagcooling ... flagmetals
             &  0_4, 0_4, 0_4, 0_4, 0_4, 0_4,                         & ! nallhighw[6]
             &  0_4, gadget_fillheader                                  ! flagentropy, fillheader

        ! note: nallhighw (= more than 2^32 particles of one type) is currently not supported
        ! by AHF anyway. as soon as it will be added to AHF, this code will also be updated.

     end if ! [ RF ]


!!! WRITE DATA INTO OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! we use the gadget particle types 0 = gas, 1 = dark matter, 4 = stars

     ! there are three cases possible: we have only dm particles, or dm and gas particles,
     ! or dm, gas and stars particles.

     if (gas) then
        if (stars) then
           ! in this case we have dm, gas and stars. careful now: ramsespart contains dm AND star 
           ! particles mixed up! dm particles are distinguished by ramsespart_age == 0

           ! allocate new, separated data arrays for dm and stars
           allocate ( dmpart_pos(1:3,1:nramsespart) ) ! [ RF ]
           allocate ( dmpart_vel(1:3,1:nramsespart) ) ! [ RF ]
           allocate ( dmpart_m(1:nramsespart) )	      ! [ RF ]
           allocate ( dmpart_id(1:nramsespart) )	! [ RF ]
           allocate ( starpart_pos(1:3,1:nramsespart) ) ! [ RF ]
           allocate ( starpart_vel(1:3,1:nramsespart) ) ! [ RF ]
           allocate ( starpart_m(1:nramsespart) )	! [ RF ]
           allocate ( starpart_id(1:nramsespart) )      ! [ RF ]
           ! P.S. I really tried to sort it somehow without allocating more memory, 
           ! but didnt figure out how...

           ! go through the arrays with ramses data; sort and copy data
           idm=0 !dm counter
           istar=0 !star counter
           do i=1,nramsespart
              if ( ramsespart_age(i) == 0 ) then
                 ! particle is dark matter
                 idm = idm + 1
                 dmpart_pos(:,idm) = ramsespart_pos(:,i)
                 dmpart_vel(:,idm) = ramsespart_vel(:,i)
                 dmpart_m(idm) = ramsespart_m(i)
                 dmpart_id(idm) = ramsespart_id(i)
              else
                 ! particle is a star
                 istar = istar + 1
                 starpart_pos(:,istar) = ramsespart_pos(:,i) ! [ RF ]
                 starpart_vel(:,istar) = ramsespart_vel(:,i) ! [ RF ]
                 starpart_m(istar) = ramsespart_m(i)	     ! [ RF ]
                 starpart_id(istar) = ramsespart_id(i)	     ! [ RF ]
              end if
           end do

	   ndmpart=idm;     ! [ RF ]
     	   nstarpart=istar; ! [ RF ]

         ! check if particle count is ok
           if ( idm /= ndmpart .or. istar /= nstarpart .or. idm + istar /= nramsespart) then
              if (icpu==0) write(*,*) &
                   & 'Error in input files: counted wrong amount of dm and star particles'
              call terminate(4)            
           end if

	   ! >> [ RF ]
	   ! write header block:
           write(output_file) &
                &  ngaspart, ndmpart, 0_4, 0_4, nstarpart, 0_4,          & ! npart[6]
                &  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,             & ! massarr[6]
                &  a, (1/a - 1d0), flagsfr, 0_4,                         & ! a, z, flagsfr, flagfeedback                 
                &  666_4, 777_4, 0_4, 0_4, 888_4, 0_4,                   & ! nall[6] - see comment above
                &  0_4, nfiles, boxsize, omega_m, omega_l, H0, 0_4, 0_4, & ! flagcooling ... flagmetals
                &  0_4, 0_4, 0_4, 0_4, 0_4, 0_4,                         & ! nallhighw[6]
                &  0_4, gadget_fillheader                                  ! flagentropy, fillheader
	   write(*,*) 'stars=',stars,' nramsespart=',nramsespart,' npart=',npart
	   write(*,*) ' ndmpart=',ndmpart,' nstarpart=',nstarpart,' ngaspart=',ngaspart
    ! << [ RF ]

    ! deallocate ramses data arrays
           deallocate (ramsespart_pos,ramsespart_vel,ramsespart_m,ramsespart_id,ramsespart_age)

           ! write particle positions ! [ changed following lines to array slices - RF ]
           write(output_file) 'POS ', 3*8*npart+8
           write(output_file) gaspart_pos, (dmpart_pos(:,col),col=1,ndmpart), (starpart_pos(:,col),col=1,nstarpart)

           ! write particle velocities
           write(output_file) 'VEL ', 3*8*npart+8
           write(output_file) gaspart_vel, (dmpart_vel(:,col),col=1,ndmpart), (starpart_vel(:,col),col=1,nstarpart)

           ! write particle ids
           write(output_file) 'ID  ', 4*npart+8
           write(output_file) gaspart_id, (dmpart_id(col),col=1,ndmpart), (starpart_id(col),col=1,nstarpart)

           ! write particle masses
           write(output_file) 'MASS', 8*npart+8
           write(output_file) gaspart_m, (dmpart_m(col),col=1,ndmpart), (starpart_m(col),col=1,nstarpart)

           ! write gas particle thermal energies
           write(output_file) 'U   ', 8*ngaspart+8
           write(output_file) gaspart_u

           ! writing data is complete. close output file and deallocate remaining data arrays.
           close(output_file)
           deallocate (gaspart_pos,gaspart_vel,gaspart_m,gaspart_id,gaspart_u)
           deallocate (dmpart_pos,dmpart_vel,dmpart_m,dmpart_id)
           deallocate (starpart_pos,starpart_vel,starpart_m,starpart_id)

        else
           ! in this case we have dark matter and gas particles, but no stars

           ! write particle positions
           write(output_file) 'POS ', 3*8*npart+8
           write(output_file) gaspart_pos, ramsespart_pos

           ! write particle velocities
           write(output_file) 'VEL ', 3*8*npart+8
           write(output_file) gaspart_vel, ramsespart_vel

           ! write particle ids
           write(output_file) 'ID  ', 4*npart+8
           write(output_file) gaspart_id, ramsespart_id

           ! write particle masses
           write(output_file) 'MASS', 8*npart+8
           write(output_file) gaspart_m, ramsespart_m

           ! write gas particle thermal energies
           write(output_file) 'U   ', 8*ngaspart+8
           write(output_file) gaspart_u

           ! writing of the data is complete. close output file and deallocate all data arrays.
           close(output_file)
           deallocate (gaspart_pos,gaspart_vel,gaspart_m,gaspart_id,gaspart_u)
           deallocate (ramsespart_pos,ramsespart_vel,ramsespart_m,ramsespart_id)

        end if
     else
        ! in this case we have only dark matter particles

        ! write particle positions
        write(output_file) 'POS ', 3*8*npart+8
        write(output_file) ramsespart_pos

        ! write particle velocities
        write(output_file) 'VEL ', 3*8*npart+8
        write(output_file) ramsespart_vel

        ! write particle ids
        write(output_file) 'ID  ', 4*npart+8
        write(output_file) ramsespart_id

        ! write particle masses
        write(output_file) 'MASS', 8*npart+8
        write(output_file) ramsespart_m

        ! writing of the data is complete. close output file and deallocate all data arrays.
        close(output_file)
        deallocate (ramsespart_pos,ramsespart_vel,ramsespart_m,ramsespart_id)

     end if

     ! update counters for *total* particle counts (needed later for the GADGET header)
     ngaspart_total  = ngaspart_total  + ngaspart
     ndmpart_total   = ndmpart_total   + ndmpart
     nstarpart_total = nstarpart_total + nstarpart
     ! careful: if we run with MPI, this will only get the total particle count PER PROCESS.
     ! this means we will have to do an mpi_allreduce later to really get the total number.

     ! loop over snapshot files ends here
  end do ! ifile = firstfile, lastfile

  ! in case we run with MPI, get the total number of particles over all processes.
#ifdef MPI
  call MPI_ALLREDUCE(ngaspart_total,mpi_reduce_buffer, &
       & 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
  ngaspart_total = mpi_reduce_buffer

  call MPI_ALLREDUCE(ndmpart_total,mpi_reduce_buffer, &
       & 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
  ndmpart_total = mpi_reduce_buffer

  call MPI_ALLREDUCE(nstarpart_total,mpi_reduce_buffer, &
       & 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
  nstarpart_total = mpi_reduce_buffer
#endif

  ! the parallelised part ends here. shut down MPI environment.
#ifdef MPI
  call MPI_FINALIZE(mpi_ierr)
#endif

  if (icpu==0) write (*,*) 'Data processing completed successfully.'


!!! UPDATE GADGET HEADER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! now we have to reopen all output files in direct access mode and write the total number
  ! of particles for the three types in the nall[6] array. use the dummy values to check
  ! whether we are writing to the correct position.

  ! this is a very small task, so there is no need to parallelise here.
  if (icpu == 0) then

     write (*,*) 'Finishing GADGET-format files...'

     do ifile=1, nfiles

        ! open the written files again, but this time in direct access mode
        write(output_file_suffix,'(i5.1)') ifile - 1
        output_filename = trim(dir_name) // '/ramses2gadget_' // dir_number(3:5) // &
             & '.' // trim(adjustl(output_file_suffix))
        output_file  = 10*ifile + 7 ! again a unique I/O unit id
        open(unit=output_file,file=output_filename,status='old',form='unformatted', &
             & action='readwrite', access='direct', recl=4)

        ! display message that we are updating header.
        write (*,*) 'Writing GADGET-format header in file ' // trim(output_filename)

        ! check whether the dummy entries for ngaspart_total, ndmpart_total and nstarpart_total
        ! are at the correct position and if yes, overwrite with correct values.
        read( output_file,rec=30 )  ngaspart_dummy
        read( output_file,rec=31 )  ndmpart_dummy
        read( output_file,rec=34 )  nstarpart_dummy

        if (ngaspart_dummy==666 .and. ndmpart_dummy==777 .and. nstarpart_dummy==888) then
           write( output_file,rec=30 )  ngaspart_total
           write( output_file,rec=31 )  ndmpart_total
           write( output_file,rec=34 )  nstarpart_total
        else
	   write(*,*) 'ngaspart_dummy=',ngaspart_dummy, 'ndmpart_dummy=',ndmpart_dummy, & ! [ RF ]
                &          ' nstarpart_dummy=',nstarpart_dummy                                 ! [ RF ]
           write(*,'(A,I4,A)') 'file#', ifile, &
                & ': Error while writing GADGET header. The nall[6] entry will be wrong.'
        end if

        close(output_file)

     end do

     ! now we are done with everything. display message and end program.
     write (*,*) 'Generating GADGET-format files completed successfully.'

  end if ! (icpu == 0)

end program ramses2gadget


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_args(icpu, dir_name, gas, interpolate)

  ! this takes care of the arguments from the program call and sets mode flags

  implicit none

  integer,           intent(in)  :: icpu
  integer                        :: n, iargc
  character(len=2)               :: arg1
  character(len=128)             :: arg2
  character(len=128),intent(out) :: dir_name
  logical,           intent(out) :: gas, interpolate

  ! retrieve arguments    
  n = iargc()
  if (n /= 2) then
     call wrongcall(icpu)
  else
     ! first argument tells us in which mode to run the code
     call getarg(1,arg1)
     if (arg1 == '-g') then
        gas = .true.
        interpolate = .false.
     else if (arg1 == '-i') then
        gas = .true.
        interpolate = .true.
     else if (arg1 == '-d') then
        gas = .false.
        interpolate = .false.
     else
        call wrongcall
     end if

     ! second argument tells us the input directory path
     call getarg(2,arg2)
     dir_name = trim(arg2)  
  end if

  return  
end subroutine read_args

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrongcall(icpu)

  integer, intent(in) :: icpu

  if (icpu == 0) then
     write (*,*) 'usage: ./ramses2gadget -mode input_dir'
     write (*,*) 'possible modes:'
     write (*,*) ' -g   process dm and gas (and stars if present); one-particle-per-cell mode'
     write (*,*) ' -i   process dm and gas (and stars if present); interpolation mode'
     write (*,*) ' -d   process only dm (ignore gas and stars if present)'
  end if
  call terminate(5)

end subroutine wrongcall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hello

  write(*,*) 'ramses2gadget version 1.0 -- (c) 2009 by Timur Doumler'
  write(*,*) '======================================================'
#ifdef MPI
  write(*,*) 'Compiled with MPI.'
#else
  write(*,*) 'Compiled without MPI.'
#endif

end subroutine hello

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inquire_file(filename)

  ! this checks if a file is present

  implicit none

  character(LEN=128), intent(in)::filename
  logical :: ok

  inquire(file=filename, exist=ok) 
  if ( .not. ok ) then
     write (*,*) 'Error in input files: ' // trim(filename)// ' not found.'
     call terminate(6)
  endif

  return
end subroutine inquire_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine terminate(ierr)

  ! this subroutine terminates the whole program
  ! TODO: for some reason, the error message is not displayed anymore!!

  integer, intent(in) :: ierr
  integer :: mpi_ierr

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
#endif

  write (*,'(A,I3,A)') 'Terminating (Error code', ierr, ').'

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
  call MPI_FINALIZE(mpi_ierr)
#endif
  stop

!!! BANG !!!

end subroutine terminate

