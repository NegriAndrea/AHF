SUBROUTINE read_header(datadir,isnap)

! Originally by Tom Theuns

USE Hdf2GADGET_io
USE hdf5_wrapper

IMPLICIT NONE

CHARACTER*3:: snap
CHARACTER*120:: datadir
INTEGER:: file_handle, isnap
LOGICAL:: file_exists

! Construct file name

WRITE (snap,'(i3.3)') isnap
snapfile = TRIM(datadir) // TRIM(dir_prefix) // snap // '/' // &
 & TRIM(file_prefix) // snap // '.0.hdf5'
snapfile = TRIM(snapfile)
INQUIRE(FILE=snapfile,EXIST=file_exists)

IF (.NOT. file_exists) THEN
	WRITE (*,*) ' Error: snap file does not exist'
	WRITE (*,'(a)') TRIM(snapfile)
	STOP
ENDIF

CALL hdf5_open_file(file_handle,snapfile,readonly=.true.)
CALL hdf5_read_attribute(file_handle,'Header/NumPart_ThisFile',NumPart_ThisFile)
CALL hdf5_read_attribute(file_handle,'Header/NumPart_Total',NumPart_Total)
CALL hdf5_read_attribute(file_handle,'Header/NumPart_Total_HighWord', &
		& NumPart_Total_HighWord)
CALL hdf5_read_attribute(file_handle,'Header/NumFilesPerSnapshot', &
		& NumFilesPerSnapshot)
CALL hdf5_read_attribute(file_handle,'Header/Flag_Sfr',Flag_Sfr)
CALL hdf5_read_attribute(file_handle,'Header/Flag_Cooling',Flag_Cooling)
CALL hdf5_read_attribute(file_handle,'Header/Flag_StellarAge',Flag_StellarAge)
CALL hdf5_read_attribute(file_handle,'Header/Flag_Metals',Flag_Metals)
CALL hdf5_read_attribute(file_handle,'Header/Flag_Feedback',Flag_Feedback)
CALL hdf5_read_attribute(file_handle,'Header/MassTable',MassTable)
CALL hdf5_read_attribute(file_handle,'Header/Time_GYR',Time)
CALL hdf5_read_attribute(file_handle,'Header/Redshift',Redshift)
CALL hdf5_read_attribute(file_handle,'Header/ExpansionFactor',ExpansionFactor)
CALL hdf5_read_attribute(file_handle,'Header/BoxSize',BoxSize)
CALL hdf5_read_attribute(file_handle,'Header/Omega0',Omega0)
CALL hdf5_read_attribute(file_handle,'Header/OmegaBaryon',OmegaBaryon)
CALL hdf5_read_attribute(file_handle,'Header/OmegaLambda',OmegaLambda)
CALL hdf5_read_attribute(file_handle,'Header/HubbleParam',HubbleParam)
CALL hdf5_close_file(file_handle)

END SUBROUTINE read_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_gadget_file(datadir,isnap,ifile,n)

USE Hdf2GADGET_io
USE hdf5_wrapper

IMPLICIT NONE
LOGICAL :: file_exists
CHARACTER*1 :: is
CHARACTER*40 :: read_string
CHARACTER*120 :: param_file='param_file', datadir
CHARACTER(len=3) :: snap, ThisFile
INTEGER :: itype,ifile,file_handle,i,isnap,Ncum,Npart_all_types_this_file,n
INTEGER*4, ALLOCATABLE :: part_ids(:)
REAL*4, ALLOCATABLE :: pos_tf(:,:),m_tf(:),vel_tf(:,:),u_tf(:)

WRITE (snap,'(i3.3)') isnap
IF (ifile .LT. 10) THEN
	WRITE (ThisFile,'(I1)') ifile
ELSE IF (ifile .LT. 100) THEN
	WRITE (ThisFile,'(I2)') ifile
ELSE
	WRITE (ThisFile,'(I3)') ifile
ENDIF

snapfile = TRIM(datadir) // TRIM(dir_prefix) //snap // '/' //TRIM(file_prefix) &
		& // snap // '.' // TRIM(ThisFile) // '.hdf5'
INQUIRE(FILE=snapfile,EXIST=file_exists)
IF (.NOT. file_exists) THEN
	WRITE (*,*) ' error: snap file does not exist'
	WRITE (*,'(a)') TRIM(snapfile)
	STOP
ENDIF

CALL hdf5_open_file(file_handle,snapfile,readonly=.true.)
CALL hdf5_read_attribute(file_handle,'Header/NumPart_ThisFile',NumPart_ThisFile)

Npart_all_types_this_file = 0
DO itype=0,5
	Npart_all_types_this_file = Npart_all_types_this_file + &
		& NumPart_ThisFile(itype)
ENDDO

ALLOCATE(pos(3,Npart_all_types_this_file),vel(3,Npart_all_types_this_file), &
	& mass(Npart_all_types_this_file),part_ind(Npart_all_types_this_file))

PRINT *,'Arrays allocated'

! Reading the particle attributes from current file for all types

Ncum = 0 ! number of particles read so far

DO itype=0,5
	IF (NumPart_ThisFile(itype) .GT. 0) THEN

! Reading coordinates

		CALL hdf5_open_file(file_handle,snapfile,readonly=.true.)
		ALLOCATE(pos_tf(3,NumPart_ThisFile(itype)))
		WRITE(is,'(i1)') itype
		is=TRIM(is)
		read_string='PartType' // is // '/Coordinates'
		CALL hdf5_read_data(file_handle, read_string, pos_tf)
		pos(:,NCum+1:Ncum+NumPart_ThisFile(itype)) = pos_tf(:,:)
		DEALLOCATE(pos_tf)

! Reading velocities

		ALLOCATE(vel_tf(3,NumPart_ThisFile(itype)))
		read_string='PartType' // is // '/Velocity'
		CALL hdf5_read_data(file_handle, read_string, vel_tf)
		vel(:,NCum+1:Ncum+NumPart_ThisFile(itype)) = vel_tf(:,:)
		DEALLOCATE(vel_tf)

! Reading mass

		IF (MassTable(itype) .EQ. 0) THEN
			ALLOCATE(m_tf(NumPart_ThisFile(itype)))
			read_string='PartType' // is // '/Mass'
			CALL hdf5_read_data(file_handle,read_string,m_tf)
			mass(Ncum+1:Ncum+NumPart_ThisFile(itype))=m_tf
			DEALLOCATE(m_tf)
		ELSE
			mass(Ncum+1:Ncum+NumPart_ThisFile(itype)) = &
				& MassTable(itype)
		ENDIF

! Reading indexes

		ALLOCATE(part_ids(NumPart_ThisFile(itype)))
		read_string='PartType' // is // '/ParticleIDs'
		CALL hdf5_read_data(file_handle,read_string,part_ids)
		part_ind(NCum+1:Ncum+NumPart_ThisFile(itype))=part_ids
		DEALLOCATE(part_ids)

! Advancing index of particles read

		Ncum = Ncum + NumPart_ThisFile(itype)
		WRITE (*,'('' read '',I8,'' particles from file '',I3,'' of a total of '',I11)') NumPart_ThisFile(itype), ifile, NumPart_Total(itype)

	ENDIF
ENDDO

! Reading internal energy

ALLOCATE(u(NumPart_ThisFile(0)))
WRITE(is,'(i1)') 0
is=TRIM(is)
read_string='PartType' // is // '/InternalEnergy'
CALL hdf5_read_data(file_handle,read_string,u)

PRINT *,'Internal energies read'

CALL hdf5_close_file(file_handle)

n=Npart_all_types_this_file

PRINT *, 'Finished reading ',Npart_all_types_this_file,' particles from file ',&
		& ifile

END SUBROUTINE get_gadget_file

