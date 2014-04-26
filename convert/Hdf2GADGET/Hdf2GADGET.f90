PROGRAM Hdf2GADGET

USE Hdf2GADGET_io

IMPLICIT NONE

LOGICAL :: wm=.FALSE.
CHARACTER*3 :: istr
CHARACTER*4 :: header
CHARACTER*30, PARAMETER :: param_file='Hdf2GADGET_param_file.dat'
CHARACTER*200 :: out_file, out_file_temp, datadir
INTEGER :: isnap,ifile,itype,i,n
INTEGER*4, PARAMETER :: n_o_b=264

OPEN(1,FILE=param_file,STATUS='OLD',FORM='FORMATTED')
READ(1,*) datadir
READ(1,*) isnap
READ(1,*) out_file
CLOSE(1)

WRITE(istr,'(i3.3)') isnap
out_file=TRIM(ADJUSTL(out_file)) // 'snapshot_' // istr // '.'
PRINT *,out_file

CALL read_header(datadir,isnap)

DO ifile=0,NumFilesPerSnapshot-1
	CALL get_gadget_file(datadir,isnap,ifile,n)

	WRITE(istr,'(i3)') ifile
	out_file_temp=TRIM(ADJUSTL(out_file)) // ADJUSTL(istr)
	OPEN(1,FILE=out_file_temp,STATUS='UNKNOWN',FORM='UNFORMATTED')

	header='HEAD'
	WRITE(1) header,n_o_b
	WRITE(1) NumPart_ThisFile,MassTable,ExpansionFactor,Redshift,Flag_Sfr, &
		& Flag_Feedback,NumPart_Total,Flag_Cooling,NumFilesPerSnapshot,&
		& BoxSize,Omega0,OmegaLambda,HubbleParam,Flag_StellarAge, &
		& Flag_Metals,NumPart_Total_HighWord,fill

! Write positions

	header='POS'
	WRITE(1) header,n_o_b
	WRITE(1) pos
	DEALLOCATE(pos)	

! Write velocities

	header='VEL'
	WRITE(1) header,n_o_b
	WRITE(1) vel
	DEALLOCATE(vel)

! Write indexes

	header='ID'
	WRITE(1) header,n_o_b
	WRITE(1) part_ind
	DEALLOCATE(part_ind)

! Write masses

	DO itype=0,5
		IF ((MassTable(itype) .EQ. 0) .AND. &
			 &(NumPart_ThisFile(itype) .GT. 0)) wm=.TRUE.
	ENDDO
	IF (wm) THEN
	header='MASS'
	WRITE(1) header,n_o_b
	WRITE(1) mass
	DEALLOCATE(mass)
	ENDIF
	
	! Write internal energies
	header='U'

	WRITE(1) header,n_o_b
	WRITE(1) u
	DEALLOCATE(u)

ENDDO

PRINT *,'Done writing GADGET2 snapshot ',out_file
END PROGRAM Hdf2GADGET
