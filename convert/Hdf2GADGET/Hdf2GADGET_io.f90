module Hdf2GADGET_io

IMPLICIT NONE
 
CHARACTER*200  :: dir_prefix = 'snapshot_', file_prefix = 'snap_', snapfile

INTEGER*4 :: NumPart_ThisFile(0:5), NumPart_Total(0:5), &
	& NumPart_Total_HighWord(0:5),NumFilesPerSnapshot,Flag_SFR, &
	& Flag_Cooling,Flag_StellarAge,Flag_Metals,Flag_Feedback,fill(16)
DOUBLE PRECISION :: MassTable(0:5), Time, Redshift, ExpansionFactor, &
	& BoxSize, Omega0, OmegaBaryon, OmegaLambda, HubbleParam
 
INTEGER*4, ALLOCATABLE :: part_ind(:),ptype(:)
REAL*4, ALLOCATABLE :: pos(:,:),vel(:,:),u(:),mass(:)
  
end module Hdf2GADGET_io
