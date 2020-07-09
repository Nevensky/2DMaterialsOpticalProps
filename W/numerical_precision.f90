MODULE precision

	USE ISO_FORTRAN_ENV, ONLY : INT32, INT64, REAL32, REAL64 
	IMPLICIT NONE
	
		PRIVATE
		INTERGER, PARAMETER ::  i4 = INT32  &
								i8 = INT64  & 
								r4 = REAL32 & 
								r8 = REAL64
		
		PUBLIC :: i4, i8, r4, r8
	
END MODULE precision