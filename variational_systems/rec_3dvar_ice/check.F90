  SUBROUTINE CHECK(STATUS, line_num)

!USE MPI
USE PNETCDF
  IMPLICIT NONE

include 'mpif.h'

    INTEGER, INTENT ( IN) :: STATUS,line_num

    IF(STATUS /= NF_NOERR) THEN
      PRINT *, TRIM(NFmpi_strerror(STATUS))," line number=",line_num
     END IF
  
 IF(STATUS /= 0) PRINT *, "error line number=",line_num
!      CALL ABORT()
!    END IF
  END SUBROUTINE CHECK
