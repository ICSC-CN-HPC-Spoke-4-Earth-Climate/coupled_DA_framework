SUBROUTINE ICE_ANAM

!-----------------------------------------------------------------------
! CALCULATE DENSITY INCREMENTS
!-----------------------------------------------------------------------


  USE SET_KND
  USE GRD_STR
  USE ICE_TRANSF
  USE RUN, ONLY : NTSTEPS
  IMPLICIT NONE

  INTEGER(I4)    :: I,J,status
  REAL(R8) :: WRKD
   !! We use the last timestep read to determine the increments

        !SIC
        DO J=1,GRD%JM
           DO I=1,GRD%IM
              IF (GRD%MSK(I,J,1)<0.9_R8 ) CYCLE
                IF (SIC_COEFF_TLAD(I,J) .NE. 0._R8) THEN
                     GRD%SIC(I,J)=GRD%SIC(I,J)*SIC_COEFF_TLAD(I,J) !FROM NOW ON SIC CAN BE USED AS PHYSICAL SPACE (RE-INITIALIZED AT THE BENINNING OF NEXT ITERATION 
                ELSE
                     GRD%SIC(I,J)=0._R8
                ENDIF
           ENDDO
        ENDDO

        !SIT
         DO J=1,GRD%JM
           DO I=1,GRD%IM
              IF (GRD%MSK(I,J,1)<0.9_R8 ) CYCLE
              IF (SIT_COEFF_TLAD(I,J) .NE. 0._R8) THEN
                   GRD%SIT(I,J)=GRD%SIT(I,J)*SIT_COEFF_TLAD(I,J)    !+GRD%TRA_SITB(I,J,NTSTEPS)
              ELSE
                   GRD%SIT(I,J)=0._R8
              ENDIF

           ENDDO
        ENDDO

END SUBROUTINE ICE_ANAM
