SUBROUTINE ICE_ANAM_AD

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE ANAMORPH ADJOINT
!                                                                      !
!-----------------------------------------------------------------------


  USE SET_KND
  USE GRD_STR
  USE ICE_TRANSF
  USE RUN, ONLY : NTSTEPS
  IMPLICIT NONE

  INTEGER(I4)    :: K,I,J
  REAL(R8)  :: WORK
  !! We use the last timestep read to determine the increments NTSTEP
  !GRD%SIC_AD is initialized in res_inc
        !SIC
        DO J=1,GRD%JM
           DO I=1,GRD%IM
              IF (GRD%MSK(I,J,1) .LT. 0.9_R8 ) CYCLE
              !WORK=GRD%SIC_AD(I,J)+GRD%SICB(I,J,NTSTEPS)
              !CALL sangoma_anamorphosis(1,qdim_sic,1,qua_ref_sic,SIC_QUANT(:,:,I,J),WORK,status)
              !GRD%SIC_AD(I,J) = (WORK-GRD%TRA_SICB(I,J,NTSTEPS))*GRD%MSK(I,J,1)
              IF (SIC_COEFF_TLAD(I,J) .NE. 0._R8) THEN
                   GRD%SIC_AD(I,J)=GRD%SIC_AD(I,J)+GRD%SIC_PHYS_AD(I,J)*SIC_COEFF_TLAD(I,J) 
                   GRD%SIC_PHYS_AD(I,J)=0._R8
               ELSE
                   GRD%SIC_AD(I,J)=0._R8
                   GRD%SIC_PHYS_AD(I,J)=0._R8
              ENDIF 
          ENDDO
        ENDDO
        
        !SIT

        DO J=1,GRD%JM
           DO I=1,GRD%IM
              IF (GRD%MSK(I,J,1) .LT. 0.9_R8  ) CYCLE
              !WORK=GRD%SIT_AD(I,J)+GRD%SITB(I,J,NTSTEPS)
              !CALL sangoma_anamorphosis(1,qdim_sit,1,qua_ref_sic,SIT_QUANT(:,:,I,J),WORK,status)
              !GRD%SIT_AD(I,J) = (WORK-GRD%TRA_SITB(I,J,NTSTEPS))*GRD%MSK(I,J,1)
              IF (SIT_COEFF_TLAD(I,J) .NE. 0._R8) THEN
                   GRD%SIT_AD(I,J)=GRD%SIT_AD(I,J)+GRD%SIT_PHYS_AD(I,J)*SIT_COEFF_TLAD(I,J)
                   GRD%SIT_PHYS_AD(I,J)=0._R8
               ELSE
                   GRD%SIT_AD(I,J)=0._R8
                   GRD%SIT_PHYS_AD(I,J)=0._R8
              ENDIF 
           ENDDO
        ENDDO

END SUBROUTINE ICE_ANAM_AD
