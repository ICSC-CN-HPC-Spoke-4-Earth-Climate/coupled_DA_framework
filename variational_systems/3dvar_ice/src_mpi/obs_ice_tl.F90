SUBROUTINE OBS_ICE_TL

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR SIC                                 !
!                                                                      !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR
 USE ICE_TRANSF
 USE TLAD_VARS
USE MYFRTPROF, ONLY : MYFRTPROF_WALL, MYFRTPROF_PRINT

 IMPLICIT NONE

 INTEGER(I4)   ::  KK,status
 REAL(KIND=R8) :: WKVAR

 CALL MYFRTPROF_WALL('OBS_SIC: SIC OBSERVATION OPERATOR',0)

if ( SIC%NO .GT. 0 ) THEN
#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(KK,status)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
CYOBS : DO KK=1,SIC%NO
  IF( .NOT. LL_4DVAR) THEN
        IF(SIC%FLC(KK).NE.1 ) CYCLE CYOBS
!        IF (SIC_COEFF_TLAD(SIC%IB(KK,1),SIC%JB(KK,1)) .NE. 0._R8) THEN
!           WKVAR =GRD%SIC(SIC%IB(KK,1),SIC%JB(KK,1)) 
!           SIC%INC(KK)=WKVAR/SIC_COEFF_TLAD(SIC%IB(KK,1),SIC%JB(KK,1))
!        ELSE
!            SIC%INC(KK)=0._R8
!        ENDIF
         SIC%INC(KK) = OSUM( SIC%PQ(KK,1:NPQ),GRD%SIC(SIC%IB(KK,1:NPQ),SIC%JB(KK,1:NPQ)) )
         
  ENDIF
ENDDO CYOBS
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

ENDIF

if ( SIT%NO .GT. 0 ) THEN
#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(KK,status)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
CYOBS2 : DO KK=1,SIT%NO
  IF( .NOT. LL_4DVAR) THEN
         IF(SIT%FLC(KK).NE.1 ) CYCLE CYOBS2
!        IF (SIT_COEFF_TLAD(SIT%IB(KK,1),SIT%JB(KK,1)) .NE. 0._R8) THEN
!                WKVAR=GRD%SIT(SIT%IB(KK,1),SIT%JB(KK,1))
!                SIT%INC(KK)=WKVAR/SIT_COEFF_TLAD(SIT%IB(KK,1),SIT%JB(KK,1))
!
!        ELSE
!                SIT%INC(KK)=0._R8
!        ENDIF
        SIT%INC(KK) = OSUM(SIT%PQ(KK,1:NPQ),GRD%SIT(SIT%IB(KK,1:NPQ),SIT%JB(KK,1:NPQ)) )

  ENDIF
ENDDO CYOBS2
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

ENDIF
 CALL MYFRTPROF_WALL('OBS_SIC: SIC OBSERVATION OPERATOR',1)

END SUBROUTINE OBS_ICE_TL
