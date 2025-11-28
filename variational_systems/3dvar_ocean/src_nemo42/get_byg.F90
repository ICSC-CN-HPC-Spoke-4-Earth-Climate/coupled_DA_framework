SUBROUTINE GET_BYG

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE VERTICAL INTEGRAL OF BOUYANCY GRADIENT                     !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR
  USE MYFRTPROF, ONLY :MYFRTPROF_WALL

  IMPLICIT NONE

  INTEGER(I4)    :: I,J,K
  REAL(R8), PARAMETER :: ZCOEFF = 9.81_R8 / 1020._R8

CALL MYFRTPROF_WALL('GET_BYG: BUOYANCY GRADIENT',0)

     CALL DENS_TL

     GRD%DNS = ZCOEFF * GRD%DNS

! BOUYANCY FORCE

    GRD%B_X(:,:,:) = 0._R8
    GRD%B_Y(:,:,:) = 0._R8

    DO J=2,GRD%JM
     DO I=2,GRD%IM
       GRD%B_X(I,J,1) = (GRD%DNS(I,J,1)-GRD%DNS(I-1,J,1))*&
       & GRD%DZ(1)*GRD%MSK(I,J,1)*GRD%MSK(I-1,J,1)
       GRD%B_Y(I,J,1) = (GRD%DNS(I,J,1)-GRD%DNS(I,J-1,1))*&
       & GRD%DZ(1)*GRD%MSK(I,J,1)*GRD%MSK(I,J-1,1)
     ENDDO
    ENDDO

    DO K=2,GRD%KM
     DO J=2,GRD%JM
      DO I=2,GRD%IM
       GRD%B_X(I,J,K) = GRD%B_X(I,J,K-1) + &
       & (GRD%DNS(I,J,K)-GRD%DNS(I-1,J,K))*GRD%DZ(K)*GRD%MSK(I,J,K)*GRD%MSK(I-1,J,K)
       GRD%B_Y(I,J,K) = GRD%B_Y(I,J,K-1) + &
       & (GRD%DNS(I,J,K)-GRD%DNS(I,J-1,K))*GRD%DZ(K)*GRD%MSK(I,J,K)*GRD%MSK(I,J-1,K)
      ENDDO
     ENDDO
    ENDDO

! VERICAL INTEGRAL OF BOUYANCY FORCE

    GRD%BX(:,:) = 0._R8
    GRD%BY(:,:) = 0._R8
    DO K=1,GRD%KM
     DO J=2,GRD%JM
      DO I=2,GRD%IM
       GRD%BX(I,J) = GRD%BX(I,J) + GRD%B_X(I,J,K)*GRD%DZ(K) * &
       & GRD%MSK(I,J,K)*GRD%MSK(I-1,J,K)
       GRD%BY(I,J) = GRD%BY(I,J) + GRD%B_Y(I,J,K)*GRD%DZ(K) * &
       & GRD%MSK(I,J,K)*GRD%MSK(I,J-1,K)
      ENDDO
     ENDDO
    ENDDO

    DO J=2,GRD%JM
     DO I=2,GRD%IM
      GRD%BX(I,J) = GRD%BX(I,J)/GRD%DX(I,J)
      GRD%BY(I,J) = GRD%BY(I,J)/GRD%DY(I,J)
     ENDDO
    ENDDO

CALL MYFRTPROF_WALL('GET_BYG: BUOYANCY GRADIENT',1)
END SUBROUTINE GET_BYG
