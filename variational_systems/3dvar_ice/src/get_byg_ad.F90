SUBROUTINE GET_BYG_AD

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE VERTICAL INTEGRAL OF BOUYANCY GRADIENT (ADJOINT)           !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR
  USE MYFRTPROF, ONLY : MYFRTPROF_WALL

  IMPLICIT NONE

  INTEGER(I4)    :: I,J,K
  REAL(R8), PARAMETER :: ZCOEFF = 9.81_R8 / 1020._R8

CALL MYFRTPROF_WALL('GET_BYG_AD: ADJOINT OF BUOYANCY GRADIENT',0)

    GRD%DNS(:,:,:) = 0._R8
    GRD%B_X(:,:,:) = 0._R8
    GRD%B_Y(:,:,:) = 0._R8

    DO J=2,GRD%JM
     DO I=2,GRD%IM
      GRD%BX(I,J) = GRD%BX(I,J)/GRD%DX(I,J)
      GRD%BY(I,J) = GRD%BY(I,J)/GRD%DY(I,J)
     ENDDO
    ENDDO

    DO K=GRD%KM,1,-1
     DO J=GRD%JM,2,-1
      DO I=GRD%IM,2,-1
       GRD%B_Y(I,J,K) = GRD%B_Y(I,J,K) + GRD%BY(I,J)*GRD%DZ(K) * GRD%MSK(I,J,K)*GRD%MSK(I,J-1,K)
       GRD%B_X(I,J,K) = GRD%B_X(I,J,K) + GRD%BX(I,J)*GRD%DZ(K) * GRD%MSK(I,J,K)*GRD%MSK(I-1,J,K)
      ENDDO
     ENDDO
    ENDDO

    DO K=GRD%KM,2,-1
     DO J=GRD%JM,2,-1
      DO I=GRD%IM,2,-1
       GRD%B_Y(I,J,K-1) = GRD%B_Y(I,J,K-1) + GRD%B_Y(I,J,K)
       GRD%DNS(I,J  ,K) = GRD%DNS(I,J  ,K) +               &
       & GRD%B_Y(I,J,K)*GRD%DZ(K)*GRD%MSK(I,J,K)*GRD%MSK(I,J-1,K)
       GRD%DNS(I,J-1,K) = GRD%DNS(I,J-1,K) -               &
       & GRD%B_Y(I,J,K)*GRD%DZ(K)*GRD%MSK(I,J,K)*GRD%MSK(I,J-1,K)
       GRD%B_Y(I,J,K)   = 0._R8

       GRD%B_X(I,J,K-1) = GRD%B_X(I,J,K-1) + GRD%B_X(I,J,K)
       GRD%DNS(I  ,J,K) = GRD%DNS(I  ,J,K) + &
       & GRD%B_X(I,J,K)*GRD%DZ(K)*GRD%MSK(I,J,K)*GRD%MSK(I-1,J,K)
       GRD%DNS(I-1,J,K) = GRD%DNS(I-1,J,K) -  &
       & GRD%B_X(I,J,K)*GRD%DZ(K)*GRD%MSK(I,J,K)*GRD%MSK(I-1,J,K)
       GRD%B_X(I,J,K)   = 0._R8
     ENDDO
    ENDDO
   ENDDO

   DO J=GRD%JM,2,-1
    DO I=GRD%IM,2,-1
     GRD%DNS(I,  J,1) = GRD%DNS(I,J  ,1) +               &
     & GRD%B_Y(I,J,1)*GRD%DZ(1)*GRD%MSK(I,J,1)*GRD%MSK(I,J-1,1)
     GRD%DNS(I,J-1,1) = GRD%DNS(I,J-1,1) -               &
     & GRD%B_Y(I,J,1)*GRD%DZ(1)*GRD%MSK(I,J,1)*GRD%MSK(I,J-1,1)
     GRD%B_Y(I,J,1)   = 0._R8

     GRD%DNS(I  ,J,1) = GRD%DNS(I  ,J,1) +               &
     & GRD%B_X(I,J,1)*GRD%DZ(1)*GRD%MSK(I,J,1)*GRD%MSK(I-1,J,1)
     GRD%DNS(I-1,J,1) = GRD%DNS(I-1,J,1) -             &
     & GRD%B_X(I,J,1)*GRD%DZ(1)*GRD%MSK(I,J,1)*GRD%MSK(I-1,J,1)
     GRD%B_X(I,J,1)   = 0._R8
    ENDDO
   ENDDO

   GRD%DNS = ZCOEFF * GRD%DNS
   CALL DENS_AD

CALL MYFRTPROF_WALL('GET_BYG_AD: ADJOINT OF BUOYANCY GRADIENT',1)
END SUBROUTINE GET_BYG_AD
