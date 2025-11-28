MODULE EOSINSITU
!
! ADAPTED FROM NEMOVAR
!
! NOTE : SPATIALLY-VARYING DEPTH (FOR PARTIAL STEPS)
! IS REPLACED BY SPATIALLY-HOMOGENEOUS DEPTH
!

USE SET_KND
USE GRD_STR

IMPLICIT NONE

CONTAINS

   SUBROUTINE EOS_INSITU(JPI,JPJ,JPK,DEP,PT,PS,PRD)
      IMPLICIT NONE
      INTEGER(I4), INTENT(IN) ::   JPI, JPJ, JPK
      REAL(R8), DIMENSION(JPI,JPJ,JPK), INTENT(IN   ) ::   PT   ! 1 : POTENTIAL TEMPERATURE  [CELCIUS]
      REAL(R8), DIMENSION(JPI,JPJ,JPK), INTENT(IN   ) ::   PS   ! 2 : SALINITY [PSU]
      REAL(R8), DIMENSION(        JPK), INTENT(IN   ) ::   DEP
      REAL(R8), DIMENSION(JPI,JPJ,JPK), INTENT(  OUT) ::   PRD   ! IN SITU DENSITY            [-]
      REAL(R8), DIMENSION(JPI,JPJ,JPK)                ::   ZWS
      !!
      INTEGER(I4) ::   JI, JJ, JK           ! DUMMY LOOP INDICES
      REAL(R8) ::   ZT , ZS , ZH , ZSR   ! LOCAL SCALARS
      REAL(R8) ::   ZR1, ZR2, ZR3, ZR4   !   -      -
      REAL(R8) ::   ZRHOP, ZE, ZBW, ZB   !   -      -
      REAL(R8) ::   ZD , ZC , ZAW, ZA    !   -      -
      REAL(R8) ::   ZB1, ZA1, ZKW, ZK0   !   -      -
      !!----------------------------------------------------------------------
      !
!CDIR NOVERRCHK
         ZWS(:,:,:) = SQRT( ABS( PS(:,:,:) ) )
         !
         DO JK = 1, JPK-1
            DO JJ = 1, JPJ
               DO JI = 1, JPI
                  ZT = PT   (JI,JJ,JK)
                  ZS = PS   (JI,JJ,JK)
                  ZH = DEP(JK)        ! DEPTH
                  ZSR= ZWS   (JI,JJ,JK)        ! SQUARE ROOT SALINITY
                  !
                  ! COMPUTE VOLUMIC MASS PURE WATER AT ATM PRESSURE
                  ZR1= ( ( ( ( 6.536332E-9_R8  *ZT - 1.120083E-6_R8 )*ZT + 1.001685E-4_R8 )*ZT   &
                     &        -9.095290E-3_R8 )*ZT + 6.793952E-2_R8 )*ZT +  999.842594_R8
                  ! SEAWATER VOLUMIC MASS ATM PRESSURE
                  ZR2= ( ( ( 5.3875E-9_R8*ZT-8.2467E-7_R8 ) *ZT+7.6438E-5_R8 ) *ZT        &
                     &                      -4.0899E-3_R8 ) *ZT+0.824493_R8
                  ZR3= ( -1.6546E-6_R8*ZT+1.0227E-4_R8 )    *ZT-5.72466E-3_R8
                  ZR4= 4.8314E-4_R8
                  !
                  ! POTENTIAL VOLUMIC MASS (REFERENCE TO THE SURFACE)
                  ZRHOP= ( ZR4*ZS + ZR3*ZSR + ZR2 ) *ZS + ZR1
                  !
                  ! ADD THE COMPRESSION TERMS
                  ZE = ( -3.508914E-8_R8*ZT-1.248266E-8_R8 ) *ZT-2.595994E-6_R8
                  ZBW= (  1.296821E-6_R8*ZT-5.782165E-9_R8 ) *ZT+1.045941E-4_R8
                  ZB = ZBW + ZE * ZS
                  !
                  ZD = -2.042967E-2_R8
                  ZC =   (-7.267926E-5_R8*ZT+2.598241E-3_R8 ) *ZT+0.1571896_R8
                  ZAW= ( ( 5.939910E-6_R8*ZT+2.512549E-3_R8 ) *ZT-0.1028859_R8 ) *ZT - 4.721788_R8
                  ZA = ( ZD*ZSR + ZC ) *ZS + ZAW
                  !
                  ZB1=   (-0.1909078_R8*ZT+7.390729_R8 )        *ZT-55.87545_R8
                  ZA1= ( ( 2.326469E-3_R8*ZT+1.553190_R8)       *ZT-65.00517_R8 ) *ZT+1044.077_R8
                  ZKW= ( ( (-1.361629E-4_R8*ZT-1.852732E-2_R8 ) *ZT-30.41638_R8 ) *ZT + 2098.925_R8 ) *ZT+190925.6_R8
                  ZK0= ( ZB1*ZSR + ZA1 )*ZS + ZKW
                  !
                  ! MASKED IN SITU DENSITY ANOMALY
                  PRD(JI,JJ,JK) = (  ZRHOP / (  1.0_R8 - ZH / ( ZK0 - ZH * ( ZA - ZH * ZB ) )  )    &
                     &   ) * GRD%MSK(JI,JJ,JK)
               END DO
            END DO
         END DO
         !
   END SUBROUTINE EOS_INSITU

   SUBROUTINE EOS_INSITU_TL(JPI,JPJ,JPK,DEP,PTEM,PSAL,PTEM_TL,PSAL_TL,PRD_TL)
      IMPLICIT NONE
      INTEGER(I4), INTENT(IN) ::   JPI, JPJ, JPK           
      REAL(R8), DIMENSION(        JPK), INTENT(IN   ) ::   DEP
      REAL(R8), DIMENSION(JPI,JPJ,JPK)                ::   ZWS 
      REAL(R8), DIMENSION(JPI,JPJ,JPK), INTENT( IN ) ::   &
         & PTEM,                 &  ! POTENTIAL TEMPERATURE
         & PSAL,                 &  ! SALINITY
         & PTEM_TL,              &  ! TL OF POTENTIAL TEMPERATURE
         & PSAL_TL                  ! TL OF SALINITY
      REAL(R8), DIMENSION(JPI,JPJ,JPK), INTENT( OUT ) ::   &
         & PRD_TL                   ! TL OF POTENTIAL DENSITY (SURFACE REFERENCED)
      !! * LOCAL DECLARATIONS
      INTEGER ::  JI, JJ, JK      ! DUMMY LOOP INDICES
      REAL(R8) ::   &             ! TEMPORARY SCALARS
         ZT, ZS, ZH, ZSR, ZR1, ZR2, ZR3, ZR4, ZRHOP, ZE, ZBW,   &
         ZB, ZD, ZC, ZAW, ZA, ZB1, ZA1, ZKW, ZK0,               &
         ZTTL, ZSTL, ZHTL, ZSRTL, ZR1TL, ZR2TL, ZR3TL,          &
         ZR4TL, ZRHOPTL, ZETL, ZBWTL,                           &
         ZBTL, ZDTL, ZCTL, ZAWTL, ZATL, ZB1TL, ZA1TL,           &
         ZKWTL, ZK0TL, ZPES, ZRDC1, ZRDC2, ZEPS,                &
         ZMASK
      !!----------------------------------------------------------------------
      ZEPS = 1.E-14
!CDIR NOVERRCHK
         ZWS(:,:,:) = SQRT( ABS( PSAL(:,:,:) ) )

         DO JK = 1, JPK-1                                 ! HORIZONTAL SLAB
            !                                             ! ===============
            DO JJ = 1, JPJ
               DO JI = 1, JPI
                  ZT = PTEM(JI,JJ,JK)
                  ZS = PSAL(JI,JJ,JK)
                  ! DEPTH
                  ZH = DEP(JK)        ! DEPTH
                  ! SQUARE ROOT SALINITY
                  ZSR= ZWS(JI,JJ,JK)
                  ! COMPUTE VOLUMIC MASS PURE WATER AT ATM PRESSURE
                  ZR1= ( ( ( ( 6.536332E-9*ZT-1.120083E-6 )*ZT+1.001685E-4)*ZT &
                     -9.095290E-3 )*ZT+6.793952E-2 )*ZT+999.842594
                  ! SEAWATER VOLUMIC MASS ATM PRESSURE
                  ZR2= ( ( ( 5.3875E-9*ZT-8.2467E-7 ) *ZT+7.6438E-5 ) *ZT   &
                     -4.0899E-3 ) *ZT+0.824493
                  ZR3= ( -1.6546E-6*ZT+1.0227E-4 ) *ZT-5.72466E-3
                  ZR4= 4.8314E-4

                  ! POTENTIAL VOLUMIC MASS (REFERENCE TO THE SURFACE)
                  ZRHOP= ( ZR4*ZS + ZR3*ZSR + ZR2 ) *ZS + ZR1

                  ! ADD THE COMPRESSION TERMS
                  ZE = ( -3.508914E-8*ZT-1.248266E-8 ) *ZT-2.595994E-6
                  ZBW= (  1.296821E-6*ZT-5.782165E-9 ) *ZT+1.045941E-4
                  ZB = ZBW + ZE * ZS

                  ZD = -2.042967E-2
                  ZC =   (-7.267926E-5*ZT+2.598241E-3 ) *ZT+0.1571896
                  ZAW= ( ( 5.939910E-6*ZT+2.512549E-3 ) *ZT-0.1028859 ) *ZT - 4.721788
                  ZA = ( ZD*ZSR + ZC ) *ZS + ZAW

                  ZB1=   (-0.1909078*ZT+7.390729 ) *ZT-55.87545
                  ZA1= ( ( 2.326469E-3*ZT+1.553190)*ZT-65.00517 ) *ZT+1044.077
                  ZKW= ( ( (-1.361629E-4*ZT-1.852732E-2 ) *ZT-30.41638 ) *ZT + 2098.925 ) *ZT+190925.6
                  ZK0= ( ZB1*ZSR + ZA1 )*ZS + ZKW

                  ! TANGENT LINEAR PART

                  ZTTL = PTEM_TL(JI,JJ,JK)
                  ZSTL = PSAL_TL(JI,JJ,JK)

                  ZSRTL= ( 1.0 / MAX( 2.*ZSR, ZEPS ) ) &
                     &   * GRD%MSK(JI,JJ,JK)                * ZSTL

                  ZR1TL= ( ( ( (  5.*6.536332E-9   * ZT &
                     &           -4.*1.120083E-6 ) * ZT &
                     &           +3.*1.001685E-4 ) * ZT &
                     &           -2.*9.095290E-3 ) * ZT &
                     &           +   6.793952E-2        ) * ZTTL

                  ZR2TL= ( ( (    4.*5.3875E-9   * ZT &
                     &           -3.*8.2467E-7 ) * ZT &
                     &           +2.*7.6438E-5 ) * ZT &
                     &           -   4.0899E-3        ) * ZTTL

                  ZR3TL= (       -2.*1.6546E-6   * ZT &
                     &           +   1.0227E-4        ) * ZTTL

                  ZRHOPTL=                                  ZR1TL &
                     &     + ZS                           * ZR2TL &
                     &     + ZSR * ZS                     * ZR3TL &
                     &     + ZR3 * ZS                     * ZSRTL &
                     &     + (  2. * ZR4 * ZS + ZR2              &
                     &        + ZR3 * ZSR           )     * ZSTL

                  ZETL = (       -2.*3.508914E-8   * ZT &
                     &           -   1.248266E-8        ) * ZTTL

                  ZBWTL= (        2.*1.296821E-6   * ZT &
                     &           -   5.782165E-9        ) * ZTTL

                  ZBTL =                                    ZBWTL &
                     &    + ZS                            * ZETL  &
                  &       + ZE                            * ZSTL

                  ZCTL = (       -2.*7.267926E-5   * ZT &
                     &           +   2.598241E-3        ) * ZTTL

                  ZAWTL= ( (      3.*5.939910E-6   * ZT &
                     &           +2.*2.512549E-3 ) * ZT &
                     &           -   0.1028859          ) * ZTTL

                  ZATL =                                    ZAWTL &
                  &      + ZD * ZS                        * ZSRTL &
                  &      + ZS                             * ZCTL  &
                  &      + ( ZD * ZSR + ZC )              * ZSTL

                  ZB1TL= (       -2.*0.1909078     * ZT &
                     &           +   7.390729           ) * ZTTL

                  ZA1TL= ( (      3.*2.326469E-3   * ZT &
                     &           +2.*1.553190    ) * ZT &
                     &           -   65.00517           ) * ZTTL

                  ZKWTL= ( ( (   -4.*1.361629E-4   * ZT &
                     &           -3.*1.852732E-2 ) * ZT &
                     &           -2.*30.41638    ) * ZT &
                     &           +   2098.925           ) * ZTTL

                  ZK0TL=                                    ZKWTL &
                     &  + ZB1 * ZS                        * ZSRTL &
                     &  + ZS  * ZSR                       * ZB1TL &
                     &  + ZS                              * ZA1TL &
                     &  + ( ZB1 * ZSR + ZA1 )             * ZSTL

                  ! MASKED IN SITU DENSITY

                  ZRDC1 = 1.0 / ( ZK0 - ZH * ( ZA - ZH * ZB ) )
                  ZRDC2 = 1.0 / ( 1.0 - ZH * ZRDC1 )

                  PRD_TL(JI,JJ,JK) = GRD%MSK(JI,JJ,JK) * ZRDC2 * &
                     &               (                              ZRHOPTL &
                     &                 - ZRDC2 * ZH * ZRDC1**2 * ZRHOP      &
                     &                   * (                        ZK0TL   &
                     &                       - ZH * (               ZATL    &
                     &                                - ZH        * ZBTL ) ) )
                    

               END DO

            END DO
            !                                             ! ===============
         END DO                                           !   END OF SLAB
         !                                                ! ===============

END SUBROUTINE EOS_INSITU_TL

SUBROUTINE EOS_INSITU_AD(JPI,JPJ,JPK,DEP,PTEM,PSAL,PTEM_AD,PSAL_AD,PRD_AD)
      IMPLICIT NONE
      INTEGER(I4), INTENT(IN) ::   JPI, JPJ, JPK
      REAL(R8), DIMENSION(        JPK), INTENT(IN   ) ::   DEP
      REAL(R8), DIMENSION(JPI,JPJ,JPK), INTENT( IN ) ::   &
         PTEM,                 &  ! POTENTIAL TEMPERATURE
         PSAL                     ! SALINITY
      REAL(R8), DIMENSION(JPI,JPJ,JPK), INTENT( INOUT ) ::   &
         PTEM_AD,              &  ! POTENTIAL TEMPERATURE
         PSAL_AD                  ! SALINITY
      REAL(R8), DIMENSION(JPI,JPJ,JPK), INTENT( INOUT ) ::   &
         PRD_AD                   ! POTENTIAL DENSITY (SURFACE REFERENCED)
      !! * LOCAL DECLARATIONS
      INTEGER ::  JI, JJ, JK      ! DUMMY LOOP INDICES
      REAL(R8) ::   &             ! TEMPORARY SCALARS
         ZT, ZS, ZH, ZSR, ZR1, ZR2, ZR3, ZR4, ZRHOP, ZE, ZBW,   &
         ZB, ZD, ZC, ZAW, ZA, ZB1, ZA1, ZKW, ZK0,               &
         ZTAD, ZSAD, ZHAD, ZSRAD, ZR1AD, ZR2AD, ZR3AD,          &
         ZR4AD, ZRHOPAD, ZEAD, ZBWAD,                           &
         ZBAD, ZDAD, ZCAD, ZAWAD, ZAAD, ZB1AD, ZA1AD,           &
         ZKWAD, ZK0AD, ZPES, ZRDC1, ZRDC2, ZEPS,                &
         ZMASK
      REAL(R8), DIMENSION(JPI,JPJ,JPK) :: ZWS

      ZTAD    = 0.0_R8
      ZSAD    = 0.0_R8
      ZHAD    = 0.0_R8
      ZSRAD   = 0.0_R8
      ZR1AD   = 0.0_R8
      ZR2AD   = 0.0_R8
      ZR3AD   = 0.0_R8
      ZR4AD   = 0.0_R8
      ZRHOPAD = 0.0_R8
      ZEAD    = 0.0_R8
      ZBWAD   = 0.0_R8
      ZBAD    = 0.0_R8
      ZDAD    = 0.0_R8
      ZCAD    = 0.0_R8
      ZAWAD   = 0.0_R8
      ZAAD    = 0.0_R8
      ZB1AD   = 0.0_R8
      ZA1AD   = 0.0_R8
      ZKWAD   = 0.0_R8
      ZK0AD   = 0.0_R8

      ZEPS = 1.E-14

      ZWS(:,:,:) = SQRT( ABS( PSAL(:,:,:) ) )

         !                                                ! ===============
         DO JK = 1, JPK-1                                 ! HORIZONTAL SLAB
            !                                             ! ===============
            DO JJ = 1, JPJ
               DO JI = 1, JPI
                  ZT = PTEM(JI,JJ,JK)
                  ZS = PSAL(JI,JJ,JK)
                  ! DEPTH
                  ZH = DEP(JK)        ! DEPTH
                  ! SQUARE ROOT SALINITY
                  ZSR= ZWS(JI,JJ,JK)
                  ! COMPUTE VOLUMIC MASS PURE WATER AT ATM PRESSURE
                  ZR1= ( ( ( ( 6.536332E-9*ZT-1.120083E-6 )*ZT+1.001685E-4)*ZT &
                     -9.095290E-3 )*ZT+6.793952E-2 )*ZT+999.842594
                  ! SEAWATER VOLUMIC MASS ATM PRESSURE
                  ZR2= ( ( ( 5.3875E-9*ZT-8.2467E-7 ) *ZT+7.6438E-5 ) *ZT   &
                     -4.0899E-3 ) *ZT+0.824493
                  ZR3= ( -1.6546E-6*ZT+1.0227E-4 ) *ZT-5.72466E-3
                  ZR4= 4.8314E-4

                  ! POTENTIAL VOLUMIC MASS (REFERENCE TO THE SURFACE)
                  ZRHOP= ( ZR4*ZS + ZR3*ZSR + ZR2 ) *ZS + ZR1

                  ! ADD THE COMPRESSION TERMS
                  ZE = ( -3.508914E-8*ZT-1.248266E-8 ) *ZT-2.595994E-6
                  ZBW= (  1.296821E-6*ZT-5.782165E-9 ) *ZT+1.045941E-4
                  ZB = ZBW + ZE * ZS

                  ZD = -2.042967E-2
                  ZC =   (-7.267926E-5*ZT+2.598241E-3 ) *ZT+0.1571896
                  ZAW= ( ( 5.939910E-6*ZT+2.512549E-3 ) *ZT-0.1028859 ) *ZT - 4.721788
                  ZA = ( ZD*ZSR + ZC ) *ZS + ZAW

                  ZB1=   (-0.1909078*ZT+7.390729 ) *ZT-55.87545
                  ZA1= ( ( 2.326469E-3*ZT+1.553190)*ZT-65.00517 ) *ZT+1044.077
                  ZKW= ( ( (-1.361629E-4*ZT-1.852732E-2 ) *ZT-30.41638 ) *ZT + 2098.925 ) *ZT+190925.6
                  ZK0= ( ZB1*ZSR + ZA1 )*ZS + ZKW

                  ZRDC1 = 1.0 / ( ZK0 - ZH * ( ZA - ZH * ZB ) )
                  ZRDC2 = 1.0 / ( 1.0 - ZH * ZRDC1 )
                  ! ============
                  ! ADJOINT PART
                  ! ============

                  ! MASKED IN SITU DENSITY

                  ZRHOPAD = ZRHOPAD + PRD_AD(JI,JJ,JK) * GRD%MSK(JI,JJ,JK)    &
                     &                                 * ZRDC2 
                  ZK0AD   = ZK0AD   - PRD_AD(JI,JJ,JK) * GRD%MSK(JI,JJ,JK)    &
                     &                                 * ZRDC2 * ZRDC2 * ZH &
                     &                                 * ZRDC1**2 * ZRHOP 
                  ZAAD    = ZAAD    + PRD_AD(JI,JJ,JK) * GRD%MSK(JI,JJ,JK)    &
                     &                                 * ZRDC2 * ZRDC2 * ZH &
                     &                                 * ZRDC1**2 * ZRHOP   &
                     &                                 * ZH
                  ZBAD    = ZBAD    - PRD_AD(JI,JJ,JK) * GRD%MSK(JI,JJ,JK)    &
                     &                                 * ZRDC2 * ZRDC2 * ZH &
                     &                                 * ZRDC1**2 * ZRHOP   &
                     &                                 * ZH * ZH
                  PRD_AD(JI,JJ,JK) = 0.0_R8

                  ZKWAD = ZKWAD + ZK0AD
                  ZSRAD = ZSRAD + ZK0AD * ZB1 * ZS
                  ZB1AD = ZB1AD + ZK0AD * ZS  * ZSR
                  ZA1AD = ZA1AD + ZK0AD * ZS
                  ZSAD  = ZSAD  + ZK0AD * ( ZB1 * ZSR + ZA1 )
                  ZK0AD = 0.0_R8

                  ZTAD  = ZTAD + ZKWAD * ( ( (-4.*1.361629E-4   * ZT &
                     &                        -3.*1.852732E-2 ) * ZT &
                     &                        -2.*30.41638    ) * ZT &
                     &                        +   2098.925           )
                  ZKWAD = 0.0_R8

                  ZTAD  = ZTAD + ZA1AD * ( ( 3.*2.326469E-3   * ZT &
                     &                      +2.*1.553190    ) * ZT &
                     &                      -   65.00517           )
                  ZA1AD = 0.0_R8

                  ZTAD  = ZTAD + ZB1AD * (-2.*0.1909078     * ZT &
                     &                    +   7.390729           )
                  ZB1AD = 0.0_R8

                  ZAWAD = ZAWAD + ZAAD
                  ZSRAD = ZSRAD + ZAAD *   ZD * ZS
                  ZCAD  = ZCAD  + ZAAD *   ZS
                  ZSAD  = ZSAD  + ZAAD * ( ZD * ZSR + ZC )
                  ZAAD  = 0.0_R8

                  ZTAD  = ZTAD + ZAWAD * ( ( 3.*5.939910E-6   * ZT &
                     &                      +2.*2.512549E-3 ) * ZT &
                     &                      -   0.1028859          )
                  ZAWAD = 0.0_R8

                  ZTAD  = ZTAD + ZCAD * (-2.*7.267926E-5   * ZT &
                     &                   +   2.598241E-3        )
                  ZCAD  = 0.0_R8

                  ZBWAD = ZBWAD + ZBAD
                  ZEAD  = ZEAD  + ZBAD * ZS
                  ZSAD  = ZSAD  + ZBAD * ZE
                  ZBAD  = 0.0_R8

                  ZTAD  = ZTAD + ZBWAD *  ( 2.*1.296821E-6   * ZT &
                     &                     -   5.782165E-9        )
                  ZBWAD = 0.0_R8

                  ZTAD  = ZTAD + ZEAD * (-2.*3.508914E-8   * ZT &
                     &                   -   1.248266E-8        )
                  ZEAD =  0.0_R8

                  ZR1AD   = ZR1AD + ZRHOPAD
                  ZR2AD   = ZR2AD + ZRHOPAD * ZS
                  ZR3AD   = ZR3AD + ZRHOPAD * ZSR * ZS
                  ZSRAD   = ZSRAD + ZRHOPAD * ZR3 * ZS
                  ZSAD    = ZSAD  + ZRHOPAD * ( 2. * ZR4 * ZS + ZR2  &
                     &                        + ZR3 * ZSR           )
                  ZRHOPAD = 0.0_R8

                  ZTAD  = ZTAD + ZR3AD * (-2.*1.6546E-6   * ZT &
                     &                    +   1.0227E-4        )
                  ZR3AD = 0.0_R8

                  ZTAD  = ZTAD + ZR2AD * ( ( ( 4.*5.3875E-9   * ZT &
                     &                        -3.*8.2467E-7 ) * ZT &
                     &                        +2.*7.6438E-5 ) * ZT &
                     &                        -   4.0899E-3        )
                  ZR2AD = 0.0_R8

                  ZTAD  = ZTAD + ZR1AD * ( ( ( ( 5.*6.536332E-9   * ZT &
                     &                          -4.*1.120083E-6 ) * ZT &
                     &                          +3.*1.001685E-4 ) * ZT &
                     &                          -2.*9.095290E-3 ) * ZT &
                     &                          +   6.793952E-2        )
                  ZR1AD = 0.0_R8

                  ZSAD  = ZSAD + ZSRAD * ( 1.0 / MAX( 2.*ZSR, ZEPS ) ) &
                     &                 * GRD%MSK(JI,JJ,JK)
                  ZSRAD = 0.0_R8

                  PSAL_AD(JI,JJ,JK) = PSAL_AD(JI,JJ,JK) + ZSAD
                  PTEM_AD(JI,JJ,JK) = PTEM_AD(JI,JJ,JK) + ZTAD
                  ZTAD = 0.0_R8
                  ZSAD = 0.0_R8
               END DO

            END DO
            !                                             ! ===============
         END DO                                           !   END OF SLAB
         !                                                ! ===============

END SUBROUTINE EOS_INSITU_AD
END MODULE EOSINSITU
