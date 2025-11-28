MODULE ICE_TRANSF

USE SET_KND
USE OBS_STR
USE GRD_STR
USE NETCDF
USE MYNETCDF
USE RUN, ONLY : IANMM
USE MPIREL, ONLY : CMPIDOM, MYPROC
USE IOUNITS , ONLY : IOUNERR,IOUNLOG
USE READFG

IMPLICIT NONE
REAL(R8), ALLOCATABLE :: SIC_quant(:,:,:,:) ! NO NEED TO EXTEND IF USED ONLY AT THE BEGINNING
REAL(R8), ALLOCATABLE :: SIT_quant(:,:,:,:) ! NO NEED TO EXTEND IF USED ONLY AT THE BEGINNING
REAL(R8), ALLOCATABLE :: SIC_COEFF_TLAD(:,:)
REAL(R8), ALLOCATABLE :: SIT_COEFF_TLAD(:,:)

REAL(R8), ALLOCATABLE :: qua_def_sic(:),qua_def_sit(:)
REAL(R8), ALLOCATABLE :: qua_ref_sic(:),qua_ref_sit(:)
INTEGER(I4) :: qdim_sic,qdim_sit
CHARACTER(LEN=300) :: CFILE_QUANTILE_SIC
CHARACTER(LEN=300) :: CFILE_QUANTILE_SIT
CONTAINS

SUBROUTINE ice_transf_ini

INTEGER(I4) :: MM
REAL(R8), ALLOCATABLE :: temp(:,:,:,:,:)

CFILE_QUANTILE_SIC="SIC_quant.nc"
CFILE_QUANTILE_SIT="SIT_quant.nc"
!month
WRITE(IOUNLOG,*) ' MONTH CHOSEN FOR ICE TRANSFORMATION'

MM=IANMM

!IF(OBS%SIC.GE.1) THEN
        WRITE(IOUNLOG,*) ' READING SIC TRANSFORMATION FROM' ,TRIM(CFILE_QUANTILE_SIC)

        CALL GETNCDIM(TRIM(CFILE_QUANTILE_SIC),"quantile",qdim_sic)
        ALLOCATE(qua_def_sic(qdim_sic),qua_ref_sic(qdim_sic))
        ALLOCATE(temp(1,qdim_sic,GRD%IM,GRD%JM,12))
        ALLOCATE(SIC_quant(1,qdim_sic,GRD%IM,GRD%JM))
        CALL READ_ANAMORF_FILE(TRIM(CFILE_QUANTILE_SIC),1,qdim_sic,GRD%IM,GRD%JM,12,temp,"ANAMORPH_Q",qua_def_sic,qua_ref_sic)
        SIC_quant(:,:,:,:)=temp(:,:,:,:,MM)
        DEALLOCATE(temp)
!ENDIF

!IF(OBS%SIT.GE.1) THEN
        WRITE(IOUNLOG,*) ' READING SIT TRANSFORMATION FROM ',TRIM(CFILE_QUANTILE_SIT)
        CALL GETNCDIM(TRIM(CFILE_QUANTILE_SIT),"quantile",qdim_sit)
        ALLOCATE(qua_def_sit(qdim_sit),qua_ref_sit(qdim_sit))
        ALLOCATE(temp(1,qdim_sit,GRD%IM,GRD%JM,12))
        ALLOCATE(SIT_quant(1,qdim_sit,GRD%IM,GRD%JM))
        CALL READ_ANAMORF_FILE(TRIM(CFILE_QUANTILE_SIT),1,qdim_sit,GRD%IM,GRD%JM,12,temp,"ANAMORPH_Q",qua_def_sit,qua_ref_sit)
        SIT_quant(:,:,:,:)=temp(:,:,:,:,MM)
!ENDIF
        DEALLOCATE(temp)


END SUBROUTINE ice_transf_ini

SUBROUTINE ICE_INI_COEFF(SIC_OR_SIT)
IMPLICIT NONE
INTEGER(I4), INTENT(IN) :: SIC_OR_SIT

IF (SIC_OR_SIT .EQ. 1) CALL ICE_INI_COEFFTL(SIC_OR_SIT)
IF (SIC_OR_SIT .EQ. 2) CALL ICE_INI_COEFFTL(SIC_OR_SIT)

WRITE(IOUNLOG,*) 'COEFFICIENT FOR TANGENT LINEAR CALCULATED!!'

END SUBROUTINE ICE_INI_COEFF


SUBROUTINE ICE_INI_COEFFTL(SIC_OR_SIT)

USE RUN, ONLY : NTSTEPS, LL_RESTART_FILE, NSUBDOMAINS

IMPLICIT NONE


INTEGER(I4), INTENT(IN) :: SIC_OR_SIT
INTEGER(I4) :: I,J,status
REAL(R8) :: X0,X1,Y0,Y1

 IF (SIC_OR_SIT .EQ. 1) THEN

 ALLOCATE(SIC_COEFF_TLAD(GRD%IM,GRD%JM))
 SIC_COEFF_TLAD=0._R8

     DO J=1,GRD%JM
     DO I=1,GRD%IM
!WRITE(IOUNLOG,*) I,J, GRD%MSK(I,J,1),NFGFILES,GRD%SICB(I,J,NFGFILES)
!        IF (GRD%MSK(I,J,1) .GT. 0.9_R8 .AND. GRD%SICB(I,J,NFGFILES) .GT. 0._R8 .AND. GRD%SICB(I,J,NFGFILES) .LT. 999999._R8 ) THEN
         IF (GRD%MSK(I,J,1) .GT. 0.9_R8 .AND. GRD%SICB(I,J,NFGFILES) .LT. 999999._R8 ) THEN

                !IF (GRD%SICB(I,J,NFGFILES) .LT. 0.95_R8 .AND. (GRD%SICB(I,J,NFGFILES) .GT. 0.05_R8)) THEN
                     X0=GRD%TRA_SICB(I,J,NFGFILES)-0.05
                   Y0=X0
                     CALL sangoma_anamorphosis(1,qdim_sic,-1,qua_ref_sic,SIC_QUANT(:,:,I,J),Y0,status)
                     X1=GRD%TRA_SICB(I,J,NFGFILES)+0.05
                     Y1=X1
                     CALL sangoma_anamorphosis(1,qdim_sic,-1,qua_ref_sic,SIC_QUANT(:,:,I,J),Y1,status)
                     SIC_COEFF_TLAD(I,J)=(Y1-Y0)/(X1-X0)

               !#ELSEIF (GRD%SICB(I,J,NFGFILES) .GE. 0.9_R8) THEN
               !      X0=GRD%SICB(I,J,NFGFILES)-0.1
               !      Y0=X0
               !      CALL sangoma_anamorphosis(1,qdim_sic,1,qua_ref_sic,SIC_QUANT(:,:,I,J),Y0,status)
               !      X1=GRD%SICB(I,J,NFGFILES)
               !      Y1=X1
               !      CALL sangoma_anamorphosis(1,qdim_sic,1,qua_ref_sic,SIC_QUANT(:,:,I,J),Y1,status)
               !      !WRITE(IOUNLOG,*) I,J,Y1,Y0,X1,X0
               !      SIC_COEFF_TLAD(I,J)=(Y1-Y0)/(X1-X0)
                 
               ! ELSE
               !      X0=GRD%SICB(I,J,NFGFILES)
               !      Y0=X0
               !      CALL sangoma_anamorphosis(1,qdim_sic,1,qua_ref_sic,SIC_QUANT(:,:,I,J),Y0,status)
               !      X1=GRD%SICB(I,J,NFGFILES)+0.1
               !      Y1=X1
               !      CALL sangoma_anamorphosis(1,qdim_sic,1,qua_ref_sic,SIC_QUANT(:,:,I,J),Y1,status)
               !      SIC_COEFF_TLAD(I,J)=(Y1-Y0)/(X1-X0)
               ! ENDIF
        ENDIF
      ENDDO
      ENDDO

CALL WRITE_VARNCDF('TLAD_'//TRIM(CMPIDOM)//'.NC',&
      & GRD%IM,GRD%JM,SIC_COEFF_TLAD,'TLAD',&
      & ' ','sic tlad')


 ENDIF



IF (SIC_OR_SIT .EQ. 2) THEN

ALLOCATE(SIT_COEFF_TLAD(GRD%IM,GRD%JM))
 SIT_COEFF_TLAD=0._R8

     DO J=1,GRD%JM
     DO I=1,GRD%IM
        IF (GRD%MSK(I,J,1) .GT. 0.9_R8 .AND. GRD%SITB(I,J,NFGFILES) .LT. 999999._R8 ) THEN
                !IF (GRD%SITB(I,J,NFGFILES) .LT. 0.99_R8 .AND. (GRD%SITB(I,J,NFGFILES) .GT. 0.01_R8)) THEN
                     X0=GRD%TRA_SITB(I,J,NFGFILES)-0.05
                     Y0=X0
                     CALL sangoma_anamorphosis(1,qdim_sit,-1,qua_ref_sit,SIT_QUANT(:,:,I,J),Y0,status)
                     X1=GRD%TRA_SITB(I,J,NFGFILES)+0.05
                     Y1=X1
                     CALL sangoma_anamorphosis(1,qdim_sit,-1,qua_ref_sit,SIT_QUANT(:,:,I,J),Y1,status)
                     SIT_COEFF_TLAD(I,J)=(Y1-Y0)/(X1-X0) !*10
               ! ELSEIF (GRD%SITB(I,J,NFGFILES) .GE. 0.99_R8) THEN
               !      X0=GRD%SITB(I,J,NFGFILES)-0.02
               !      Y0=X0
               !      CALL sangoma_anamorphosis(1,qdim_sit,1,qua_ref_sit,SIT_QUANT(:,:,I,J),Y0,status)
               !      X1=GRD%SITB(I,J,NFGFILES)
               !      Y1=X1
               !      CALL sangoma_anamorphosis(1,qdim_sit,1,qua_ref_sit,SIT_QUANT(:,:,I,J),Y1,status)
               !      SIT_COEFF_TLAD(I,J)=(Y1-Y0)/(X1-X0)
               ! ELSE
               !      X0=GRD%SITB(I,J,NFGFILES)
               !      Y0=X0
               !      CALL sangoma_anamorphosis(1,qdim_sit,1,qua_ref_sit,SIT_QUANT(:,:,I,J),Y0,status)
               !      X1=GRD%SITB(I,J,NFGFILES)+0.02
               !      Y1=X1
               !      CALL sangoma_anamorphosis(1,qdim_sit,1,qua_ref_sit,SIT_QUANT(:,:,I,J),Y1,status)
               !      SIT_COEFF_TLAD(I,J)=(Y1-Y0)/(X1-X0)
               ! ENDIF
        ENDIF
      ENDDO
      ENDDO

CALL WRITE_VARNCDF('TLAD_SIT_'//TRIM(CMPIDOM)//'.NC',&
      & GRD%IM,GRD%JM,SIT_COEFF_TLAD,'TLAD',&
      & ' ','sit tlad')

  ENDIF



END SUBROUTINE ICE_INI_COEFFTL

SUBROUTINE sisort(n, veca)

! Sorts a vector veca(1:n) into ascending numerical order, by
! straight insertion.
!
! For large vectors, this routine will not be efficient.

  IMPLICIT NONE

  INTEGER, INTENT(in) :: n
  REAL, INTENT(inout) :: veca(n)

  INTEGER :: i, j, k
  REAL :: tmpa
  LOGICAL :: eflag

  DO j = 2, n

     eflag = .false.

     tmpa = veca(j)

     sortloop: DO i = j-1, 1, -1
        k = i

        IF(veca(i) <= tmpa) THEN
           eflag = .true.
           EXIT sortloop
        END IF

        veca(i+1) = veca(i)
     ENDDO sortloop
     IF (.NOT.eflag) k=0

     veca(k+1) = tmpa

  ENDDO

END SUBROUTINE sisort


SUBROUTINE sisort2(n, veca, vecb)

! Sorts a vector veca(1:n) into ascending numerical order, by
! straight insertion. vecb is sorted in the same order.
!
! For large vectors, this routine will not be efficient.

  IMPLICIT NONE

  INTEGER, INTENT(in) :: n
  REAL, INTENT(inout) :: veca(n), vecb(n)

  INTEGER :: i, j, k
  REAL :: tmpa, tmpb
  LOGICAL :: eflag

  DO j = 2, n

     tmpa = veca(j)
     tmpb = vecb(j)

     sortloop: DO i = j-1, 1, -1
        k = i

        IF(veca(i) <= tmpa) THEN
           eflag = .true.
           EXIT sortloop
        END IF

        veca(i+1) = veca(i)
        vecb(i+1) = vecb(i)
     ENDDO sortloop
     IF (.NOT.eflag) k=0

     veca(k+1) = tmpa
     vecb(k+1) = tmpb

  ENDDO

END SUBROUTINE sisort2


END MODULE ICE_TRANSF
