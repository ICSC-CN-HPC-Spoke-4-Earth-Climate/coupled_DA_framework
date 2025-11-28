PROGRAM MAIN

USE SET_KND
USE MPI
USE PNETCDF
USE SORTING

IMPLICIT NONE

!include <mpif.h>

INTEGER(MPI_OFFSET_KIND) ::  NX, NY, NZ, MSD, NMR,NXTEMP,NYTEMP
INTEGER(MPI_OFFSET_KIND) :: DX_rec,DY_rec,DX_send,DY_send,DZ
INTEGER(MPI_OFFSET_KIND) :: DX_send_diag,DX_send_dx,DX_send_dy
INTEGER(MPI_OFFSET_KIND) :: DY_send_diag,DY_send_dx,DY_send_dy
INTEGER(MPI_OFFSET_KIND), ALLOCATABLE, DIMENSION(:) :: X_ST, X_EN, Y_ST, Y_EN, NDIR, XCST, XCEN
INTEGER(MPI_OFFSET_KIND), ALLOCATABLE, DIMENSION(:) :: X2_ST, X2_EN, Y2_ST, Y2_EN, ND1, ND2, YMST, YMEN
INTEGER :: info,STAT

INTEGER(I4) :: JI, JJ, NY2, ND,var_xtype
REAL(R4), ALLOCATABLE, DIMENSION(:,:,:)  :: SA_TE,VERYTEMP
REAL(R4), ALLOCATABLE, DIMENSION(:,:,:,:)  :: TE_SA_overlap_rec,TE_SA_overlap_send
REAL(R4), ALLOCATABLE, DIMENSION(:,:,:,:)  :: TE_SA_overlap_send_diag,TE_SA_overlap_send_dx,TE_SA_overlap_send_dy
!REAL(R8), ALLOCATABLE :: LON(:,:), LAT(:,:), DEP(:)
REAL(R4) :: ZZ, ZPI
CHARACTER(LEN=99) :: CFA

INTEGER(I4) :: X_DIMID,Y_DIMID,Z_DIMID,DIMIDS3(3),IERR=0,FLUSH_IND=0
INTEGER(I4) :: J_VAR_ATT,VAR_ATT_number,status(MPI_status_size)
INTEGER(MPI_OFFSET_KIND) :: count_rec(3),count_send(3),start_in(3),start_out(3)
INTEGER(MPI_OFFSET_KIND) :: number_rec,number_rec_2D,number_send,number_send_2D,start_in_send(3),start_in_rec(3)
INTEGER(MPI_OFFSET_KIND) :: number_send_diag,number_send_diag_2D,start_in_send_diag(3),count_send_diag(3)
INTEGER(MPI_OFFSET_KIND) :: number_send_dx,number_send_dx_2D,start_in_send_dx(3),count_send_dx(3)
INTEGER(MPI_OFFSET_KIND) :: number_send_dy,number_send_dy_2D,start_in_send_dy(3),count_send_dy(3)


INTEGER(I4) :: NCID,NCID_sec_file
INTEGER(I4) :: VARID,VARID_first_file_SA,VARID_first_file_TE,VARID_first_file_SSH
INTEGER(I4) :: NCUOTSAL,NCUOTTEM,NCUOTSSH,VARIDOUT_SA,VARIDOUT_TE,VARIDOUT_SSH
INTEGER(I8) :: number
CHARACTER(LEN=99) :: C_NAME_VAR_ATT

REAL(R4) :: COEF_REC1,COEF_REC2,COEF_REC3,COEF_REC4,REL_POSX,REL_POSY

!CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: processor_name
INTEGER(I4) :: rank, size

integer :: new_rank,colour,size_new_comm,new_group,world_group,new_comm
integer :: new_rank2,colour2,size_new_comm2,new_group2,world_group2,new_comm2

INTEGER(I4), ALLOCATABLE, DIMENSION(:)  :: MP_SD,MP_SDST
INTEGER(I4), ALLOCATABLE, DIMENSION(:,:)  :: MPI_MAP
INTEGER(I4) :: DOMAINS,NXDOM,NYDOM,JX,JY,COUNTX,COUNTY,NUM_PROCS_ONMAP,ROWX,COLMY
INTEGER(I4), ALLOCATABLE, DIMENSION(:)  :: PROC_XMED,PROC_YMED,PROC_XMED_ORDERED,PROC_YMED_ORDERED,PROC_DIAG_SW,PROC_DIAG_NE
INTEGER(I4), ALLOCATABLE, DIMENSION(:)  :: PROC_XINI,PROC_YINI,PROC_XEND,PROC_YEND,PROC_XINI_ORDERED,PROC_YINI_ORDERED,PROC_XEND_ORDERED,PROC_YEND_ORDERED
INTEGER(I4), ALLOCATABLE, DIMENSION(:)  :: PROC_WE,PROC_EA,PROC_NO,PROC_SO,PROC_POSX_IN_MAP,PROC_POSY_IN_MAP
INTEGER(I4), ALLOCATABLE, DIMENSION(:)  :: MAP_XMED,MAP_YMED,NPROC_TEMP_SORT
INTEGER(I4), ALLOCATABLE, DIMENSION(:)  :: MAP_XINI,MAP_YINI,MAP_XEND,MAP_YEND
LOGICAL, ALLOCATABLE, DIMENSION(:)  :: WE_BOUNDARY,EA_BOUNDARY,NO_BOUNDARY,SO_BOUNDARY
LOGICAL :: EXIST_SSH

call MPI_INIT ( ierr )
call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, size, ierr)
!call MPI_Get_processor_name(processor_name,namelen,ierr)
call MPI_COMM_group(MPI_COMM_WORLD,world_group,ierr)


CALL MPI_INFO_CREATE(info, ierr)
!call MPI_Info_set (info, "nc_header_read_chunk_size", "8388608",ierr);
call MPI_Info_set(info, "cb_buffer_size", "16777216",ierr)


!---define grid---------!


OPEN(UNIT=99,FILE='mpp_conf.dat',STATUS='OLD')
READ(99,*) DOMAINS

!PLEASE READ FILE IN THIS WAY FOR XCUT DOMAIN
!READ(99,*) DOMAINS,NXDOM,NYDOM
!GOFS CONDITION
NXDOM=1
NYDOM=DOMAINS


ALLOCATE( MP_SD(DOMAINS), MP_SDST(DOMAINS), X_ST(DOMAINS), X_EN(DOMAINS), &
& Y_ST(DOMAINS), Y_EN(DOMAINS) )

DO JJ=1,DOMAINS
    READ(99,*) MP_SD(JJ),MP_SDST(JJ),X_ST(JJ),X_EN(JJ),Y_ST(JJ),Y_EN(JJ)
    ENDDO
CLOSE(99)

NY=MAXVAL(Y_EN)
NX=MAXVAL(X_EN)

ALLOCATE(PROC_XMED(DOMAINS),PROC_YMED(DOMAINS))
ALLOCATE(PROC_XINI(DOMAINS),PROC_YINI(DOMAINS),PROC_XEND(DOMAINS),PROC_YEND(DOMAINS))
ALLOCATE(PROC_XMED_ORDERED(DOMAINS),PROC_YMED_ORDERED(DOMAINS))
ALLOCATE(PROC_XINI_ORDERED(DOMAINS),PROC_YINI_ORDERED(DOMAINS))
ALLOCATE(PROC_XEND_ORDERED(DOMAINS),PROC_YEND_ORDERED(DOMAINS))
ALLOCATE(PROC_WE(DOMAINS),PROC_EA(DOMAINS),PROC_NO(DOMAINS),PROC_SO(DOMAINS))
ALLOCATE(PROC_POSX_IN_MAP(DOMAINS),PROC_POSY_IN_MAP(DOMAINS))
ALLOCATE(MAP_XMED(NXDOM),MAP_YMED(NYDOM),MAP_XINI(NXDOM),MAP_YINI(NYDOM),MAP_XEND(NXDOM),MAP_YEND(NYDOM))
ALLOCATE(MPI_MAP(NXDOM,NYDOM))
ALLOCATE(NPROC_TEMP_SORT(DOMAINS),PROC_DIAG_SW(DOMAINS),PROC_DIAG_NE(DOMAINS))
ALLOCATE(WE_BOUNDARY(DOMAINS),EA_BOUNDARY(DOMAINS),NO_BOUNDARY(DOMAINS),SO_BOUNDARY(DOMAINS))


PROC_XMED(:)= INT((X_ST(:)+ X_EN(:))/2)
PROC_YMED(:)= INT((Y_ST(:)+ Y_EN(:))/2)
PROC_XINI(:)= INT(X_ST(:))
PROC_YINI(:)= INT(Y_ST(:))
PROC_XEND(:)= INT(X_EN(:))
PROC_YEND(:)= INT(Y_EN(:))


NPROC_TEMP_SORT(:)=REAL(PROC_XMED(:),R4)
CALL SORT(NPROC_TEMP_SORT,DOMAINS)
PROC_XMED_ORDERED(:)=INT(NPROC_TEMP_SORT(:),I4)

NPROC_TEMP_SORT(:)=REAL(PROC_XINI(:),R4)
CALL SORT(NPROC_TEMP_SORT,DOMAINS)
PROC_XINI_ORDERED(:)=INT(NPROC_TEMP_SORT(:),I4)

NPROC_TEMP_SORT(:)=REAL(PROC_XEND(:),R4)
CALL SORT(NPROC_TEMP_SORT,DOMAINS)
PROC_XEND_ORDERED(:)=INT(NPROC_TEMP_SORT(:),I4)

NPROC_TEMP_SORT(:)=REAL(PROC_YMED(:),R4)
CALL SORT(NPROC_TEMP_SORT,DOMAINS)
PROC_YMED_ORDERED(:)=INT(NPROC_TEMP_SORT(:),I4)

NPROC_TEMP_SORT(:)=REAL(PROC_YINI(:),R4)
CALL SORT(NPROC_TEMP_SORT,DOMAINS)
PROC_YINI_ORDERED(:)=INT(NPROC_TEMP_SORT(:),I4)

NPROC_TEMP_SORT(:)=REAL(PROC_YEND(:),R4)
CALL SORT(NPROC_TEMP_SORT,DOMAINS)
PROC_YEND_ORDERED(:)=INT(NPROC_TEMP_SORT(:),I4)

DEALLOCATE(NPROC_TEMP_SORT)

MAP_XMED(1)=PROC_XMED_ORDERED(1)
MAP_YMED(1)=PROC_YMED_ORDERED(1)

MAP_XINI(1)=PROC_XINI_ORDERED(1)
MAP_YINI(1)=PROC_YINI_ORDERED(1)

MAP_XEND(1)=PROC_XEND_ORDERED(1)
MAP_YEND(1)=PROC_YEND_ORDERED(1)


COUNTX=1
COUNTY=1
DO JJ=2,DOMAINS
	IF (PROC_XMED_ORDERED(JJ-1) .NE. PROC_XMED_ORDERED(JJ)) THEN
		COUNTX=COUNTX+1
		IF (COUNTX .GT. NXDOM) THEN
            WRITE(*,*)    "DOMAINS READ > NXDOM !! please check mpp_conf.dat "
            STOP
        ENDIF
        MAP_XMED(COUNTX)=PROC_XMED_ORDERED(JJ)
        MAP_XINI(COUNTX)=PROC_XINI_ORDERED(JJ)
        MAP_XEND(COUNTX)=PROC_XEND_ORDERED(JJ)
	ENDIF
	IF (PROC_YMED_ORDERED(JJ-1) .NE. PROC_YMED_ORDERED(JJ)) THEN
		COUNTY=COUNTY+1
		IF (COUNTY .GT. NYDOM) THEN
            WRITE(*,*) "DOMAINS READ > NYDOM !! please check mpp_conf.dat "
            STOP
        ENDIF
        MAP_YMED(COUNTY)=PROC_YMED_ORDERED(JJ)
        MAP_YINI(COUNTY)=PROC_YINI_ORDERED(JJ)
        MAP_YEND(COUNTY)=PROC_YEND_ORDERED(JJ)
	ENDIF
ENDDO

IF((COUNTX .NE. NXDOM) .OR. (COUNTY .NE. NYDOM)) THEN
   WRITE(*,*) "DOMAINS READ < NXDOM,NYDOM !!! please check mpp_conf.dat"
    STOP
ENDIF


PROC_POSX_IN_MAP(:)=-1
PROC_POSY_IN_MAP(:)=-1
MPI_MAP(:,:)=-1
	
DO JJ=1,DOMAINS
	DO JX=1,NXDOM
		if (PROC_XMED(JJ) .EQ. MAP_XMED(JX)) THEN
		PROC_POSX_IN_MAP(JJ)=JX
		CONTINUE
		ENDIF
	ENDDO
	IF(PROC_POSX_IN_MAP(JJ) .EQ. -1) THEN
        WRITE(*,*) "PROC=JJ-1->",JJ-1," MISMATCH X DIVISIONS AND DOMAINS!!! please check mpp_conf.dat"
        STOP
	ENDIF
	DO JY=1,NYDOM
		IF (PROC_YMED(JJ) .EQ. MAP_YMED(JY)) THEN
		PROC_POSY_IN_MAP(JJ)=JY
		CONTINUE
		ENDIF
	ENDDO
	IF(PROC_POSY_IN_MAP(JJ) .EQ. -1) THEN
            WRITE(*,*) "PROC=JJ-1->",JJ-1," MISMATCH Y DIVISION AND DOMAINS!!! please check mpp_conf.dat"
        STOP
    ENDIF

	MPI_MAP(PROC_POSX_IN_MAP(JJ),PROC_POSY_IN_MAP(JJ))=JJ
ENDDO


PROC_WE(:)=-1
PROC_EA(:)=-1
PROC_NO(:)=-1
PROC_SO(:)=-1
PROC_DIAG_SW(:)=-1
PROC_DIAG_NE(:)=-1

DO JJ=1,DOMAINS
	!..NEXT EA PROC (INCREASE X)
	IF(PROC_POSX_IN_MAP(JJ).LT.NXDOM) &
		PROC_EA(JJ)=MPI_MAP(PROC_POSX_IN_MAP(JJ)+1,PROC_POSY_IN_MAP(JJ))
	!..NEXT WE PROC (DECREASE X)
	IF(PROC_POSX_IN_MAP(JJ).GT.1) &
		PROC_WE(JJ)=MPI_MAP(PROC_POSX_IN_MAP(JJ)-1,PROC_POSY_IN_MAP(JJ))
    !..NEXT NO PROC (INCREASE Y) !NO BOUNDARY CONDITION
	IF(PROC_POSY_IN_MAP(JJ).LT.NYDOM) &
	&	PROC_NO(JJ)=MPI_MAP(PROC_POSX_IN_MAP(JJ),PROC_POSY_IN_MAP(JJ)+1)
	!..NEXT SO PROC (DECREASE Y)
	IF(PROC_POSY_IN_MAP(JJ).GT.1) &
	&	PROC_SO(JJ)=MPI_MAP(PROC_POSX_IN_MAP(JJ),PROC_POSY_IN_MAP(JJ)-1)
    !..NEXT DIAG SOUTH-WEST PROC (DECREASE X and Y)
	IF(PROC_POSY_IN_MAP(JJ).GT.1 .AND. PROC_POSX_IN_MAP(JJ).GT.1) &
    & PROC_DIAG_SW(JJ)=MPI_MAP(PROC_POSX_IN_MAP(JJ)-1,PROC_POSY_IN_MAP(JJ)-1)
    !..NEXT DIAG NORTH-EAST PROC (INCREASE X and Y)
	IF(PROC_POSX_IN_MAP(JJ).LT.NXDOM .AND. PROC_POSY_IN_MAP(JJ).LT.NYDOM) &
    & PROC_DIAG_NE(JJ)=MPI_MAP(PROC_POSX_IN_MAP(JJ)+1,PROC_POSY_IN_MAP(JJ)+1)

    !! THEY CAN STILL BE -1 !!
ENDDO

WE_BOUNDARY(:)=.FALSE.
EA_BOUNDARY(:)=.FALSE.
NO_BOUNDARY(:)=.FALSE.
SO_BOUNDARY(:)=.FALSE.

DO JJ=1,DOMAINS
    IF (PROC_POSX_IN_MAP(JJ) .EQ. 1_I4 ) WE_BOUNDARY(JJ)=.TRUE.
    IF (PROC_POSX_IN_MAP(JJ) .EQ. NXDOM) EA_BOUNDARY(JJ)=.TRUE.
    IF (PROC_POSY_IN_MAP(JJ) .EQ. 1_I4 ) SO_BOUNDARY(JJ)=.TRUE.
    IF (PROC_POSY_IN_MAP(JJ) .EQ. NYDOM) NO_BOUNDARY(JJ)=.TRUE.
ENDDO



IF( rank .EQ. 0) THEN
    WRITE(*,*) " LAYOUT OF PROCESSORS AS READ, LAND is -1"
    WRITE(*,*) ""
    11210 FORMAT(A)
    11211 FORMAT(A,I3,A)
    11212 FORMAT(A,I4,A,I4,A)
    DO JY=NYDOM,1,-1
        DO JX=1,NXDOM
            WRITE(*,11210,ADVANCE="no") " -----------------"
        ENDDO
        WRITE(*,*) ""
        DO JX=1,NXDOM
                WRITE(*,11211,ADVANCE="no") "       ",MPI_MAP(JX,JY),"       |"
        ENDDO
        WRITE(*,*) ""
        DO JX=1,NXDOM
            WRITE(*,11212,ADVANCE="no") " (x=",MAP_XMED(JX),",y=",MAP_YMED(JY),") |"
       !     WRITE(*,11212,ADVANCE="no") " (x=",MAP_XINI(JX),",y=",MAP_YINI(JY),") |"
       !     WRITE(*,11212,ADVANCE="no") " (x=",MAP_XEND(JX),",y=",MAP_YEND(JY),") |"
        ENDDO
        WRITE(*,*) ""
        IF(MPI_MAP(1,JY) .EQ. -1_I4) THEN
            WRITE(*,11210,ADVANCE="no") "                 |"
        ELSE IF ( PROC_WE(MPI_MAP(1,JY)) .NE. -1_I4) THEN
            WRITE(*,11210,ADVANCE="no") "<< WE            |"
        ELSE 
		WRITE(*,11210,ADVANCE="no") "                 |"
        ENDIF
        DO JX=2,NXDOM-1
            WRITE(*,11210,ADVANCE="no") "                 |"
        ENDDO
        IF(MPI_MAP(NXDOM,JY) .EQ. -1_I4) THEN
            WRITE(*,11210,ADVANCE="no") "                 |"
        ELSE IF ( PROC_EA(MPI_MAP(NXDOM,JY)) .NE. -1_I4) THEN
            WRITE(*,11210,ADVANCE="no") "            EA >>|"
        ELSE
            WRITE(*,11210,ADVANCE="no") "                 |"
        ENDIF

        WRITE(*,*) ""
    ENDDO
    DO JX=1,NXDOM
        WRITE(*,11210,ADVANCE="no") " -----------------"
    ENDDO
    WRITE(*,*) ""
    WRITE(*,*) ""

ENDIF

NUM_PROCS_ONMAP=0
DO JX=1,NXDOM
DO JY=1,NYDOM
    IF(MPI_MAP(JX,JY) .NE. -1) NUM_PROCS_ONMAP=NUM_PROCS_ONMAP+1
ENDDO
ENDDO

IF(NUM_PROCS_ONMAP .NE. DOMAINS) STOP  "SOME DOMAINS ARE REPLICATED!!! please check mpp_conf.dat"

DO JJ=1,DOMAINS
	IF (MPI_MAP(PROC_POSX_IN_MAP(JJ),PROC_POSY_IN_MAP(JJ)) .EQ. -1) THEN
    WRITE(*,*)"ERROR IN MAPPING THE DOMAIN NUM =",JJ
    STOP
    ENDIF
ENDDO

DEALLOCATE(PROC_XMED_ORDERED,PROC_YMED_ORDERED,PROC_XINI_ORDERED,PROC_YINI_ORDERED,PROC_XEND_ORDERED,PROC_YEND_ORDERED)!,PROC_POSX_IN_MAP,PROC_POSY_IN_MAP)
DEALLOCATE(MAP_XMED,MAP_YMED)

!---- domain read---!


!OPEN(21,FILE='merge_domains.dat',STATUS='OLD')
!READ(21,*) DOMAINS, NX, NY, NZ

WRITE(*,*) DOMAINS,size
if (DOMAINS /= size) stop "ERROR!! In this version, # of procs must be equal to # domains (T,S,SSH) !"


!ALLOCATE( X_ST(DOMAINS), X_EN(DOMAINS), Y_ST(DOMAINS), Y_EN(DOMAINS) )
ALLOCATE( X2_ST(DOMAINS), X2_EN(DOMAINS), Y2_ST(DOMAINS), Y2_EN(DOMAINS) )
ALLOCATE( YMST(DOMAINS), YMEN(DOMAINS), XCST(DOMAINS), XCEN(DOMAINS) )

!DO JJ=1,DOMAINS
!  READ(21,*) X_ST(JJ), X_EN(JJ), Y_ST(JJ), Y_EN(JJ)
!ENDDO

DO JJ=1,DOMAINS
    IF( .NOT. WE_BOUNDARY(JJ)) THEN
        IF(PROC_WE(JJ) .NE. -1) THEN
            X2_ST(JJ)=X_EN(PROC_WE(JJ))
        ELSE
            ROWX=PROC_POSX_IN_MAP(JJ)-1
            DO  JY=1,NYDOM
                IF(MPI_MAP(ROWX,JY) .NE. -1) THEN
                    X2_ST(JJ)=X_EN(MPI_MAP(ROWX,JY))
                    CONTINUE
                ENDIF
            ENDDO
        ENDIF
    ELSE
        X2_ST(JJ)=X_ST(JJ)
    ENDIF
    IF( .NOT. EA_BOUNDARY(JJ)  ) THEN
        IF(PROC_EA(JJ) .NE. -1) THEN
            X2_EN(JJ)=X_ST(PROC_EA(JJ))
        ELSE
            ROWX=PROC_POSX_IN_MAP(JJ)+1
            DO  JY=1,NYDOM
                IF(MPI_MAP(ROWX,JY) .NE. -1) THEN
                    X2_EN(JJ)=X_ST(MPI_MAP(ROWX,JY))
                    CONTINUE
                ENDIF
            ENDDO
        ENDIF
    ELSE
        X2_EN(JJ)=X_EN(JJ)
    ENDIF
    IF( .NOT. SO_BOUNDARY(JJ) ) THEN
        IF(PROC_SO(JJ) .NE. -1 ) THEN
            Y2_ST(JJ)=Y_EN(PROC_SO(JJ))
        ELSE
            COLMY=PROC_POSY_IN_MAP(JJ)-1
            DO  JX=1,NXDOM
                IF(MPI_MAP(JX,COLMY) .NE. -1) THEN
                    Y2_ST(JJ)=Y_EN(MPI_MAP(JX,COLMY))
                    CONTINUE
                ENDIF
            ENDDO
        ENDIF
    ELSE
        Y2_ST(JJ)=Y_ST(JJ)
    ENDIF
    IF( .NOT. NO_BOUNDARY(JJ) ) THEN
        IF(PROC_NO(JJ) .NE. -1 ) THEN
            Y2_EN(JJ)=Y_ST(PROC_NO(JJ))
        ELSE
            COLMY=PROC_POSY_IN_MAP(JJ)+1
            DO  JX=1,NXDOM
                IF(MPI_MAP(JX,COLMY) .NE. -1) THEN
                    Y2_EN(JJ)=Y_ST(MPI_MAP(JX,COLMY))
                    CONTINUE
                ENDIF
            ENDDO
        ENDIF
    ELSE
        Y2_EN(JJ)=Y_EN(JJ)
    ENDIF
ENDDO


DO JJ=1,DOMAINS
    IF( .NOT. EA_BOUNDARY(JJ)) THEN
        XCST(JJ)=X2_EN(JJ)
        XCEN(JJ)=X_EN(JJ)
    ELSE
        XCST(JJ)=X_EN(JJ)
        XCEN(JJ)=X_EN(JJ)
    ENDIF
    IF( .NOT. NO_BOUNDARY(JJ) ) THEN
        YMST(JJ)=Y2_EN(JJ)
        YMEN(JJ)=Y_EN(JJ)
    ELSE
        YMST(JJ)=Y_EN(JJ)
        YMEN(JJ)=Y_EN(JJ)
    ENDIF
ENDDO


! open one file to define variables
if (rank==0) THEN
!WRITE(*,*)WE_BOUNDARY
!WRITE(*,*)EA_BOUNDARY
!WRITE(*,*)SO_BOUNDARY
!WRITE(*,*)NO_BOUNDARY
ENDIF

if (rank==0) print*, "----Reconstructiong one by one S , T , ssh----"
FLUSH_IND=FLUSH_IND+1
CALL FLUSH(FLUSH_IND)

!---rules for colours- first splitting----
!colour=0
!if (rank >= INT(size/3) .AND. rank < INT(2*size/3) ) colour=1
!if (rank >= INT(2*size/3) ) colour=2
!--- splitting

DO colour=0,2

!call MPI_COMM_SPLIT(MPI_COMM_WORLD,colour,rank,new_comm,ierr)
!call MPI_COMM_group(new_comm,new_group,ierr)
!call MPI_COMM_SIZE(new_comm,size_new_comm,ierr)
!call MPI_COMM_RANK(new_comm,new_rank,ierr)

size_new_comm=size
new_rank=rank
new_comm=MPI_COMM_WORLD

write(*,*) "new_size",size_new_comm," for color ",colour," old rank=", rank,", local rank=",new_rank
FLUSH_IND=FLUSH_IND+1
CALL FLUSH(FLUSH_IND)



WRITE(CFA,'(A)') 'ANINCR.NC.0000'
CALL CHECK( NFmpi_OPEN(MPI_COMM_WORLD,CFA, NF_NOWRITE,info, NCID),__LINE__ )
CALL CHECK( NFmpi_INQ_DIMID(NCID, 'z',Z_DIMID),__LINE__  )
CALL CHECK( NFmpi_INQ_DIMLEN(NCID, Z_DIMID,NZ),__LINE__  )

EXIST_SSH=.TRUE.
STAT=NFmpi_INQ_VARID(NCID, 'INCSSHEIG', VARID)
IF (STAT /= NF_NOERR) THEN
write(*,*) "SSH does not seem to be present"
EXIST_SSH=.FALSE.
ENDIF

if (0==colour ) THEN
            !-----SALINITY-----
    CALL CHECK( NFMPI_CREATE(new_comm,"ANINCR_S_MERGED.nc", NF_64BIT_OFFSET, info, NCUOTSAL),__LINE__  )
    CALL CHECK( NFmpi_DEF_DIM(NCUOTSAL, 'x',NX,  X_DIMID),__LINE__  )
    CALL CHECK( NFmpi_DEF_DIM(NCUOTSAL, 'y',NY, Y_DIMID),__LINE__  )
    CALL CHECK( NFmpi_DEF_DIM(NCUOTSAL, 'z',NZ,  Z_DIMID),__LINE__  )
    CALL CHECK( NFmpi_INQ_VARID(NCID, 'INCSALINE', VARID),__LINE__  )
    CALL  CHECK( nfmpi_inq_vartype(NCID, VARID,var_xtype),__LINE__  )
    DIMIDS3=  (/X_DIMID, Y_DIMID, Z_DIMID /)
    CALL CHECK( NFmpi_DEF_VAR(NCUOTSAL, 'INCSALINE',var_xtype,3, DIMIDS3, VARIDOUT_SA),__LINE__  )
    CALL CHECK( NFmpi_INQ_varnatts(NCID,VARID,VAR_ATT_number),__LINE__  )
    DO J_VAR_ATT=1,VAR_ATT_number
        CALL CHECK(  NFmpi_inq_attname(NCID,VARID,J_VAR_ATT,C_NAME_VAR_ATT),__LINE__  )
        CALL CHECK(  NFmpi_copy_att(NCID,VARID,C_NAME_VAR_ATT,NCUOTSAL,VARIDOUT_SA),__LINE__  )
    END DO
    call CHECK( NFmpi_enddef(NCUOTSAL),__LINE__  )

ELSEIF (1==colour ) THEN
            !------TEMPERATURE-----
    CALL CHECK( NFmpi_CREATE(new_comm,"ANINCR_T_MERGED.nc",NF_64BIT_OFFSET,info,NCUOTTEM),__LINE__  )
    CALL CHECK( NFmpi_DEF_DIM(NCUOTTEM, 'x',NX,  X_DIMID),__LINE__  )
    CALL CHECK( NFmpi_DEF_DIM(NCUOTTEM, 'y',NY, Y_DIMID),__LINE__  )
    CALL CHECK( NFmpi_DEF_DIM(NCUOTTEM, 'z',NZ,  Z_DIMID),__LINE__  )
    CALL CHECK( NFmpi_INQ_VARID(NCID, 'INCTEMPER', VARID),__LINE__  )
    CALL  CHECK( nfmpi_inq_vartype(NCID, VARID,var_xtype),__LINE__  )
    DIMIDS3=  (/X_DIMID, Y_DIMID, Z_DIMID /) !pay attention!! DIMIDS3 is the same for SALINE and TEMP since have been just defined
    CALL CHECK( NFmpi_DEF_VAR(NCUOTTEM, 'INCTEMPER',var_xtype, 3,DIMIDS3, VARIDOUT_TE),__LINE__  )
    CALL CHECK( NFmpi_INQ_varnatts(NCID,VARID, VAR_ATT_number),__LINE__  )
    DO J_VAR_ATT=1,VAR_ATT_number
        CALL CHECK(  NFmpi_inq_attname(NCID,VARID,J_VAR_ATT,C_NAME_VAR_ATT),__LINE__  )
        CALL CHECK(  NFmpi_copy_att(NCID,VARID,C_NAME_VAR_ATT,NCUOTTEM,VARIDOUT_TE),__LINE__  )
    END DO
    call CHECK( NFmpi_enddef(NCUOTTEM),__LINE__  )

ELSEIF (2==colour .AND. EXIST_SSH ) THEN
            !------TEMPERATURE-----
    CALL CHECK( NFmpi_CREATE(new_comm,"ANINCR_SSH_MERGED.nc",NF_64BIT_OFFSET,info,NCUOTSSH),__LINE__  )
    CALL CHECK( NFmpi_DEF_DIM(NCUOTSSH, 'x',NX,  X_DIMID),__LINE__  )
    CALL CHECK( NFmpi_DEF_DIM(NCUOTSSH, 'y',NY, Y_DIMID),__LINE__  )
    CALL CHECK( NFmpi_INQ_VARID(NCID, 'INCSSHEIG', VARID),__LINE__  )
    CALL  CHECK( nfmpi_inq_vartype(NCID, VARID,var_xtype),__LINE__  )
    DIMIDS3=  (/X_DIMID, Y_DIMID,Z_DIMID /) !pay attention!! DIMIDS3 is the same for SALINE and TEMP since have been just defined
CALL CHECK( NFmpi_DEF_VAR(NCUOTSSH, 'INCSSHEIG',var_xtype, 2,DIMIDS3(1:2), VARIDOUT_SSH),__LINE__  )
    CALL CHECK( NFmpi_INQ_varnatts(NCID,VARID, VAR_ATT_number),__LINE__  )
    DO J_VAR_ATT=1,VAR_ATT_number
        CALL CHECK(  NFmpi_inq_attname(NCID,VARID,J_VAR_ATT,C_NAME_VAR_ATT),__LINE__  )
        CALL CHECK(  NFmpi_copy_att(NCID,VARID,C_NAME_VAR_ATT,NCUOTSSH,VARIDOUT_SSH),__LINE__  )
    END DO
    call CHECK( NFmpi_enddef(NCUOTSSH),__LINE__  )


ENDIF

CALL CHECK( NFmpi_CLOSE(NCID),__LINE__  )! closing file for definitions

!---------
    !----------
    JI=new_rank+1   !---using both JI and rank, rank for IF/ELSE CONDITION, JI for extracting values from arrays
    !----------
    !----------

    !opening files to be read
    WRITE(CFA,'(A,I4.4)') 'ANINCR.NC.',JI-1
    WRITE(*,*) ' OPENING FILE ', TRIM(CFA)
    CALL CHECK( NFmpi_OPEN(MPI_COMM_SELF,TRIM(CFA), NF_NOWRITE,info,NCID),__LINE__  )
    !-----

!-----COPYING non-overlapping regions
    if (0==rank)  print *, "---------- merging NON-OVERLAPING DATA "


    DX_rec=X2_EN(JI)-X2_ST(JI)+1
    DY_rec=Y2_EN(JI)-Y2_ST(JI)+1
    DZ=NZ


    ALLOCATE( SA_TE(X_EN(JI)-X2_ST(JI)+1,Y_EN(JI)-Y2_ST(JI)+1,DZ) ) !USING THE SAME ARRAY FOR SALINITY OR TEMPERATURE
    SA_TE(:,:,:)=0._R4

    start_in_rec=(/X2_ST(JI)-X_ST(JI)+1_MPI_OFFSET_KIND,Y2_ST(JI)-Y_ST(JI)+1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/) !array defining the first point to be read
    start_out=(/X2_ST(JI),Y2_ST(JI),1_MPI_OFFSET_KIND/)  !array defining the first point to be written
    count_rec=(/DX_rec,DY_rec,DZ/)

!if (0==rank) THEN
!print *,rank, "YMST ",YMST
!print *,rank, "YMEM ",YMEN
!print *,rank, "Y_ST ",Y_ST
!print *,rank, "Y_EN ",Y_EN
!print *,rank, "Y2_ST ",Y2_ST
!print *,rank, "Y2_EN ",Y2_EN
!print *,rank, "XCST ",XCST
!print *,rank, "XCEN",XCEN
!print *,rank, "XST ",X_ST
!print *,rank, "X_EN ",X_EN
!print *,rank, "X2_ST ",X2_ST
!print *,rank, "X2_EN ",X2_EN
!
!ENDIF


if (0 == colour ) THEN
    !---salitnity----
CALL CHECK( NFmpi_INQ_VARID(NCID, 'INCSALINE', VARID_first_file_SA),__LINE__  )
CALL CHECK(NFmpi_GET_VARA_real_all(NCID,VARID_first_file_SA,start_in_rec,count_rec,SA_TE(1:DX_rec,1:DY_rec,:)),__LINE__  )
!    CALL CHECK( nfmpi_put_vara_real_all(NCUOTSAL,VARIDOUT_SA,start_out,count_rec,SA_TE),__LINE__  )
ELSE IF (1 == colour ) THEN
    !---temperature---
CALL CHECK( NFmpi_INQ_VARID(NCID, 'INCTEMPER', VARID_first_file_TE),__LINE__  )
CALL CHECK(  NFmpi_GET_VARA_real_all(NCID,VARID_first_file_TE,start_in_rec,count_rec,SA_TE(1:DX_rec,1:DY_rec,:)),__LINE__  )
!    CALL CHECK(  nfmpi_put_vara_real_all(NCUOTTEM,VARIDOUT_TE,start_out,count_rec,SA_TE),__LINE__  )
ELSE IF (2 == colour .AND. EXIST_SSH ) THEN
!---temperature---
CALL CHECK( NFmpi_INQ_VARID(NCID, 'INCSSHEIG', VARID_first_file_SSH),__LINE__  )
CALL CHECK(  NFmpi_GET_VARA_real_all(NCID,VARID_first_file_SSH,start_in_rec(1:2),count_rec(1:2),SA_TE(1:DX_rec,1:DY_rec,1)),__LINE__  )
!    CALL CHECK(  nfmpi_put_vara_real_all(NCUOTTEM,VARIDOUT_TE,start_out,count_rec,SA_TE),__LINE__  )
ENDIF
!    DEALLOCATE(SA_TE)
!
!!
!-----COPYING SINGLE DX OVERLAPPING REGION
    if (0==rank) WRITE(*,*) '---------- NOW merging SINGLE OVERLAPPING DX DATA'
    !DZ=NZ

    !--ALMOST every process receive and send data with different size
    !receive data
    DX_rec=X_EN(JI)-X2_EN(JI)+1
    DY_rec=Y2_EN(JI)-Y2_ST(JI)+1
    start_in_rec=(/XCST(JI)-X_ST(JI)+1_MPI_OFFSET_KIND,Y2_ST(JI)-Y_ST(JI)+1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    !start_in_rec=(/X2_ST(JI)-X_ST(JI)+1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    count_rec=(/DX_rec,DY_rec,DZ/)
    number_rec=DX_rec*DY_rec*DZ
    number_rec_2D=DX_rec*DY_rec
    !send data
    DX_send=X2_ST(JI)-X_ST(JI)+1
    DY_send=Y2_EN(JI)-Y2_ST(JI)+1
    start_in_send=(/1_MPI_OFFSET_KIND,Y2_ST(JI)-Y_ST(JI)+1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    !start_in_send=(/1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    count_send=(/DX_send,DY_send,DZ/)
    number_send=DX_send*DY_send*DZ
    number_send_2D=DX_send*DY_send


    !---each process allocate both rec and send even if process 0 and NDA would never use both
    ALLOCATE( TE_SA_overlap_rec(DX_rec,DY_rec,DZ,2),stat=IERR )
    ALLOCATE( TE_SA_overlap_send(DX_send,DY_send,DZ,1),stat=IERR )
    TE_SA_overlap_rec(:,:,:,:)=0._R4
    TE_SA_overlap_send(:,:,:,:)=0._R4
    if (IERR/=0) print *," not allocated rank",rank

    !the following routine works well since comunicatori is MPI_COMM_SELF
    !if it were MPI_COMM_WORLD,  it would have failed (since get_var is not for all processes)

    ! read data to be iterpolated with data from rank+1
    if ( .NOT. EA_BOUNDARY(JI) ) THEN
if (0==colour) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SA,start_in_rec,count_rec, TE_SA_overlap_rec(:,:,:,1)),__LINE__  )
if (1==colour) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_TE,start_in_rec,count_rec, TE_SA_overlap_rec(:,:,:,1)),__LINE__  )
if (2==colour .AND. EXIST_SSH ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SSH,start_in_rec(1:2),count_rec(1:2), TE_SA_overlap_rec(:,:,1,1)),__LINE__  )
    ENDIF
    ! read data to be sent to rank-1
    if ( .NOT. WE_BOUNDARY(JI) ) THEN
if (0==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SA,start_in_send,count_send,TE_SA_overlap_send(:,:,:,1)),__LINE__  )
if (1==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_TE,start_in_send,count_send,TE_SA_overlap_send(:,:,:,1)),__LINE__  )
if (2==colour .AND. EXIST_SSH) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SSH,start_in_send(1:2),count_send(1:2),TE_SA_overlap_send(:,:,1,1)),__LINE__  )
    ENDIF


if ( .NOT. WE_BOUNDARY(JI) .AND. PROC_WE(JI) .NE. -1 ) THEN
    if (0==colour) call MPI_Send(TE_SA_overlap_send(:,:,:,1), number_send, MPI_REAL4,PROC_WE(JI)-1,11,new_comm,IERR)
    if (1==colour) call MPI_Send(TE_SA_overlap_send(:,:,:,1), number_send, MPI_REAL4,PROC_WE(JI)-1,12,new_comm,IERR)
    if (2==colour .AND. EXIST_SSH ) call MPI_Send(TE_SA_overlap_send(:,:,1,1), number_send_2D, MPI_REAL4,PROC_WE(JI)-1,13,new_comm,IERR)
    ENDIF

if ( .NOT.EA_BOUNDARY(JI) .AND. PROC_EA(JI) .NE. -1 ) THEN
    if (0==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,2), number_rec, MPI_REAL4,PROC_EA(JI)-1,11,new_comm,status,IERR)
    if (1==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,2), number_rec, MPI_REAL4,PROC_EA(JI)-1,12,new_comm,status,IERR)
    if (2==colour .AND. EXIST_SSH) call MPI_Recv(TE_SA_overlap_rec(:,:,1,2), number_rec_2D, MPI_REAL4,PROC_EA(JI)-1,13,new_comm,status,IERR)
ENDIF



!----INTERPOLATE RESULTS IN OVERLAPPING AREA---
!---interpolation does not cover the last rank

    if (0==new_rank) print*," INTERPOLATE RESULTS IN OVERLAPPING AREA..."


! for SSH also the rest of the layers are interpolated although being zero
if ( .NOT.EA_BOUNDARY(JI) .AND. PROC_EA(JI) .NE. -1 ) THEN

           ZPI = 2._R4 * ASIN( 1._R4 ) ! pi=3.1415..
            DO JJ=1,DX_rec
                ZZ = 0.5_R4 * (  1._R4 - COS( ZPI * REAL( JJ-1 ) / REAL(DX_rec-1) )  )
                ZZ = 1._R4 - ZZ
                TE_SA_overlap_rec(JJ,:,:,2) =  ZZ*TE_SA_overlap_rec(JJ,:,:,1) + (1._R4-ZZ)*TE_SA_overlap_rec(JJ,:,:,2)
            ENDDO
            SA_TE((X2_EN(JI)-X2_ST(JI)+1):(X_EN(JI)-X2_ST(JI)+1),1:(Y2_EN(JI)-Y2_ST(JI)+1),:)=TE_SA_overlap_rec(:,:,:,2)
ENDIF

!

DEALLOCATE( TE_SA_overlap_rec,TE_SA_overlap_send )

call MPI_BARRIER(MPI_COMM_WORLD,IERR)



!-----COPYING SINGLE DY OVERLAPPING REGION
    if (0==rank) WRITE(*,*) '---------- NOW merging OVERLAPPING DY DATA'
    DZ=NZ

    !--ALMOST every process receive and send data with different size
    !receive data
    DX_rec=X2_EN(JI)-X2_ST(JI)+1
    DY_rec=Y_EN(JI)-Y2_EN(JI)+1
    start_in_rec=(/X2_ST(JI)-X_ST(JI)+1_MPI_OFFSET_KIND,Y2_EN(JI)-Y_ST(JI)+1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    count_rec=(/DX_rec,DY_rec,DZ/)
    number_rec=DX_rec*DY_rec*DZ
    number_rec_2D=DX_rec*DY_rec
    print *," rank",JI,start_in_rec,count_rec
!send data
    DX_send=X2_EN(JI)-X2_ST(JI)+1
    DY_send=Y2_ST(JI)-Y_ST(JI)+1
    start_in_send=(/X2_ST(JI)-X_ST(JI)+1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    count_send=(/DX_send,DY_send,DZ/)
    number_send=DX_send*DY_send*DZ
    number_send_2D=DX_send*DY_send


    !---each process allocate both rec and send even if process 0 and NDA would never use both
    ALLOCATE( TE_SA_overlap_rec(DX_rec,DY_rec,DZ,2),stat=IERR )
    ALLOCATE( TE_SA_overlap_send(DX_send,DY_send,DZ,1),stat=IERR )
    TE_SA_overlap_rec(:,:,:,:)=0._R4
    TE_SA_overlap_send(:,:,:,:)=0._R4
    if (IERR/=0) print *," not allocated rank",rank

    !the following routine works well since comunicatori is MPI_COMM_SELF
    !if it were MPI_COMM_WORLD,  it would have failed (since get_var is not for all processes)

    ! read data to be iterpolated with data from rank+1
    if (.NOT. NO_BOUNDARY(JI) ) THEN
if (0==colour) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SA,start_in_rec,count_rec, TE_SA_overlap_rec(:,:,:,1)),__LINE__  )
if (1==colour) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_TE,start_in_rec,count_rec, TE_SA_overlap_rec(:,:,:,1)),__LINE__  )
if (2==colour .AND. EXIST_SSH) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SSH,start_in_rec(1:2),count_rec(1:2), TE_SA_overlap_rec(:,:,1,1)),__LINE__  )
    ENDIF
    ! read data to be sent to rank-1
    if (.NOT. SO_BOUNDARY(JI)) THEN
if (0==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SA,start_in_send,count_send,TE_SA_overlap_send(:,:,:,1)),__LINE__  )
if (1==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_TE,start_in_send,count_send,TE_SA_overlap_send(:,:,:,1)),__LINE__  )
if (2==colour .AND. EXIST_SSH) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SSH,start_in_send(1:2),count_send(1:2),TE_SA_overlap_send(:,:,1,1)),__LINE__  )
    ENDIF



if ( .NOT. SO_BOUNDARY(JI) .AND. PROC_SO(JI) .NE. -1 ) THEN
    if (0==colour) call MPI_Send(TE_SA_overlap_send(:,:,:,1), number_send, MPI_REAL4,PROC_SO(JI)-1,13,new_comm,IERR)
    if (1==colour) call MPI_Send(TE_SA_overlap_send(:,:,:,1), number_send, MPI_REAL4,PROC_SO(JI)-1,14,new_comm,IERR)
    if (2==colour .AND. EXIST_SSH) call MPI_Send(TE_SA_overlap_send(:,:,1,1), number_send_2D, MPI_REAL4,PROC_SO(JI)-1,15,new_comm,IERR)
ENDIF

if ( .NOT. NO_BOUNDARY(JI) .AND. PROC_NO(JI) .NE. -1) THEN
    if (0==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,2), number_rec, MPI_REAL4,PROC_NO(JI)-1,13,new_comm,status,IERR)
    if (1==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,2), number_rec, MPI_REAL4,PROC_NO(JI)-1,14,new_comm,status,IERR)
    if (2==colour .AND. EXIST_SSH) call MPI_Recv(TE_SA_overlap_rec(:,:,1,2), number_rec_2D, MPI_REAL4,PROC_NO(JI)-1,15,new_comm,status,IERR)
ENDIF



!----INTERPOLATE RESULTS IN OVERLAPPING AREA---
!---interpolation does not cover the last rank

    if (0==new_rank) print*," INTERPOLATE RESULTS IN OVERLAPPING AREA..."

if (.NOT. NO_BOUNDARY(JI) .AND. PROC_NO(JI) .NE. -1) THEN

        ZPI = 2._R4 * ASIN( 1._R4 ) ! pi=3.1415..

        DO JJ=1,DY_rec

            ZZ = 0.5_R4 * (  1._R4 - COS( ZPI * REAL( JJ-1 ) / REAL(DY_rec-1) )  )
            ZZ = 1._R4 - ZZ

                TE_SA_overlap_rec(:,JJ,:,2) =  ZZ*TE_SA_overlap_rec(:,JJ,:,1) + (1._R4-ZZ)*TE_SA_overlap_rec(:,JJ,:,2)
        ENDDO

        SA_TE(1:(X2_EN(JI)-X2_ST(JI)+1),(Y2_EN(JI)-Y2_ST(JI)+1):(Y_EN(JI)-Y2_ST(JI)+1),:)=TE_SA_overlap_rec(:,:,:,2)

ENDIF

!

DEALLOCATE( TE_SA_overlap_rec,TE_SA_overlap_send )


call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!
!
!
!
!-----COPYING DOUBLE DX DY OVERLAPPING REGION
    if (0==rank) WRITE(*,*) '---------- NOW merging OVERLAPPING DX-DY DATA'
    DZ=NZ

    !--ALMOST every process receive and send data with different size
    !receive data
    DX_rec=X_EN(JI)-X2_EN(JI)+1
    DY_rec=Y_EN(JI)-Y2_EN(JI)+1
    start_in_rec=(/X2_EN(JI)-X_ST(JI)+1_MPI_OFFSET_KIND,Y2_EN(JI)-Y_ST(JI)+1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    count_rec=(/DX_rec,DY_rec,DZ/)
    number_rec=DX_rec*DY_rec*DZ
    number_rec_2D=DX_rec*DY_rec
    !send data diagonal
    DX_send_diag=X2_ST(JI)-X_ST(JI)+1
    DY_send_diag=Y2_ST(JI)-Y_ST(JI)+1
    start_in_send_diag=(/1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    count_send_diag=(/DX_send_diag,DY_send_diag,DZ/)
    number_send_diag=DX_send_diag*DY_send_diag*DZ
    number_send_diag_2D=DX_send_diag*DY_send_diag
    !send data dx
    DX_send_dx=X2_ST(JI)-X_ST(JI)+1
    DY_send_dx=Y_EN(JI)-Y2_EN(JI)+1
    start_in_send_dx=(/1_MPI_OFFSET_KIND,Y2_EN(JI)-Y_ST(JI)+1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    count_send_dx=(/DX_send_dx,DY_send_dx,DZ/)
    number_send_dx=DX_send_dx*DY_send_dx*DZ
    number_send_dx_2D=DX_send_dx*DY_send_dx
    !send data dy
    DX_send_dy=X_EN(JI)-X2_EN(JI)+1
    DY_send_dy=Y2_ST(JI)-Y_ST(JI)+1
    start_in_send_dy=(/X2_EN(JI)-X_ST(JI)+1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND,1_MPI_OFFSET_KIND/)
    count_send_dy=(/DX_send_dy,DY_send_dy,DZ/)
    number_send_dy=DX_send_dy*DY_send_dy*DZ
    number_send_dy_2D=DX_send_dy*DY_send_dy

    !---each process allocate both rec and send even if process 0 and NDA would never use both
    ALLOCATE( TE_SA_overlap_rec(DX_rec,DY_rec,DZ,4),stat=IERR )
    ALLOCATE( TE_SA_overlap_send_diag(DX_send_diag,DY_send_diag,DZ,1),&
    & TE_SA_overlap_send_dx(DX_send_dx,DY_send_dx,DZ,1),TE_SA_overlap_send_dy(DX_send_dy,DY_send_dy,DZ,1),stat=IERR )
    TE_SA_overlap_rec(:,:,:,:)=0._R4
    TE_SA_overlap_send_diag(:,:,:,:)=0._R4
    TE_SA_overlap_send_dx(:,:,:,:)=0._R4
    TE_SA_overlap_send_dy(:,:,:,:)=0._R4
    if (IERR/=0) print *," not allocated rank",rank

    !the following routine works well since comunicatori is MPI_COMM_SELF
    !if it were MPI_COMM_WORLD,  it would have failed (since get_var is not for all processes)

    !
    if ( .NOT. EA_BOUNDARY(JI) .AND. .NOT. NO_BOUNDARY(JI) ) THEN
        if (0==colour) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SA,start_in_rec,count_rec, TE_SA_overlap_rec(:,:,:,1)),__LINE__  )
        if (1==colour) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_TE,start_in_rec,count_rec, TE_SA_overlap_rec(:,:,:,1)),__LINE__  )
        if (2==colour .AND. EXIST_SSH) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SSH,start_in_rec(1:2),count_rec(1:2), TE_SA_overlap_rec(:,:,1,1)),__LINE__  )
    ENDIF
    !
    if ( .NOT. WE_BOUNDARY(JI) .AND. .NOT. SO_BOUNDARY(JI) ) THEN
        if (0==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SA,start_in_send_diag,count_send_diag,TE_SA_overlap_send_diag(:,:,:,1)),__LINE__  )
        if (1==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_TE,start_in_send_diag,count_send_diag,TE_SA_overlap_send_diag(:,:,:,1)),__LINE__  )
        if (2==colour .AND. EXIST_SSH) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SSH,start_in_send_diag(1:2),count_send_diag(1:2),TE_SA_overlap_send_diag(:,:,1,1)),__LINE__  )
    ENDIF

    !
    if ( .NOT. WE_BOUNDARY(JI) ) THEN
        if (0==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SA,start_in_send_dx,count_send_dx,TE_SA_overlap_send_dx(:,:,:,1)),__LINE__  )
        if (1==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_TE,start_in_send_dx,count_send_dx,TE_SA_overlap_send_dx(:,:,:,1)),__LINE__  )
        if (2==colour .AND. EXIST_SSH) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SSH,start_in_send_dx(1:2),count_send_dx(1:2),TE_SA_overlap_send_dx(:,:,1,1)),__LINE__  )
    ENDIF
	
    !
    if ( .NOT. SO_BOUNDARY(JI) ) THEN
        if (0==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SA,start_in_send_dy,count_send_dy,TE_SA_overlap_send_dy(:,:,:,1)),__LINE__  )
        if (1==colour ) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_TE,start_in_send_dy,count_send_dy,TE_SA_overlap_send_dy(:,:,:,1)),__LINE__  )
        if (2==colour .AND. EXIST_SSH) CALL CHECK(NFmpi_GET_VARA_REAL_all(NCID,VARID_first_file_SSH,start_in_send_dy(1:2),count_send_dy(1:2),TE_SA_overlap_send_dy(:,:,1,1)),__LINE__  )
    ENDIF

    if (.NOT. WE_BOUNDARY(JI) .AND. .NOT. SO_BOUNDARY(JI) .AND. (PROC_DIAG_SW(JI)) .NE. -1 ) THEN
    if (0==colour) call MPI_Send(TE_SA_overlap_send_diag(:,:,:,1), number_send_diag, MPI_REAL4,PROC_DIAG_SW(JI)-1,111,new_comm,IERR)
    if (1==colour) call MPI_Send(TE_SA_overlap_send_diag(:,:,:,1), number_send_diag, MPI_REAL4,PROC_DIAG_SW(JI)-1,112,new_comm,IERR)
    if (2==colour .AND. EXIST_SSH) call MPI_Send(TE_SA_overlap_send_diag(:,:,1,1), number_send_diag_2D, MPI_REAL4,PROC_DIAG_SW(JI)-1,113,new_comm,IERR)
    ENDIF

    if ( .NOT. EA_BOUNDARY(JI) .AND.  .NOT. NO_BOUNDARY(JI) .AND. (PROC_DIAG_NE(JI)) .NE. -1 ) THEN
        if (0==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,2), number_rec, MPI_REAL4,PROC_DIAG_NE(JI)-1,111,new_comm,status,IERR)
        if (1==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,2), number_rec, MPI_REAL4,PROC_DIAG_NE(JI)-1,112,new_comm,status,IERR)
        if (2==colour .AND. EXIST_SSH) call MPI_Recv(TE_SA_overlap_rec(:,:,1,2), number_rec_2D, MPI_REAL4,PROC_DIAG_NE(JI)-1,113,new_comm,status,IERR)
    ENDIF

call MPI_BARRIER(MPI_COMM_WORLD,IERR)

    if (.NOT. WE_BOUNDARY(JI) .AND. PROC_WE(JI) .NE. -1 ) THEN
        if (0==colour) call MPI_Send(TE_SA_overlap_send_dx(:,:,:,1), number_send_dx, MPI_REAL4,PROC_WE(JI)-1,114,new_comm,IERR)
        if (1==colour) call MPI_Send(TE_SA_overlap_send_dx(:,:,:,1), number_send_dx, MPI_REAL4,PROC_WE(JI)-1,115,new_comm,IERR)
        if (2==colour .AND. EXIST_SSH) call MPI_Send(TE_SA_overlap_send_dx(:,:,1,1), number_send_dx_2D, MPI_REAL4,PROC_WE(JI)-1,116,new_comm,IERR)
    ENDIF

    if (.NOT. EA_BOUNDARY(JI) .AND. PROC_EA(JI) .NE. -1  ) THEN
        if (0==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,3), number_rec, MPI_REAL4,PROC_EA(JI)-1,114,new_comm,status,IERR)
        if (1==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,3), number_rec, MPI_REAL4,PROC_EA(JI)-1,115,new_comm,status,IERR)
        if (2==colour .AND. EXIST_SSH) call MPI_Recv(TE_SA_overlap_rec(:,:,1,3), number_rec_2D, MPI_REAL4,PROC_EA(JI)-1,116,new_comm,status,IERR)
    ENDIF

call MPI_BARRIER(MPI_COMM_WORLD,IERR)


    if (.NOT. SO_BOUNDARY(JI) .AND. PROC_SO(JI) .NE. -1  ) THEN
        if (0==colour) call MPI_Send(TE_SA_overlap_send_dy(:,:,:,1), number_send_dy, MPI_REAL4,PROC_SO(JI)-1,117,new_comm,IERR)
        if (1==colour) call MPI_Send(TE_SA_overlap_send_dy(:,:,:,1), number_send_dy, MPI_REAL4,PROC_SO(JI)-1,118,new_comm,IERR)
        if (2==colour .AND. EXIST_SSH) call MPI_Send(TE_SA_overlap_send_dy(:,:,1,1), number_send_dy_2D, MPI_REAL4,PROC_SO(JI)-1,119,new_comm,IERR)
    ENDIF

    if ( .NOT. NO_BOUNDARY(JI) .AND. PROC_NO(JI) .NE. -1 ) THEN
        if (0==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,4), number_rec, MPI_REAL4,PROC_NO(JI)-1,117,new_comm,status,IERR)
        if (1==colour) call MPI_Recv(TE_SA_overlap_rec(:,:,:,4), number_rec, MPI_REAL4,PROC_NO(JI)-1,118,new_comm,status,IERR)
        if (2==colour .AND. EXIST_SSH) call MPI_Recv(TE_SA_overlap_rec(:,:,1,4), number_rec_2D, MPI_REAL4,PROC_NO(JI)-1,119,new_comm,status,IERR)
    ENDIF



!----INTERPOLATE RESULTS IN OVERLAPPING AREA---
!---interpolation does not cover the last rank

    if (0==new_rank) print*," INTERPOLATE RESULTS IN OVERLAPPING DX-DY AREA..."

if (.NOT. NO_BOUNDARY(JI) .AND. .NOT. EA_BOUNDARY(JI) .AND. PROC_NO(JI) .NE. -1 .AND. PROC_EA(JI) .NE. -1 .AND. PROC_DIAG_NE(JI) .NE. -1) THEN

        ZPI = 2._R4 * ASIN( 1._R4 ) ! pi=3.1415..

        DO JX=1,DX_rec
        DO JY=1,DY_rec

!LINEAR
!            COEF_REC1= (1._R4-REAL(JX-1)/REAL(DX_rec-1))*(1._R4-REAL(JY-1)/REAL(DY_rec-1))
!            COEF_REC2= (REAL(JX-1)/REAL(DX_rec-1))*(REAL(JX-1)/REAL(DX_rec-1))
!            COEF_REC3= (REAL(JX-1)/REAL(DX_rec-1))*(1._R4-REAL(JY-1)/REAL(DY_rec-1))
!            COEF_REC4= (1._R4-REAL(JX-1)/REAL(DX_rec-1))*(REAL(JY-1)/REAL(DY_rec-1))
!COSIN
             REL_POSY=REAL(JY-1)/REAL(DY_rec-1)
             REL_POSX=REAL(JX-1)/REAL(DX_rec-1)

             COEF_REC1=0.5_R4*(1._R4+COS(ZPI*REL_POSX))*0.5_R4*(1._R4+COS(ZPI*REL_POSY))
             COEF_REC2=0.5_R4*(1._R4-COS(ZPI*REL_POSX))*0.5_R4*(1._R4-COS(ZPI*REL_POSY))
             COEF_REC3=0.5_R4*(1._R4-COS(ZPI*REL_POSX))*0.5_R4*(1._R4+COS(ZPI*REL_POSY))
             COEF_REC4=0.5_R4*(1._R4+COS(ZPI*REL_POSX))*0.5_R4*(1._R4-COS(ZPI*REL_POSY))
!EQUAL TO
!             COEF_REC1=0.5_R4*(1._R4-COS(ZPI*(1._R4- REL_POSX)))*0.5_R4*(1._R4-COS(ZPI*(1._R4-REL_POSY)))
!             COEF_REC2=0.5_R4*(1._R4-COS(ZPI*REL_POSX))*0.5_R4*(1._R4-COS(ZPI*REL_POSY))
!             COEF_REC3=0.5_R4*(1._R4-COS(ZPI*REL_POSX))*0.5_R4*(1._R4-COS(ZPI*(1._R4-REL_POSY)))
!             COEF_REC4=0.5_R4*(1._R4-COS(ZPI*(1._R4- REL_POSX)))*0.5_R4*(1._R4-COS(ZPI*(REL_POSY)))

            TE_SA_overlap_rec(JX,JY,:,2)= COEF_REC1*TE_SA_overlap_rec(JX,JY,:,1)+COEF_REC2*TE_SA_overlap_rec(JX,JY,:,2)+ &
            & COEF_REC3*TE_SA_overlap_rec(JX,JY,:,3)+COEF_REC4*TE_SA_overlap_rec(JX,JY,:,4)

        ENDDO
        ENDDO

        SA_TE((X2_EN(JI)-X2_ST(JI)+1):(X_EN(JI)-X2_ST(JI)+1),(Y2_EN(JI)-Y2_ST(JI)+1):(Y_EN(JI)-Y2_ST(JI)+1),:)=TE_SA_overlap_rec(:,:,:,2)

ENDIF

!

DEALLOCATE( TE_SA_overlap_rec,TE_SA_overlap_send_diag,TE_SA_overlap_send_dx,TE_SA_overlap_send_dy )




CALL CHECK( NFmpi_CLOSE(NCID),__LINE__  )

! WRITE(*,*) ' OVERLAP ', JI, ' DOMAINS : ',ND1(JI),ND2(JI),' INDEXES : ',YMST(JI),YMEN(JI)

call MPI_BARRIER(MPI_COMM_WORLD,IERR)




if (0==new_rank) print*," ...done"


!write(*,*) "rank,",rank," ",count_rec," start ",start_out," size ",(X_EN(JI)-X2_ST(JI)+1)," ",&
!& Y_EN(JI)-Y2_ST(JI)+1," ",DZ

    if ( 0== colour) THEN
        start_out=(/X2_ST(JI),Y2_ST(JI),1_MPI_OFFSET_KIND/)
        count_rec=(/X_EN(JI)-X2_ST(JI)+1_MPI_OFFSET_KIND,Y_EN(JI)-Y2_ST(JI)+1_MPI_OFFSET_KIND,DZ/)
        CALL CHECK(nfmpi_put_vara_real_all(NCUOTSAL,VARIDOUT_SA,start_out,count_rec,SA_TE(:,:,:)),__LINE__  )

        DO JX=1,NXDOM
        DO JY=1,NYDOM
            IF ( MPI_MAP(JX,JY)== -1) THEN
                NXTEMP=MAP_XEND(JX)-MAP_XINI(JX)+1_MPI_OFFSET_KIND
                NYTEMP=MAP_YEND(JY)-MAP_YINI(JY)+1_MPI_OFFSET_KIND
                ALLOCATE(VERYTEMP(NXTEMP,NYTEMP,DZ))
                VERYTEMP(:,:,:)=0._R4
                if (0==rank) WRITE(*,*) 'FILLING with ZEROS,  x=(',MAP_XINI(JX),'',MAP_XEND(JX),') y=(',MAP_YINI(JY),'',MAP_YEND(JY),')'
                start_out=(/MAP_XINI(JX),MAP_YINI(JY),1/)
                count_rec=(/NXTEMP,NYTEMP,DZ/)
                CALL CHECK(nfmpi_put_vara_real_all(NCUOTSAL,VARIDOUT_SA,start_out,count_rec,VERYTEMP(:,:,:)),__LINE__  )
                DEALLOCATE(VERYTEMP)
            ENDIF
        ENDDO
        ENDDO
        CALL CHECK( NFmpi_close(NCUOTSAL),__LINE__  )
    ENDIF
    IF ( 1== colour) THEN
        start_out=(/X2_ST(JI),Y2_ST(JI),1_MPI_OFFSET_KIND/)
        count_rec=(/X_EN(JI)-X2_ST(JI)+1_MPI_OFFSET_KIND,Y_EN(JI)-Y2_ST(JI)+1_MPI_OFFSET_KIND,DZ/)
        CALL CHECK(nfmpi_put_vara_real_all(NCUOTTEM,VARIDOUT_TE,start_out,count_rec,SA_TE(:,:,:)),__LINE__  )
        DO JX=1,NXDOM
        DO JY=1,NYDOM
            IF ( MPI_MAP(JX,JY)== -1) THEN
                NXTEMP=MAP_XEND(JX)-MAP_XINI(JX)+1_MPI_OFFSET_KIND
                NYTEMP=MAP_YEND(JY)-MAP_YINI(JY)+1_MPI_OFFSET_KIND
                ALLOCATE(VERYTEMP(NXTEMP,NYTEMP,DZ))
                VERYTEMP(:,:,:)=0._R4
                if (0==rank) WRITE(*,*) 'FILLING with ZEROS, AREA  x=(',MAP_XINI(JX),'',MAP_XEND(JX),') y=(',MAP_YINI(JY),'',MAP_YEND(JY),')'
                start_out=(/MAP_XINI(JX),MAP_YINI(JY),1/)
                count_rec=(/NXTEMP,NYTEMP,DZ/)
                CALL CHECK(nfmpi_put_vara_real_all(NCUOTSAL,VARIDOUT_SA,start_out,count_rec,VERYTEMP(:,:,:)),__LINE__  )
                DEALLOCATE(VERYTEMP)
            ENDIF
        ENDDO
        ENDDO
        CALL CHECK( NFmpi_close(NCUOTTEM),__LINE__  )

    ENDIF
     IF ( 2== colour .AND. EXIST_SSH) THEN
        start_out=(/X2_ST(JI),Y2_ST(JI),1_MPI_OFFSET_KIND/)
        count_rec=(/X_EN(JI)-X2_ST(JI)+1_MPI_OFFSET_KIND,Y_EN(JI)-Y2_ST(JI)+1_MPI_OFFSET_KIND,DZ/)
        CALL CHECK(nfmpi_put_vara_real_all(NCUOTSSH,VARIDOUT_SSH,start_out(1:2),count_rec(1:2),SA_TE(:,:,1)),__LINE__  )
        DO JX=1,NXDOM
        DO JY=1,NYDOM
            IF ( MPI_MAP(JX,JY)== -1) THEN
                NXTEMP=MAP_XEND(JX)-MAP_XINI(JX)+1_MPI_OFFSET_KIND
                NYTEMP=MAP_YEND(JY)-MAP_YINI(JY)+1_MPI_OFFSET_KIND
                ALLOCATE(VERYTEMP(NXTEMP,NYTEMP,1))
                VERYTEMP(:,:,1)=0._R4
                if (0==rank) WRITE(*,*) 'FILLING with ZEROS,  x=(',MAP_XINI(JX),'',MAP_XEND(JX),') y=(',MAP_YINI(JY),'',MAP_YEND(JY),')'
                start_out=(/MAP_XINI(JX),MAP_YINI(JY),1/)
                count_rec=(/NXTEMP,NYTEMP,DZ/)
                CALL CHECK(nfmpi_put_vara_real_all(NCUOTSAL,VARIDOUT_SA,start_out(1:2),count_rec(1:2),VERYTEMP(:,:,1)),__LINE__  )
                DEALLOCATE(VERYTEMP)
            ENDIF
        ENDDO
        ENDDO
        CALL CHECK( NFmpi_close(NCUOTSSH),__LINE__  )

    ENDIF

if (0==rank) write(*,*) "DONE!"
call MPI_BARRIER(MPI_COMM_WORLD,IERR)



DEALLOCATE(SA_TE)

!call MPI_COMM_FREE(new_comm,ierr)

ENDDO

DEALLOCATE(WE_BOUNDARY,EA_BOUNDARY,SO_BOUNDARY,NO_BOUNDARY)
DEALLOCATE(PROC_WE,PROC_EA,PROC_SO,PROC_NO)

DEALLOCATE(X_ST,X_EN,Y_ST,Y_EN, XCST, XCEN,PROC_DIAG_SW,PROC_DIAG_NE)
DEALLOCATE(X2_ST,X2_EN,Y2_ST,Y2_EN,YMEN,YMST)
DEALLOCATE(MAP_XINI,MAP_YINI)
DEALLOCATE(MAP_XEND,MAP_YEND)

call MPI_Finalize(ierr)

END PROGRAM MAIN
