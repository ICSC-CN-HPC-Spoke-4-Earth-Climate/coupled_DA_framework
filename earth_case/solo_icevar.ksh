#!/bin/ksh
#-----------------------------------------------------------------------
# Run_3Dvar - Running 3D-Var Assimilation System
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

set -x

#--- 1. Initialize main variables

# If I am submitted I must know where am I, otherwise I am called...
export SUBMISSION_APART=1



#--- Environmental Variables affecting computational resources
export KMP_STACKSIZE="200M"
export MP_TASK_AFFINITY=-1
export KMP_AFFINITY=compact


EXPDIR=$1
EXPNAME=$2
WORKINGDIR=$3
END_DATE=$4'00'
VAR_LENGTH=$5
STATIC_DA_ICE=$6
IBDIR=$7
EOFICE_DIR=${8}
VAR_EXEDIR=${9}
NPROC_ICEVAR=${10}
NPROC_FRC=${11}
SIC_ARC=${12}
OSISAFARCHIV=${SIC_ARC}/OSISAF
HADISST=${SIC_ARC}/HADISST
PIOMASPATH=path/inputdata/SIT/PIOMAS/ARCHIVE/
GIOMASPATH=path/inputdata/SIT/GIOMAS/ARCHIVE/
SIT_ARC=${13}
CRYOSMOSARCHIV=${SIT_ARC}/CryoSat_SMOS_merged/
CRYO2ARCHIV=${SIT_ARC}/CryoSat_2/
SMOSARCHIV=${SIT_ARC}/SMOS/
DTG=${14}
PREV_DTG=${15}

export OMP_NUM_THREADS=$NPROC_ICEVAR

NAME_VAREXE=var_3d_juno_mpi_ice.x #
echo ""
echo "-----------------From 3DVAR script  ---------------"
echo " EXPDIR     :: " ${EXPDIR}
echo " EXPNAME    :: " ${EXPNAME}
echo " WORKDIR    :: " ${WORKINGDIR}
echo " END_DATE   :: " ${END_DATE}
echo " VAR_LENGTH :: " ${VAR_LENGTH}
echo " STATIC_DA  :: " ${STATIC_DA_ICE}
echo " IB_DIR     :: " ${IBDIR}
echo " EOFICE_DIR    :: " ${EOFICE_DIR}
echo " VAR_EXEDIR :: " ${VAR_EXEDIR}
echo " NPROC_ICEVAR  :: " ${NPROC_ICEVAR}
echo " NPROC_FRC  :: " ${NPROC_FRC}
echo " SIC_ARC    :: " ${SIC_ARC}
echo " OSISAFARCHIV :: " ${OSISAFARCHIV}
echo " SIT_ARC    :: " ${SIT_ARC}
echo " CRYOSMOSARCHIV :: " ${CRYOSMOSARCHIV}
echo " CRYO2ARCHIV :: " ${CRYO2ARCHIV}
echo " SMOSARCHIV :: " ${SMOSARCHIV}
echo " DTG        :: " ${DTG}
echo " PREV_DTG   :: " ${PREV_DTG}
echo "__________________________END________________________"

echo "ICE-3Dvar RUN" >whichstep
echo "ICE-3Dvar RUN : for ${DTG} on `date +%Y%m%d%H%M`" >>timestamp_run

. ${SIC_ARC}/../STATIC_Tools/myfunctions
#MYCONF=103
EXPER=${EXPNAME}

[[ -s $EXPDIR/config ]] || \
  { echo "No config file, exiting!" ; exit 2 ; }

if [ $DTG -gt $END_DATE ];then
	echo "Nothing to do for DTG ($DTG) > than END_DATE ($END_DATE)"
	exit -100
fi

WORKDIR=${WORKINGDIR}/assim_ice
rm -r ${WORKINGDIR}/assim_ice
mkdir -p ${WORKDIR}
#[[ -d $WORKDIR ]] && rm -rf $WORKDIR/*
#[[ -d $WORKDIR ]] || mkdir -p $WORKDIR

cd $WORKDIR || \
  { echo "Working dir ($WORKDIR) not found, aborting" ; exit 5 ; }

#--- 2. Prepare variables

ADDDTG="path/inputdata/STATIC_Tools/dateincr -h"

set -- `splitdate $DTG`
YYYY=$1 ; MM=$2 ; DD=$3 ; HH=$4 ; YMD=$5

HALF_TW=` expr $VAR_LENGTH / 2 `

START=`$ADDDTG $DTG -$HALF_TW`
set -- `splitdate $START`
SYYYY=$1 ; SMM=$2 ; SDD=$3 ; SHH=$4 ; SYMD=$5

END=`$ADDDTG $DTG +$HALF_TW`
set -- `splitdate $END`
EYYYY=$1 ; EMM=$2 ; EDD=$3 ; EHH=$4 ; EYMD=$5

#--- 3. Fetch ancillary files

SEASON=`whichseason $MM`

ln -fs $STATIC_DA_ICE/*.nc . || exit -2


NNMPI=1

if [ $NNMPI -gt 1 ]; then
        ip=0
        while [ $ip -lt $NNMPI ]; do
                nn=`printf "%04d" $ip`
                # cp $SCRDIR/domains_${NPROC_FRC}.dat.$nn domains.dat.$nn || \
                cp $STATIC_DA_ICE/OTHR/domains_${NPROC_FRC}.dat domains.dat.$nn || \
                { echo "Layout file not found, aborting" ; exit 3 ; }

                ip=$(( $ip + 1 ))
        done
        cp $STATIC_DA_ICE/OTHR/mpp_conf_${NPROC_FRC}_${NNMPI}.dat mpp_conf.dat || \
        { echo "Mpp Conf file not found, aborting" ; exit 3 ; }
        cp $STATIC_DA_ICE/OTHR/merge_domains_${NPROC_FRC}_${NNMPI}.dat merge_domains.dat || \
        { echo "Merge file not found, aborting" ; exit 3 ; }
fi

cp $STATIC_DA_ICE/OTHR/domains_${NPROC_FRC}.dat domains.dat || \
        { echo "Layout file not found, aborting" ; exit 3 ; }

PREV_YMD=`echo $PREV_DTG | cut -c 1-8`
PREV_WDIR=${WORKINGDIR}/${PREV_DTG}/run/
DTG_YMD=`echo $DTG | cut -c 1-8`
DTG_WDIR=${MISFIT_DIR}/${PREV_YMD}
echo "PREV_YMD" ${PREV_YMD}
echo "PREV_WDIR" ${PREV_WDIR}
echo "DTG_YMD" ${DTG_YMD}

# Shuffle observations
d1=`echo $START | cut -c 1-8`
d2=`echo $END   | cut -c 1-8`
id=$d1

eyyhad=$(echo $d2 | cut -c 1-4)
emmhad=$(echo $d2 | cut -c 5-6)
eddhad=$(echo $d2 | cut -c 7-8)


if [ -e ${HADISST}/HadISST_ice.nc ];then
  cdo -O -seldate,${eyyhad}-${emmhad}-01T00:00:00,${eyyhad}-${emmhad}-27T00:00:00 ${HADISST}/HadISST_ice.nc HADISST_${eyyhad}${emmhad}${eddhad}.nc
  ln -fs ${HADISST}/HadISST_ice_std.nc HadISST_ice_std.nc
fi
ncks -F -d time,$emmhad,$emmhad ${PIOMASPATH}/htot.H${eyyhad}.nc PIOMAS_${eyyhad}${emmhad}${eddhad}.nc
ln -fs ${PIOMASPATH}/piomas_std.nc piomas_std.nc

ncks -F -d time,$emmhad,$emmhad ${GIOMASPATH}/heff.H${eyyhad}.nc GIOMAS_${eyyhad}${emmhad}${eddhad}.nc
ln -fs ${GIOMASPATH}/giomas_std.nc giomas_std.nc



while [ $id -le $d2 ]; do
   yyosi=$(echo $id | cut -c 1-4)

   if [ -e ${OSISAFARCHIV}/${yyosi}/ice_sh_${id}.nc.gz ] || [ -e ${OSISAFARCHIV}/${yyosi}/ice_sh_${id}.nc ] ; then

     [[ -e ${OSISAFARCHIV}/${yyosi}/ice_sh_${id}.nc.gz ]] && cp ${OSISAFARCHIV}/${yyosi}/ice_sh_${id}.nc* . && gunzip  ice_sh_${id}.nc*
     [[ -e ${OSISAFARCHIV}/${yyosi}/ice_nh_${id}.nc.gz ]] && cp ${OSISAFARCHIV}/${yyosi}/ice_nh_${id}.nc* . && gunzip  ice_nh_${id}.nc*
     [[ -e ${OSISAFARCHIV}/${yyosi}/ice_sh_${id}.nc ]] && cp ${OSISAFARCHIV}/${yyosi}/ice_sh_${id}.nc .
     [[ -e ${OSISAFARCHIV}/${yyosi}/ice_nh_${id}.nc ]] && cp ${OSISAFARCHIV}/${yyosi}/ice_nh_${id}.nc .

     ln -sf ice_sh_${id}.nc OSISAF_sh_${id}.nc
     ln -sf ice_nh_${id}.nc OSISAF_nh_${id}.nc
    fi

    #if exists ice2 overwrite
   if [ -e ${OSISAFARCHIV}/${yyosi}/ice2_sh_${id}.nc.gz ] || [ -e ${OSISAFARCHIV}/${yyosi}/ice2_sh_${id}.nc ] ; then

     [[ -e ${OSISAFARCHIV}/${yyosi}/ice2_sh_${id}.nc.gz ]] && cp ${OSISAFARCHIV}/${yyosi}/ice2_sh_${id}.nc* . && gunzip  ice2_sh_${id}.nc*
     [[ -e ${OSISAFARCHIV}/${yyosi}/ice2_nh_${id}.nc.gz ]] && cp ${OSISAFARCHIV}/${yyosi}/ice2_nh_${id}.nc* . && gunzip  ice2_nh_${id}.nc*
     [[ -e ${OSISAFARCHIV}/${yyosi}/ice2_sh_${id}.nc ]] && cp ${OSISAFARCHIV}/${yyosi}/ice2_sh_${id}.nc .
     [[ -e ${OSISAFARCHIV}/${yyosi}/ice2_nh_${id}.nc ]] && cp ${OSISAFARCHIV}/${yyosi}/ice2_nh_${id}.nc .

     ln -sf ice2_sh_${id}.nc OSISAF_sh_${id}.nc
     ln -sf ice2_nh_${id}.nc OSISAF_nh_${id}.nc
    fi

fief=$(ls ${CRYOSMOSARCHIV}/${yyosi}/W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_${id}_*_r_v206_01_l4sit.nc)
[[ -e "CRYO2-SMOS_MERGED_nh_${id}.nc" ]] && rm CRYO2-SMOS_MERGED_nh_${id}.nc
[[ -e $fief ]] && cp $fief CRYO2-SMOS_MERGED_nh_${id}.nc #&& ncks -F -d xc,,,8 -d yc,,,5 CRYO2-SMOS_MERGED_nh_${id}.nc -O CRYO2-SMOS_MERGED_nh_${id}.nc

id=`date +%Y%m%d -d "$id +1 day"`
[[ $? -ne 0 ]] && id=`date +%Y%m%d -d "$id 12:00 +1 day"`
done

#--- 4. Prepare namelist and fetch aux files

SEASON=`whichseason $MM`


EOF_FILE=ICEOF_CICE_${MM}.nc #EOFC_${MM}.nc
EOF_FILE_LOC=EOF.nc
EOFNC=$EOF_FILE_LOC

ln -sf $EOFICE_DIR/$EOF_FILE $EOF_FILE_LOC

   #creating OPA_FIRSTGUESS file -->
   #_______________________________________________________________________
   start_date=${SYYYY}"-"${SMM}"-"${SDD}
   start_date1=${SYYYY}${SMM}${SDD}
   #$( date -d "$start_date " +%Y%m%d )    # rewrite in YYYYMMDD format

   #echo $start_date $start_date1
   end_date=${EYYYY}"-"${EMM}"-"${EDD}
   end_date1=${EYYYY}${EMM}${EDD}
   #echo $end_date $end_date1
   i=1
   run_date1=${start_date1}
   end_loop=${EYYYY}${EMMMM}${EDD}
   echo -ne "ncea " > merge.ksh
   #EXPNAME=NEMO-CICE-PROTO07
   #PREV_WDIR=path/archive/${EXPNAME}/ice/hist/
   while [ "${run_date1}" -lt "${end_date1}" ]; do
             ttem=$run_date1
   	     run_date1=`date +%Y%m%d -d "$ttem +1 day"`
             [[ $? -ne 0 ]] && run_date1=`date +%Y%m%d -d "$ttem 12:00 +1 day"`
             #[[ $? -ne 0 ]] && run_date1=$( date -d "$run_date " +%Y%m%d )

	    #         echo $start_date1 $end_date1
   #         echo $run_date1 $end_date1
   
            DTGtemp=$run_date1'00'
            set -- `splitdate $DTGtemp`
            YYYYtemp=$1 ; MMtemp=$2 ; DDtemp=$3 ; HHtemp=$4 ; YMDtemp=$5

            ln -fs ${PREV_WDIR}/${EXPNAME}.cice.h.${YYYYtemp}-${MMtemp}-${DDtemp}.nc ${EXPNAME}.cice.h.${YYYYtemp}-${MMtemp}-${DDtemp}.nc
            file_name=${EXPNAME}.cice.h.${YYYYtemp}-${MMtemp}-${DDtemp}.nc
            echo -ne "${file_name}  " >> merge.ksh
            i=$(( i + 1 ))
   done
   echo -ne "tempor_icef1.nc && " >> merge.ksh
   echo -n "cdo selvar,aice,hi tempor_icef1.nc OPA_FIRSTGUESS.nc && rm tempor_icef1.nc && " >> merge.ksh
   echo -n "ncrename -d time,time_counter -d nj,y -d ni,x -v hi,icethic -v aice,soicecov OPA_FIRSTGUESS.nc" >> merge.ksh
   chmod +x merge.ksh
   ./merge.ksh
   rm ${EXPNAME}.cice.h.*
   ln -fs OPA_FIRSTGUESS.nc OPA_FIRSTGUESS_ice.nc
#EXPNAME=NEMO-CICE-PROTO09
#_______________________________________________________________________
#--- The Namelist

if [ $NNMPI -gt 1 ]; then
	USEMPI=.TRUE.
else
	USEMPI=.FALSE.
fi

#fr=$SCRDIR/ratio_bgerr.dat
REDNMC=1
#REDNMC=`grep ${YMD} $fr | awk '{print $2};'`
#echo ${YMD}
#echo $REDNMC

ANAMORF=1
#EXPDIR=/users_home/csp/aspect/CESM2/Exps/NEMO-CICE-PROTO09
VAR_NAMELIST=${EXPDIR}/namelist_3dvar_ice
#EXPDIR=/users_home/csp/aspect/CESM2/Exps/NEMO-CICE-PROTO07

sed -e "s/_EOF_FILE_/$EOFNC/g" \
  -e "s/_USEMPI_/$USEMPI/g" -e "s/_NNMPI_/$NNMPI/g" \
   $VAR_NAMELIST > var_3d_nml || \
  { echo "No 3DVAR namelist, aborting!" ; exit 12 ; }

#--- 5. Run the model

cp ${VAR_EXEDIR}/${NAME_VAREXE} ./MASTER || \
  { echo "No 3DVAR executable, aborting!" ; exit 12 ; }

chmod +x MASTER

time mpiexec.hydra  -n 1 -ppn 1 ./MASTER -d ${YMD} -t 00 -w $VAR_LENGTH > \
  3DVAR_ICE_${DTG}.out 2>&1 || \
  { echo "3DVAR returned non-zero status, aborting!" ; exit 12 ; }
ls ANINCR.NC* || \
  { echo "3DVAR failed, aborting!" ; exit 12 ; }


if [ $NNMPI -gt 1 ]; then
	mpiexec.hydra -n $NNMPI  $MERGED_EXE || exit -12
	if [ $ANAMORF -eq 1 ]; then
	ncrename -v INCICETHI,INCICETHI_in_10m ANINCR_MERGED.NC
	ncap2 -s "INCICETHI=float(INCICETHI_in_10m*10.)" ANINCR_MERGED.NC ANINCR.NC
	else	
	ln -sf ANINCR_MERGED.NC ANINCR.NC
	fi
else

if [ $ANAMORF -eq 1 ]; then
        ncrename -v INCLEADFR,INCICECOV -v INCICETHI,INCICETHI_in_10m ANINCR.NC
        ncap2 -s "INCICETHI=float(INCICETHI_in_10m*10.)" ANINCR.NC -O ANINCR.NC
fi
fi

n=0
for f in LOG_3DVAR*; do
  grep "NaN" $f > /dev/null && n=$(( $n + 1 ))
done

if [ $n -gt 0 ]; then
   echo "3DVAR returned NaN"
   exit -999
   echo "$DTG : 3DVAR returned NaN" >> ${EXPDIR}/assim.failed
   OLDDTG=`$ADDDTG $PREV_DTG -$VAR_LENGTH`
   cd $SCRDIR
#   cat > $SCRDIR/config <<EOF
#${PREV_DATE}
#${OLDDTG}
#EOF
#   $SCHEDULER < run_orca
   exit 0
fi

mv var_3d_nml var_3d_nml.3DVAR
cd ${WORKINGDIR}
mkdir -p ${YMD}00 && mv -f assim_ice ${YMD}00/ && mkdir -p assim_ice

exit 0

#--- End
#-----------------------------------------------------------------------
