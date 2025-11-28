#!/bin/bash
#-----------------------------------------------------------------------
# Run_3Dvar - Running 3D-Var Assimilation System
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

#echo "Hello! from PHY-3Dvar"

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
STATIC_DA=$6
MISFIT_DIR=$7
RFC_DIR=$8
IBDIR=$9
EOF_DIR=${10}
VAR_EXEDIR=${11}
NPROC_VAR=${12}
NPROC_FRC=${13}
DTG=${14}
PREV_DTG=${15}
MERGED_EXE=${VAR_EXEDIR}/merge.x

export OMP_NUM_THREADS=$NPROC_VAR

NAME_VAREXE="var_3d_juno_mpi.x"
echo ""
echo "-----------------From 3DVAR script  ---------------"
echo " EXPDIR     :: " ${EXPDIR}
echo " EXPNAME    :: " ${EXPNAME}
echo " WORKDIR    :: " ${WORKINGDIR}
echo " END_DATE   :: " ${END_DATE}
echo " VAR_LENGTH :: " ${VAR_LENGTH}
echo " STATIC_DA  :: " ${STATIC_DA}
echo " MISFIT_DIR :: " ${MISFIT_DIR}
echo " RFC_DIR    :: " ${RFC_DIR}
echo " IB_DIR     :: " ${IBDIR}
echo " EOF_DIR    :: " ${EOF_DIR}
echo " VAR_EXEDIR :: " ${VAR_EXEDIR}
echo " NPROC_VAR  :: " ${NPROC_VAR}
echo " NPROC_FRC  :: " ${NPROC_FRC}
echo " DTG        :: " ${DTG}
echo " PREV_DTG   :: " ${PREV_DTG}
echo " MERGED_EXE :: " ${MERGED_EXE}
echo "__________________________END________________________"



echo "3Dvar RUN" >whichstep
echo "3Dvar RUN : for ${DTG} on `date +%Y%m%d%H%M`" >>timestamp_run

. ${EXPDIR}/myfunctions
MYCONF=103
EXPER=${EXPNAME}

[[ -s $EXPDIR/config ]] || \
  { echo "No config file, exiting!" ; exit 2 ; }

#set -- `cat ${EXPDIR}/config`
#DTG=$1
#PREV_DATE=$2

if [ $DTG -gt $END_DATE ];then
	echo "Nothing to do for DTG ($DTG) > than END_DATE ($END_DATE)"
	exit -100
fi

WORKDIR=${WORKINGDIR}/assim_phy

cd $WORKDIR || \
  { echo "Working dir ($WORKDIR) not found, aborting" ; exit 5 ; }

#--- 2. Prepare variables

ADDINCR="path/inputdata/STATIC_Tools/dateincr -h"

set -- `splitdate $DTG`
YYYY=$1 ; MM=$2 ; DD=$3 ; HH=$4 ; YMD=$5

HALF_TW=` expr $VAR_LENGTH / 2 `

START=`$ADDINCR $DTG -$HALF_TW`
set -- `splitdate $START`
SYYYY=$1 ; SMM=$2 ; SDD=$3 ; SHH=$4 ; SYMD=$5

END=`$ADDINCR $DTG +$HALF_TW`
set -- `splitdate $END`
EYYYY=$1 ; EMM=$2 ; EDD=$3 ; EHH=$4 ; EYMD=$5

#--- 3. Fetch ancillary files

SEASON=`whichseason $MM`

cp  $STATIC_DA/*.nc .
cp  $STATIC_DA/*.dat .
NNMPI=6
#NPROC_FRC=330 #396 #To be changed according to ocean component

if [ $NNMPI -gt 1 ]; then
	ip=0
	while [ $ip -lt $NNMPI ]; do
		nn=`printf "%04d" $ip`
		# cp $SCRDIR/domains_${NPROC_FRC}.dat.$nn domains.dat.$nn || \
		cp $STATIC_DA/OTHR/domains_${NPROC_FRC}.dat domains.dat.$nn || \
	        { echo "Layout file not found, aborting" ; exit 3 ; }

		ip=$(( $ip + 1 ))
	done

        cp $STATIC_DA/OTHR/mpp_conf_${NPROC_FRC}_${NNMPI}.dat mpp_conf.dat || \
        { echo "Mpp Conf file not found, aborting" ; exit 3 ; }
#	cp $STATIC_DA/OTHR/merge_domains_${NPROC_FRC}_${NNMPI}.dat merge_domains.dat || \
#        { echo "Merge file not found, aborting" ; exit 3 ; }
fi

cp $STATIC_DA/OTHR/domains_${NPROC_FRC}.dat domains.dat || \
        { echo "Layout file not found, aborting" ; exit 3 ; }

PREV_YMD=`echo $PREV_DTG | cut -c 1-8`
DTG_YMD=`echo $DTG | cut -c 1-8`
DTG_WDIR=${MISFIT_DIR}/${PREV_YMD}00/run
echo "PREV_YMD" ${PREV_YMD}
echo "PREV_WDIR" ${PREV_WDIR}
echo "DTG_YMD" ${DTG_YMD}
cp ${DTG_WDIR}/*MIS_*NC .
cp ${DTG_WDIR}/OBS_TAB.DAT OBS_TAB.DAT

if [ $NNMPI -gt 1 ]; then
        USEMPI=.TRUE.
else
        USEMPI=.FALSE.
fi

fr=${STATIC_DA}/ratio_bgerr_full.dat #ratio_bgerr_new.dat
REDNMC=`grep ${YMD} $fr | awk '{print $2};'`


VAR_NAMELIST=${EXPDIR}/namelist_3dvar
cp ${VAR_EXEDIR}/merge.x .
cp ${VAR_EXEDIR}/shuffle_obs_nemo42.x shuffle_obs.x
SHUFFLE_OBS=${WORKDIR}/shuffle_obs.x

# Shuffle observations
echo " Running Shufle obs ... "

 $SHUFFLE_OBS $NPROC_FRC || exit -12
#--------------------------------------------------------------
#recontruct

module purge
module purge

eval "$(command /juno/opt/anaconda/3-2022.10/bin/conda 'shell.bash' 'hook' 2> /dev/null)"

conda activate cmcc_rebuild
loop_date=$PREV_DTG
while [ "$loop_date" -le "$DTG" ]; do
#ff=$(ls -1 ${DTG_WDIR}/${EXPNAME}_1d*_grid_T_*.nc  | head -1)
ff=$(ls -1 ${DTG_WDIR}/${EXPNAME}_1d*_grid_T_${loop_date:0:8}*.nc  | head -1)

mpirun -n 66 python /work/cmcc/c-glors/GIT/py_nemo_rebuild/src/py_nemo_rebuild/nemo_rebuild.py -i $ff
loop_date=$(date -d "${loop_date:0:8} ${loop_date:8:2}:00 + 1 day" +%Y%m%d%H)

done



conda deactivate

module purge
module purge
module load gcc-12.2.0 R/4.2.2

module load oneapi-2022.1.0/compiler-rt/2022.1.0 intel-2021.6.0/2021.6.0 impi-2021.6.0/2021.6.0 intel-2021.6.0/cdo-threadsafe intel-2021.6.0/nco intel-2021.6.0/ncview intel-2021.6.0/impi-2021.6.0/netcdf-c-threadsafe  intel-2021.6.0/impi-2021.6.0/netcdf-cxx-threadsafe intel-2021.6.0/impi-2021.6.0/netcdf-fortran-threadsafe intel-2021.6.0/impi-2021.6.0/parallel-netcdf  intel-2021.6.0/curl intel-2021.6.0/impi-2021.6.0/hdf5-threadsafe intel-2021.6.0/impi-2021.6.0/xios intel-2021.6.0/magics-threadsafe intel-2021.6.0/eccodes-threadsafe intel-2021.6.0/cmake oneapi-2022.1.0/tbb/2021.6.0 oneapi-2022.1.0/mkl/2022.1.0 intel-2021.6.0/tk/8.6.11-wa3ya

ffbasenam=$(basename $ff)


ddtemp=${PREV_YMD}
time_counter=1
loop_date=$PREV_DTG

while [ "$loop_date" -le "$DTG" ]; do
#do
tempname=$(ls -1 ${DTG_WDIR}/${EXPNAME}_1d_*_grid_T_${loop_date:0:8}-${loop_date:0:8}.nc)
ncks -F -d time_counter,1 -v zos  $tempname  -O ${DTG_WDIR}/../../SSH/${EXPNAME}_SSH_${loop_date:0:8}.nc #|| exit -6
loop_date=$(date -d "${loop_date:0:8} ${loop_date:8:2}:00 + 1 day" +%Y%m%d%H)
#ddtemp=$(date +%Y%m%d -d "${ddtemp} +1 day")
done

#conda activate base
conda activate /work/cmcc/c-glors/.conda/cesm
test_date=$(date +%Y%m%d -d "${DTG_YMD} -10 day")
ddtemp1=$(date +%Y%m%d -d "${test_date} -12 day")
ddtemp2=$(date +%Y%m%d -d "${test_date} +2 day")

#ddtemp1=$(date +%Y%m%d -d "${DTG_YMD} -12 day")
#ddtemp2=$(date +%Y%m%d -d "${DTG_YMD} +2 day")
python ${EXPDIR}/calc_bias.py -s ${ddtemp1}  -e ${ddtemp2} --arch-gofs-mesh "path/inputdata/STATIC_DA_NEMO42/assim/"  --dir-sla-unbias "path/${EXPNAME}/SSH/" > out_py 2>&1

mv bias.nc path/${EXPNAME}/BIAS/SLA_BIAS_OMG_${DTG_YMD}.nc

conda deactivate


module purge
module purge
module load gcc-12.2.0 R/4.2.2

module load oneapi-2022.1.0/compiler-rt/2022.1.0 intel-2021.6.0/2021.6.0 impi-2021.6.0/2021.6.0 intel-2021.6.0/cdo-threadsafe intel-2021.6.0/nco intel-2021.6.0/ncview intel-2021.6.0/impi-2021.6.0/netcdf-c-threadsafe  intel-2021.6.0/impi-2021.6.0/netcdf-cxx-threadsafe intel-2021.6.0/impi-2021.6.0/netcdf-fortran-threadsafe intel-2021.6.0/impi-2021.6.0/parallel-netcdf  intel-2021.6.0/curl intel-2021.6.0/impi-2021.6.0/hdf5-threadsafe intel-2021.6.0/impi-2021.6.0/xios intel-2021.6.0/magics-threadsafe intel-2021.6.0/eccodes-threadsafe intel-2021.6.0/cmake oneapi-2022.1.0/tbb/2021.6.0 oneapi-2022.1.0/mkl/2022.1.0 intel-2021.6.0/tk/8.6.11-wa3ya


export WORK=/path
export DATA=/path

export CESMDATAROOT=/data/inputs//


#--- 4. Prepare namelist and fetch aux files

SEASON=`whichseason $MM`

ln -sf $RFC_DIR/Corrad_R025L75_onx_${SEASON}.nc Corrad_onx.nc
ln -sf $RFC_DIR/Corrad_R025L75_ony_${SEASON}.nc Corrad_ony.nc


ln -sf ${IBDIR}/IB_error_${YYYY}.nc IB_error.nc

EOF_OPT=CTRL03BR_R05L75
EOF_FILE=EOF_${EOF_OPT}_${SEASON}.nc
EOF_FILE_LOC=EOF.nc
EOFNC=$EOF_FILE_LOC

ln -fs $EOF_DIR/$EOF_FILE $EOF_FILE_LOC


ncra -v thetao,so ${DTG_WDIR}/${EXPNAME}_1d*_grid_T_*-????????.nc -O BACKGROUND.nc
ncrename -v thetao,votemper -v so,vosaline BACKGROUND.nc

#--- The Namelist

sed -e "s/_NSEED_/$MEMB/g" -e "s/_NCONF_/$MYCONF/g" \
  -e "s/_NCPU_/$NPROC_FRC/g" -e "s/_EOF_FILE_/$EOFNC/g" \
  -e "s/_USEMPI_/$USEMPI/g" -e "s/_NNMPI_/$NNMPI/g" \
  -e "s/_REDNMC_/$REDNMC/g" \
   $VAR_NAMELIST > var_3d_nml || \
  { echo "No 3DVAR namelist, aborting!" ; exit 12 ; }


#--- 5. Run the model

cp ${VAR_EXEDIR}/${NAME_VAREXE} ./MASTER || \
  { echo "No 3DVAR executable, aborting!" ; exit 12 ; }

chmod +x MASTER

time mpirun  -n 6  ./MASTER -d ${YMD} -t 00 -w $VAR_LENGTH > \
  3DVAR_${DTG}.out 2>&1 || \
  { echo "3DVAR returned non-zero status, aborting!" ; exit 12 ; }
ls ANINCR.NC* || \
  { echo "3DVAR failed, aborting!" ; exit 12 ; }

if [ $NNMPI -gt 1 ]; then
mpirun -n 6	$MERGED_EXE || exit -12
mv ANINCR_T_MERGED.nc ANINCR_MERGED.NC
ncks ANINCR_S_MERGED.nc -A ANINCR_MERGED.NC
ncks ANINCR_SSH_MERGED.nc -A ANINCR_MERGED.NC

ncap2 -F -s "INCTEMPER(53,:,:)=INCTEMPER(53,:,:)*0.8;INCSALINE(53,:,:)=INCSALINE(53,:,:)*0.8"  ANINCR_MERGED.NC -O ANINCR_TMP.NC
ncap2 -F -s "INCTEMPER(54,:,:)=INCTEMPER(54,:,:)*0.1;INCSALINE(54,:,:)=INCSALINE(54,:,:)*0.1"  ANINCR_MERGED.NC -O ANINCR_TMP.NC
ncap2 -F -s "INCTEMPER(55:75,:,:)=0.;INCSALINE(55:75,:,:)=0."  ANINCR_TMP.NC -O ANINCR.NC && rm -f ANINCR_TMP.NC
fi


n=0
for f in LOG_3DVAR*; do
  grep "NaN" $f > /dev/null && n=$(( $n + 1 ))
done

if [ $n -gt 0 ]; then
   echo "3DVAR returned NaN"
   exit -999
   echo "$DTG : 3DVAR returned NaN" >> ${EXPDIR}/assim_phy.failed
   OLDDTG=`$ADDDTG $PREV_DATE -$VAR_LENGTH`
   cd $SCRDIR
   exit 0
fi

mv var_3d_nml var_3d_nml.3DVAR
cd ${WORKINGDIR}
mkdir -p ${YMD}00 && mv assim_phy ${YMD}00/ && mkdir -p assim_phy

exit 0

#--- End
#-----------------------------------------------------------------------
