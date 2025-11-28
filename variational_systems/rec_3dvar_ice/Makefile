######################## 
#module load intel19 impi19 intel19/szip impi19/hdf5-threadsafe/1.10.5 impi19/netcdf-threadsafe/C_4.7.0-F_4.4.6_CXX_4.3.0 impi19/parallel-netcdf
# ATHENA intel compilers
#
FC=mpiifort
FFLAGS= -qopt-zmm-usage=high -fp-model source -I${PNETCDF}/include
LINK=  -L${PARALLEL_NETCDF}/lib -lpnetcdf
#
#
#
#
OBJECTS= set_knd.o check.o sort.o
#
merge_mpi_pnet.exe: merge_mpi_pnet.F90 ${OBJECTS} 
	${FC} $(FFLAGS) merge_mpi_pnet.F90 $(OBJECTS) $(LINK) -o merge_mpi_pnet.exe

%.o:%.F90
	${FC} $(FFLAGS) -c $<
clean:
	rm -f $(OBJECTS) $(DEVELOBJS)
