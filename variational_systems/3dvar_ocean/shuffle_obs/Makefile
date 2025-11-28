F90 = ifort -convert big_endian -O2 -i4 -r8 -fno-alias -assume byterecl 
DEBUG= -g -traceback -check bounds
L90 = $(F90)
NCD = $(NETCDF_FORTRAN_DIR)
BIN = shuffle_obs_nemo42.x
BIN2 = shuffle_obs_pre_nemo42.x
EXTRALIBS = 

LIB = setfgflayout.o  abor1.o io_obs.o mod_obsdef.o obs_str.o set_dim.o setfgfiles.o \
check.o mod_allocobs.o mod_allocobs2.o mod_readfg.o readobs103.o handle_err.o set_knd.o

LIB2 = abor1.o io_obs.o mod_obsdef.o obs_str.o set_dim.o setfgfiles.o \
check.o mod_allocobs.o mod_allocobs2.o mod_readfg.o set_knd.o \
setfgflayout.o readobs103_pre.o mod_phinterp.o int_par_ins_2.o mod_mynetcdf.o index_sort.o handle_err.o

%.o:	%.F90
	$(F90) -c -o $@ $< -I$(NCD)/include

all:	$(BIN) $(BIN2)

$(BIN):	$(LIB)
	$(L90) -o $(BIN) $(LIB) -L$(NCD)/lib -lnetcdf -lnetcdff $(EXTRALIBS)

$(BIN2):	$(LIB2)
	$(L90) -o $(BIN2) $(LIB2) -L$(NCD)/lib -lnetcdf -lnetcdff $(EXTRALIBS)

clean:	
	rm -f *.o *.mod $(BIN) $(BIN2)

include dependencies.inc
