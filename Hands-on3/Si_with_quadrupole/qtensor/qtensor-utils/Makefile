FC = mpifort

H5_ROOT = /Users/lbkhoa/softwares/libraries/hdf5-1.14.6-gcc-mpi

FCFLAGS = -I$(H5_ROOT)/include

LDFLAGS = -L$(H5_ROOT)/lib

HDF5_DIR = $(H5_ROOT)/lib

HDF5_LIB = $(HDF5_DIR)/libhdf5_fortran.a $(HDF5_DIR)/libhdf5_hl_fortran.a \
		   $(HDF5_DIR)/libhdf5_hl.a $(HDF5_DIR)/libhdf5.a

H5_LIB   = $(H5_ROOT)/lib/libh5fortran.a

LIB = -lpthread -ldl -lm -lz

LIBS = $(H5_LIB) $(HDF5_LIB) $(LIB)

all: qtensor

qtensor: qtensor.o
	@$(FC) -o $@ $< $(LDFLAGS)  $(LIBS)
	@echo "gfortran -o qtensor qtensor.o"

qtensor.o: qtensor.f90
	@$(FC) $(FCFLAGS) -c $<
	@echo "gfortran -c qtensor.f90"

clean:
	rm -f qtensor qtensor.o
