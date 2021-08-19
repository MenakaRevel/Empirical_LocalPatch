# Fortran Compiler
FC  = ifort
#FC  = mpif90

# Compiler options
#FFLAGS = -O3 -warn all -assume byterecl -heap-arrays -nogen-interfaces
FFLAGS = -O3 -assume byterecl -heap-arrays -nogen-interfaces -mcmodel=large -fPIC -shared-intel -g -traceback -lpthread -mkl=parallel -parallel
# FFLAGS = -O3 -assume byterecl -heap-arrays -nogen-interfaces -mcmodel=large -fPIC -shared-intel -g -traceback -lpthread -mkl=parallel -qopenmp

## with netCDF
INC = -I/opt/intel/include -I/opt/local/include -I/opt/local/netcdf-fortran-4.4.5/include
LIB = -L/opt/intel/lib -L/opt/intel/mkl/lib -L/opt/local/lib -L/opt/local/netcdf-fortran-4.4.5/lib -lnetcdff

TARGETS=\
bin2nc \
remove_trend \
remove_season \
standardize \
semivariance \
lpara\
lparaMS

# Suffix
.SUFFIXES : .f90
.f90:
	$(FC) ${FFLAGS} ${INC} ${LIB} $^ -o $@

# commands
all: $(TARGETS)
# clean
clean:
	rm -rf  $(TARGETS)
