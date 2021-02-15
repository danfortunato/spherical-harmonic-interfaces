-include ../../Make.inc

CFLAGS = -fopenmp -O3 -march=native -mtune=skylake -ffast-math -std=gnu99 -fpic -D_GNU_SOURCE
LDFLAGS = -flat_namespace -undefined warning
MEXFLAGS = -largeArrayDims -DR2008OO CFLAGS='$(CFLAGS)' LDFLAGS='$(LDFLAGS)'
LIB = -lshtns_omp -lfftw3_omp -lfftw3 -lm
LIBDIR = -L$(SHTNS_DIR) -L/usr/local/lib -L/usr/local/opt/fftw/lib
INCDIR = -I$(SHTNS_DIR) -I/usr/local/opt/fftw/include

.PHONY: mwrap mex clean mwrapclean

default: mwrap mex

mwrap: mex.c

# Rebuild mex.c and backend.m via MWrap
mex.c: backend.mw
	$(MWRAP) -mex backends.shtns.mex -c mex.c -mb -c99complex backend.mw

mex: mex.mex*

# Rebuild mex.mex* via MEX
mex.mex*: mex.c
	$(MEX) mex.c $(INCDIR) $(MEXFLAGS) $(LIBDIR) $(LIB) -output mex

# Remove the MEX file:
clean:
	rm -f *.mex*

# Remove the MEX interface, MATLAB caller, and MEX file:
# Note: You will need MWrap to rebuild the deleted files!
mwrapclean:
	rm -f backend.m mex.c *.mex*
