-include ../../Make.inc

CFLAGS = -fopenmp -O3 -march=native -std=gnu11 -fpic
MEXFLAGS = -largeArrayDims -DR2008OO CFLAGS='$(CFLAGS)'
LIB = -lfasttransforms -lm -lmpfr -lblas -lfftw3 -lfftw3_threads
LIBDIR = -L$(FASTTRANSFORMS_DIR) -L/usr/local/opt/openblas/lib -L/usr/local/opt/fftw/lib -L/usr/local/opt/mpfr/lib -L/usr/local/opt/gmp/lib
INCDIR = -I$(FASTTRANSFORMS_DIR)/src -I/usr/local/opt/openblas/include -I/usr/local/opt/fftw/include -I/usr/local/opt/mpfr/include -I/usr/local/opt/gmp/include

.PHONY: mwrap mex clean mwrapclean

default: mwrap mex

mwrap: mex.c

# Rebuild mex.c and backend.m via MWrap
mex.c: backend.mw
	$(MWRAP) -mex backends.fasttransforms.mex -c mex.c -mb -c99complex backend.mw

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
