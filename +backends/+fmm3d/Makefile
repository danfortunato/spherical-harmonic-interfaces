-include ../../Make.inc

FFLAGS = -fallow-argument-mismatch -fpic -O2
MEXFLAGS = -largeArrayDims -lm -lgfortran -DR2008OO

# Fortran sources
FSRCS = spharmrouts.f yrecursion.f
OBJS  = $(FSRCS:.f=.o)

.PHONY: mwrap mex clean mwrapclean

default: mwrap mex

mwrap: mex.c

.f.o: Makefile
	$(FC) $(FFLAGS) -c $<

# Rebuild mex.c and backend.m via MWrap
mex.c: backend.mw
	$(MWRAP) -mex backends.fmm3d.mex -c mex.c -mb -c99complex backend.mw

mex: mex.mex*

# Rebuild mex.mex* via MEX
mex.mex*: mex.c $(OBJS)
	$(MEX) mex.c $(OBJS) $(MEXFLAGS) -output mex

# Remove the MEX file:
clean:
	rm -f *.o *.mex*

# Remove the MEX interface, MATLAB caller, and MEX file:
# Note: You will need MWrap to rebuild the deleted files!
mwrapclean:
	rm -f backend.m mex.c *.o *.mex*
