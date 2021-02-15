-include make.inc

.PHONY: all fasttransforms fmm3d shtns mwrap mex clean

default: all

all: fasttransforms fmm3d shtns

fasttransforms:
	make -C +backends/+fasttransforms/

fmm3d:
	make -C +backends/+fmm3d/

shtns:
	make -C +backends/+shtns/

# Rebuild interfaces via MWrap:
mwrap:
	make -C +backends/+fasttransforms/ mwrap
	make -C +backends/+fmm3d/ mwrap
	make -C +backends/+shtns/ mwrap

# Rebuild MEX files via MEX:
mex:
	make -C +backends/+fasttransforms/ mex
	make -C +backends/+fmm3d/ mex
	make -C +backends/+shtns/ mex

# Remove the MEX files:
clean:
	make -C +backends/+fasttransforms/ clean
	make -C +backends/+fmm3d/ clean
	make -C +backends/+shtns/ clean

# Remove the MEX interface, MATLAB caller, and MEX files:
# Note: You will need MWrap to rebuild the deleted files!
mwrapclean:
	make -C +backends/+fasttransforms/ mwrapclean
	make -C +backends/+fmm3d/ mwrapclean
	make -C +backends/+shtns/ mwrapclean
