
gibbs_sampler_impl: gibbs_sampler_impl.cpp matwrap.hpp
	mex -largeArrayDims \
	CXXFLAGS="-g -fPIC -Wall -O3 -pthread -fexceptions -fno-omit-frame-pointer -fopenmp" \
	gibbs_sampler_impl.cpp \
	-I../../../src \
	-L../../../release/src/graphlab \
	-L../../../release/src/graphlab/extern/metis/GKlib \
	-L../../../release/src/graphlab/extern/metis/libmetis \
	-L../../../release/demoapps/pgibbs \
	-lpgibbs_pic \
	-lgraphlab_pic \
	-lgomp	



