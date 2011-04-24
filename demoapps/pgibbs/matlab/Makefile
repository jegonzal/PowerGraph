
gibbs_sampler_impl: src/gibbs_sampler_impl.cpp src/matwrap.hpp
	mex -largeArrayDims \
	CXXFLAGS="-g -fPIC -Wall -O3 -pthread -fexceptions -fno-omit-frame-pointer -fopenmp" \
	src/gibbs_sampler_impl.cpp \
	-I../../../src \
	-L../../../release/src/graphlab \
	-L../../../release/src/graphlab/extern/metis/GKlib \
	-L../../../release/src/graphlab/extern/metis/libmetis \
	-L../../../release/apps/pgibbs \
	-lpgibbs_pic \
	-lgraphlab_metis_pic -lgraphlab_GKlib_pic  \
	-lgraphlab_pic -lgraphlab_util_pic  -lgomp	



