project(GraphicalModels)
subdirs(factors)

ADD_CXXTEST(factor_test.cxx)
add_graphlab_executable(loopybp_denoise loopybp_denoise.cpp)
requires_imagemagick(loopybp_denoise)



add_graphlab_executable(mplp_denoise mplp_denoise.cpp)
requires_eigen(mplp_denoise)


