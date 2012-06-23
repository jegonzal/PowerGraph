project(GraphicalModels)
subdirs(factors)

ADD_CXXTEST(factor_test.cxx)

add_graphlab_executable(loopybp_denoise loopybp_denoise.cpp)
requires_imagemagick(loopybp_denoise)


add_graphlab_executable(structured_prediction structured_prediction.cpp)
requires_eigen(structured_prediction)

add_graphlab_executable(synthetic_image_data synthetic_image_data.cpp)
requires_imagemagick(synthetic_image_data)




add_graphlab_executable(mplp_denoise mplp_denoise.cpp)
requires_eigen(mplp_denoise)
requires_imagemagick(mplp_denoise)

