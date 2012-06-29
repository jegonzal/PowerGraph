project(GraphicalModels)
subdirs(factors)

ADD_CXXTEST(factor_test.cxx)

add_graphlab_executable(loopybp_denoise loopybp_denoise.cpp)
requires_opencv(loopybp_denoise)


add_graphlab_executable(lbp_structured_prediction lbp_structured_prediction.cpp)
requires_eigen(lbp_structured_prediction)

add_graphlab_executable(synthetic_image_data synthetic_image_data.cpp)
requires_opencv(synthetic_image_data)




add_graphlab_executable(mplp_denoise mplp_denoise.cpp)
requires_eigen(mplp_denoise)
requires_opencv(mplp_denoise)

