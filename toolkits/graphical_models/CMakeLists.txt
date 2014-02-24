project(GraphicalModels)


add_graphlab_executable(dd dd_main.cpp)
requires_eigen(dd)

add_graphlab_executable(lbp_structured_prediction lbp_structured_prediction.cpp)
requires_eigen(lbp_structured_prediction)

add_graphlab_executable(profile_lbp_synthetic profile_lbp_synthetic.cpp)
requires_eigen(profile_lbp_synthetic)


add_graphlab_executable(profile_lbp_synthetic2 profile_lbp_synthetic2.cpp)
requires_eigen(profile_lbp_synthetic2)


add_graphlab_executable(synthetic_image_data synthetic_image_data.cpp)
requires_opencv(synthetic_image_data)




add_graphlab_executable(mplp_denoise mplp_denoise.cpp)
requires_eigen(mplp_denoise)
requires_opencv(mplp_denoise)


subdirs(factors)
