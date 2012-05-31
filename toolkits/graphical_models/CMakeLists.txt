project(GraphicalModels)
subdirs(factors)

ADD_CXXTEST(factor_test.cxx)
add_graphlab_executable(loopybp_denoise loopybp_denoise.cpp)
# add_dependencies(loopybp_denoise imagemagick)
# target_link_libraries(loopybp_denoise  Magick++)


add_graphlab_executable(mplp_denoise mplp_denoise.cpp)
# add_dependencies(mplp_denoise imagemagick eigen)
# target_link_libraries(mplp_denoise  Magick++)

