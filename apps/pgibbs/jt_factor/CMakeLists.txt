project(GraphLab)


add_executable(pgibbs_jt_alchemy pgibbs_jt_alchemy.cpp)
# TARGET_LINK_LIBRARIES(pgibbs_jt_alchemy tcmalloc)

add_executable(test_jt test_jt.cpp)

add_executable(alchemy_image_denoise alchemy_image_denoise.cpp)

