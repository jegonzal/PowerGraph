project(GraphLab)

#add_executable(npdenoise nonparam_loopybp_denoising.cpp)
#target_link_libraries(npdenoise tcmalloc)


add_executable(particle_bp particle_bp.cpp)
target_link_libraries(particle_bp tcmalloc)


add_executable(particle_bp_itgmm particle_bp_itgmm.cpp image_compare.cpp)
target_link_libraries(particle_bp_itgmm tcmalloc itpp)
#set_target_properties(particle_bp_itgmm PROPERTIES COMPILE_FLAGS "-ffast-math")


add_executable(discrete_bp_itgmm discrete_bp_itgmm.cpp image_compare.cpp)
target_link_libraries(discrete_bp_itgmm tcmalloc itpp)

add_executable(kbp kbp.cpp image_compare.cpp)
target_link_libraries(kbp tcmalloc itpp)


add_executable(particle_bp_make3d particle_bp_make3d.cpp)
target_link_libraries(particle_bp_make3d tcmalloc itpp)

add_executable(discrete_bp_make3d discrete_bp_make3d.cpp)
target_link_libraries(discrete_bp_make3d tcmalloc itpp)

add_executable(kbp_make3d kbp_make3d.cpp)
target_link_libraries(kbp_make3d tcmalloc itpp)

#add_executable(pbp_protein pbp_protein.cpp)
#target_link_libraries(pbp_protein tcmalloc itpp)

#add_executable(kbp_protein kbp_protein.cpp)
#target_link_libraries(kbp_protein tcmalloc itpp)

#add_executable(kbplinear kbplinear.cpp image_compare.cpp)
#target_link_libraries(kbplinear tcmalloc itpp)

#add_executable(kbplinear_protein kbplinear_protein.cpp)
#target_link_libraries(kbplinear_protein tcmalloc itpp)
