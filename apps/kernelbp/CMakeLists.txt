project(GraphLab)

#add_executable(npdenoise nonparam_loopybp_denoising.cpp)
#target_link_libraries(npdenoise tcmalloc)


add_executable(particle_bp particle_bp.cpp)
target_link_libraries(particle_bp tcmalloc)


add_executable(particle_bp_itgmm particle_bp_itgmm.cpp)
target_link_libraries(particle_bp_itgmm tcmalloc itpp)
#set_target_properties(particle_bp_itgmm PROPERTIES COMPILE_FLAGS "-ffast-math")
