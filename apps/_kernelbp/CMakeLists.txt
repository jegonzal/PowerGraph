project(GraphLab)

add_executable(npdenoise nonparam_loopybp_denoising.cpp)
target_link_libraries(npdenoise tcmalloc)


add_executable(particle_bp particle_bp.cpp)
target_link_libraries(npdenoise tcmalloc)
