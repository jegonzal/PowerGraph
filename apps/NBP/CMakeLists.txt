project(GraphLab)


add_executable(NBP NBP.cpp BallTreeDensity.cpp
cpp/BallTreeClass.cc cpp/BallTreeDensityClass.cc)
target_link_libraries(NBP tcmalloc)
target_link_libraries(NBP itpp)
add_executable(NBPD denoise.cpp BallTreeDensity.cpp
cpp/BallTreeClass.cc cpp/BallTreeDensityClass.cc)
target_link_libraries(NBP tcmalloc)
target_link_libraries(NBPD itpp_debug)
