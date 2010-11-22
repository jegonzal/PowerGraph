project(GraphLab)

# Shared memory CoEM for multiclass inference
add_executable(multicoem  multicoemapp.cpp)
target_link_libraries(multicoem tcmalloc)


