project(GraphLab)
add_graphlab_executable(kcores kcores.cpp ../shared/mmio.cpp)
add_graphlab_executable(kcores_async kcores_async.cpp ../shared/mmio.cpp)
