project(GraphLab)

# add_graphlab_executable(hello_world hello_world.cpp)
add_graphlab_executable(control-tower templates.cpp pilot.cpp graphlabjs.cpp)
requires_v8(control-tower)
