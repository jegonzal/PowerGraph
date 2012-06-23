project(GraphLab)

# add_graphlab_executable(hello_world hello_world.cpp)
add_graphlab_executable(graphlabjs graphlabjs.cpp)
requires_v8(graphlabjs)
