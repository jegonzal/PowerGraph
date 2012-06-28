project(GraphLab)

add_graphlab_executable(control-tower templates.cpp pilot.cpp graphlabjs.cpp)
requires_v8(control-tower)

copy_file(test.js)
