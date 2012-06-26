project(ComputerVision)

add_graphlab_executable(grabcut grabcut.cpp)
#requires_eigen(imagestitch)
requires_opencv(grabcut)

