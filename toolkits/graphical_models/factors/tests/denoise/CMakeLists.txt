project(GraphLab)

add_graphlab_executable(denoise denoise.cpp)
requires_opencv(denoise)

file(COPY noisy_img.png 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY denoised_gm.png
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

