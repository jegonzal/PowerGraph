project(ComputerVision)

add_graphlab_executable(grabcut grabcut.cpp)
requires_opencv(grabcut)

add_graphlab_executable(stitching stitching.cpp)
requires_opencv(stitching)

add_graphlab_executable(stitching_detailed stitching_detailed.cpp)
requires_opencv(stitching_detailed)

add_graphlab_executable(stitch stitch_main.cpp)
requires_opencv(stitch)
requires_eigen(stitch)

add_graphlab_executable(stitch_full stitch_full_main.cpp)
requires_opencv(stitch_full)
requires_eigen(stitch_full)
