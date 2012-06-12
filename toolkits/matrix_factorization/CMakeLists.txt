project(GraphLab)
# include(CheckCXXSourceCompiles)



# Build als
add_graphlab_executable(als 
  als.cpp  eigen_serialization.cpp)
requires_eigen(als) # build and attach eigen


add_graphlab_executable(make_synthetic_als_data
  make_synthetic_als_data.cpp)
requires_eigen(make_synthetic_als_data) # build and attach eigen


