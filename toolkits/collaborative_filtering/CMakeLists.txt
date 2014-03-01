project(GraphLab)
# include(CheckCXXSourceCompiles)



# Build als
add_graphlab_executable(als als.cpp)
requires_eigen(als) # build and attach eigen

add_graphlab_executable(sparse_als sparse_als.cpp)
requires_eigen(sparse_als) # build and attach eigen

add_graphlab_executable(wals wals.cpp)
requires_eigen(wals) # build and attach eigen

add_graphlab_executable(sgd sgd.cpp)
requires_eigen(sgd) # build and attach eigen

add_graphlab_executable(biassgd biassgd.cpp)
requires_eigen(biassgd) # build and attach eigen

add_graphlab_executable(svdpp svdpp.cpp)
requires_eigen(svdpp) # build and attach eigen

add_graphlab_executable(svd svd.cpp)
requires_eigen(svd) # build and attach eigen


add_graphlab_executable(nmf nmf.cpp)
requires_eigen(nmf) # build and attach eigen

add_graphlab_executable(make_synthetic_als_data
  make_synthetic_als_data.cpp)
requires_eigen(make_synthetic_als_data) # build and attach eigen

add_graphlab_executable(adpredictor adpredictor.cpp)

# add_graphlab_executable(warp_nmf warp_nmf.cpp)
# requires_eigen(warp_sgd) # build and attach eigen
# 
# add_graphlab_executable(warp_als_coord warp_als_coord.cpp)
# requires_eigen(warp_als_coord) # build and attach eigen
