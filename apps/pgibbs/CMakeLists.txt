project(GraphLab)

# subdirs(jt_factor factor mrf)

add_executable(make_denoise_alchemy 
  make_denoise_alchemy.cpp 
  image.cpp
  factorized_model.cpp
  )



add_executable(pgibbs 
  pgibbs.cpp
  image.cpp
  chromatic_sampler.cpp 
  util.cpp 
  mrf.cpp
  factorized_model.cpp
  )



