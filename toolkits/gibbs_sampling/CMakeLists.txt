project(GraphLab)




# add_library(pgibbs STATIC
#     image.cpp
#     chromatic_sampler.cpp 
#     util.cpp 
#     factorized_model.cpp
#     mrf.cpp
#     junction_tree.cpp
#     jt_splash_sampler.cpp
#     pgibbs_tls.cpp
#     global_variables.cpp
#     )

# add_library(pgibbs_pic STATIC
#     image.cpp
#     chromatic_sampler.cpp 
#     util.cpp 
#     factorized_model.cpp
#     mrf.cpp
#     junction_tree.cpp
#     jt_splash_sampler.cpp
#     pgibbs_tls.cpp
#     global_variables.cpp
#     )

# target_link_libraries(pgibbs_pic 
#   graphlab_pic)

# target_link_libraries(pgibbs 
#   graphlab)

# get_property(pgibbs_flags TARGET pgibbs PROPERTY COMPILE_FLAGS)
# set_target_properties(pgibbs_pic
#   PROPERTIES COMPILE_FLAGS "${pgibbs_flags} -fPIC")




# add_graphlab_executable(make_denoise_alchemy 
#   make_denoise_alchemy.cpp 
#   )
# target_link_libraries(make_denoise_alchemy pgibbs)





# add_graphlab_executable(sampler sampler.cpp)
# target_link_libraries(sampler pgibbs)

