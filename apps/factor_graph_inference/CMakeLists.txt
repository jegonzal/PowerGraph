project(GraphLab)


add_executable(alchemy_to_bp_parser alchemy_to_bp_parser.cpp)

add_executable(factor_bp factor_bp.cpp)
TARGET_LINK_LIBRARIES(factor_bp tcmalloc)

add_executable(factor_gibbs factor_gibbs.cpp)
TARGET_LINK_LIBRARIES(factor_gibbs tcmalloc)

add_executable(belief_diff belief_diff.cpp)
