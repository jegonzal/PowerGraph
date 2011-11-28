project(GraphLab)
link_libraries(boost_date_time)
add_graphlab_executable(stringparser stringparser.cpp formatconverter.cpp ../shared/mmio.cpp)
