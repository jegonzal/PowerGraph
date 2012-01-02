project(GraphLab)
 link_libraries(boost_date_time)
# Files were not added
# add_graphlab_executable(stringparser stringparser.cpp formatconverter.cpp)
add_graphlab_executable(stringzipparser stringzipparser.cpp
formatconverter.cpp ../shared/mmio.cpp)
add_graphlab_executable(undirectedgraph undirectedgraph.cpp
formatconverter.cpp ../shared/mmio.cpp)
add_graphlab_executable(graphbuilder graphbuilder.cpp
formatconverter.cpp ../shared/mmio.cpp)
