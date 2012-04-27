project(GraphLab)
 link_libraries(boost_date_time)
# Files were not added
# add_graphlab_executable(stringparser stringparser.cpp formatconverter.cpp)
add_graphlab_executable(stringzipparser stringzipparser.cpp formatconverter.cpp ../shared/mmio.cpp)
add_graphlab_executable(texttokenparser texttokenparser.cpp ../shared/mmio.cpp)
add_graphlab_executable(undirectedgraph undirectedgraph.cpp formatconverter.cpp ../shared/mmio.cpp)
#add_graphlab_executable(undirectedgraph2 undirectedgraph2.cpp formatconverter.cpp ../shared/mmio.cpp)
add_graphlab_executable(graphbuilder graphbuilder.cpp
formatconverter.cpp ../shared/mmio.cpp)
add_graphlab_executable(matmat2bin matmat2bin.cpp ../shared/mmio.cpp)
add_graphlab_executable(matmat2bin2 matmat2bin2.cpp ../shared/mmio.cpp)
add_graphlab_executable(reversemap reversemap.cpp ../shared/mmio.cpp)
add_graphlab_executable(kddcup_output_builder kddcup_output_builder.cpp ../shared/mmio.cpp)
add_graphlab_executable(kdd_test_splitter kdd_test_splitter.cpp ../shared/mmio.cpp)
add_graphlab_executable(kdd_train_data_parser kdd_train_data_parser.cpp ../shared/mmio.cpp formatconverter.cpp)
add_graphlab_executable(kdd_usr_itm_feature_parser kdd_usr_itm_feature_parser.cpp ../shared/mmio.cpp)
add_graphlab_executable(kdd_linear_model_builder kdd_linear_model_builder.cpp
../shared/mmio.cpp)
add_graphlab_executable(kdd_track2_parser kdd_track2_parser.cpp
../shared/mmio.cpp)
add_graphlab_executable(kdd_track2_feature_parser
kdd_track2_feature_parser.cpp ../shared/mmio.cpp)
