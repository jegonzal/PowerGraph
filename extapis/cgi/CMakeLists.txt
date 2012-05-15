project(GraphLab)

# add executable
add_graphlab_executable(dispatcher rapidjson.cpp json_message.cpp process.cpp dispatcher.cpp)
add_dependencies(dispatcher jsonsl)
target_link_libraries(dispatcher ${JSONSL_LIB_DIR}/libjsonsl.a)
# test files
copy_files(*.tsv)