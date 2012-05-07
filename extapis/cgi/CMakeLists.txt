project(GraphLab)

# add executable
add_graphlab_executable(dispatcher process.cpp dispatcher.cpp json_message.cpp)
target_link_libraries(dispatcher ${JSONSL_LIB_DIR}/libjsonsl.a)
# test files
copy_files(*.tsv)