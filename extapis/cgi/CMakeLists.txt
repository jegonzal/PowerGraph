project(GraphLab)

# add executable
add_graphlab_executable(dispatcher rapidjson.cpp json_message.cpp process.cpp dispatcher.cpp)
# test files
copy_files(*.tsv)

subdirs(demoapps)