project(GraphLab)

# add executable
add_graphlab_executable(dispatcher rapidjson.cpp dispatcher.cpp json_message.cpp process.cpp)
# test files
copy_files(*.tsv)
copy_files(*.rb)

subdirs(demoapps)