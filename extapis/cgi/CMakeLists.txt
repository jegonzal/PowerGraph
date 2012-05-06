project(GraphLab)

  # add executable
  add_graphlab_executable(dispatcher process.cpp dispatcher.cpp json_message.cpp)
  # test files
  copy_files(*.tsv)