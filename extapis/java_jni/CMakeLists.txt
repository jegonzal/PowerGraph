project(GraphLab)

if(JNI_REALLY_FOUND AND ANT_FOUND)
  add_custom_target(java_app
    COMMAND ant dist
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/my-app/src
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/my-app/bin
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/my-app/lib
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/my-app/docs
    COMMAND cp lib/* ${CMAKE_CURRENT_BINARY_DIR}/my-app/lib
    COMMAND cp dist/lib/graphlab.jar ${CMAKE_CURRENT_BINARY_DIR}/my-app/lib
    COMMAND cp -R dist/docs ${CMAKE_CURRENT_BINARY_DIR}/my-app/
    COMMAND cp -R eclipse/ ${CMAKE_CURRENT_BINARY_DIR}/my-app/
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating new Eclipse project")
  add_dependencies(java_app graphlabjni)
else()
  message(STATUS
    "Ant not detected. You will not be able to generate Eclipse project files from CMake.") 
endif()

