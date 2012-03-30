project(GraphLab)
# NOTE: do not link tcmalloc! Does not like Java.
add_jni_library(graphlabjni
  org_graphlab_Updater.cpp
  org_graphlab_Aggregator.cpp
  org_graphlab_Context.cpp
  org_graphlab_Core.cpp
  java_any.cpp)

# copy to ${SRC}/extapis/java_jni/lib when done
get_target_property(graphlabjni_location graphlabjni LOCATION)
add_custom_command(
  TARGET graphlabjni
  POST_BUILD
  COMMAND cp ${graphlabjni_location}
          ${CMAKE_SOURCE_DIR}/extapis/java_jni/lib/
  COMMAND cd ${CMAKE_SOURCE_DIR}/extapis/java_jni && ant
  COMMENT "Copying libgraphlabjni to extapis/java_jni/lib")