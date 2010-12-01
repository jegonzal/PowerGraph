project(GraphLab)
link_libraries("")
include_directories(/usr/lib/jvm/java-6-openjdk/include)
# NOTE: do not link tcmalloc! Does not like Java.
add_jni_library(graphlabjni graphlabjni.cpp)
