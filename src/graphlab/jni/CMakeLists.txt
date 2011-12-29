project(GraphLab)
# NOTE: do not link tcmalloc! Does not like Java.
add_jni_library(graphlabjni graphlabjni.cpp)
