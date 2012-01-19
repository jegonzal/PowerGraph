project(GraphLab)
# NOTE: do not link tcmalloc! Does not like Java.
add_jni_library(graphlabjni org_graphlab_Updater.cpp org_graphlab_Context.cpp org_graphlab_Core.cpp)
