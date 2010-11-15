#include <jni.h>
#include <string.h>
#include <cctype>
 

#include <string>

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>


// Just dummy edgedata and vertexdata
struct edge_data {}; 
struct vertex_data {};


//! The type of graph used in this program
typedef graphlab::graph<vertex_data, edge_data> jni_graph;


/**
 * The collection of graphlab types restricted to the graph type used
 * in this program.
 */
typedef graphlab::types<jni_graph> gl_types;

/**
 * The Page rank update function
 */
void jni_update(gl_types::iscope &scope,
                     gl_types::icallback &scheduler,
                     gl_types::ishared_data* shared_data) {

        std::cout << "Jni update: " << scope.vertex() << std::endl;
}  


 


 #ifdef __cplusplus
extern "C" {
#endif
 
gl_types::core core;
 
 
/*
 * Class:     graphlab_wrapper_GraphLabJNIWrapper
 * Method:    initGraphlab
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_initGraphlab
  (JNIEnv * env, jobject obj) {
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
    logger(LOG_INFO, "JNI Graphlab starting.");
 
    // Setup the parser
    graphlab::command_line_options
      clopts("JNI options. TODO.");

   // Set the engine options
   // TODO.
   core.set_engine_options(clopts);
  
 }

/*
 * Class:     graphlab_wrapper_GraphLabJNIWrapper
 * Method:    createGraph
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_createGraph
  (JNIEnv * env, jobject obj, jint numvertices) {
       // Create the dummy vertices
       jni_graph & graph = core.graph();
       for(int i=0; i<numvertices; i++) {
           graph.add_vertex(vertex_data());
       }
  }

/*
 * Class:     graphlab_wrapper_GraphLabJNIWrapper
 * Method:    addEdges
 * Signature: (I[I)V
 */
JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_addEdges
    (JNIEnv * env, jobject obj, jint fromVertex, jintArray toVertices) {
    jni_graph & graph = core.graph();
    jsize sz = env->GetArrayLength(toVertices);
    jboolean isCopy = false;
    jint * arr = env->GetIntArrayElements(toVertices, &isCopy);
    for(int i=0; i<sz; i++) {
        graph.add_edge(fromVertex, arr[i]);
        std::cout << " Added edge: " << fromVertex << " to " << arr[i] << std::endl;
    }
    //ReleaseIntArrayElements(env, toVertices, arr, JNI_ABORT);
}
    

/*
 * Class:     graphlab_wrapper_GraphLabJNIWrapper
 * Method:    runGraphlab
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_runGraphlab
  (JNIEnv * env, jobject obj) {
    core.start(); 
}

/*
 * Class:     graphlab_wrapper_GraphLabJNIWrapper
 * Method:    schedule
 * Signature: ([I[I)V
 */
JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_schedule
  (JNIEnv * env, jobject obj, jintArray vertices, jintArray funcs) {
  // TODO: support multiple update functions
    jboolean isCopy = false;
    jsize sz = env->GetArrayLength(vertices);
    jint * arr = env->GetIntArrayElements(vertices, &isCopy);
    for(int i=0; i<sz; i++) {
        core.add_task(gl_types::update_task(arr[i], jni_update), 1.0);
    }
 }

/*
 * Class:     graphlab_wrapper_GraphLabJNIWrapper
 * Method:    scheduleAll
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_scheduleAll
  (JNIEnv * env, jobject obj, jint funcid) {
        core.add_task_to_all(jni_update, 1.0);
  }

JNIEXPORT jint JNICALL Java_graphlab_test_JniTest_dummy
  (JNIEnv * env, jobject obj, jint i) {
    return i+1;
}

#ifdef __cplusplus
}
#endif