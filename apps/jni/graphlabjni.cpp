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

static JavaVM *jvm = NULL;
JNIEnv * cachedEnv;
jobject cachedObj;
jmethodID wrapperMethodID;
std::vector<JNIEnv *> envs;

/**
 * The Page rank update function
 */
void jni_update(gl_types::iscope &scope,
                     gl_types::icallback &scheduler,
                     gl_types::ishared_data* shared_data) {
        jint vertexid = scope.vertex();
        jint functionid = 0;
        
        int threadid = graphlab::thread::thread_id();
        JNIEnv * jenv;
        if (envs[threadid] == NULL) {
            int res = jvm->AttachCurrentThread((void **)&jenv, NULL);
            std::cout << "Attached to JVM: " << res << std::endl;
            assert(res>=0);
            envs[threadid] = jenv;
        }
        jenv = envs[threadid];
        
        jintArray result = (jintArray) jenv->CallObjectMethod(cachedObj, wrapperMethodID, vertexid, functionid);
        jsize task_sz = jenv->GetArrayLength(result)/2;
      //  std::cout << "New tasks: " << task_sz << std::endl; 
        jboolean isCopy = false;
        jint * arr = jenv->GetIntArrayElements(result, &isCopy);
        for(int i=0; i<task_sz; i++) {
            // TODO: support multiple tasks, priority
            scheduler.add_task(gl_types::update_task(arr[i], jni_update), 1.0);
        }
        jenv->ReleaseIntArrayElements(result, arr, JNI_ABORT);
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

JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setScheduler
  (JNIEnv * env, jobject obj, jstring schedulertype) {
    const char *str = env->GetStringUTFChars(schedulertype, 0);
    core.set_scheduler_type(std::string(str));
  }

/*
 * Class:     graphlab_wrapper_GraphLabJNIWrapper
 * Method:    setScopeType
 * Signature: (Ljava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setScopeType
  (JNIEnv * env, jobject obj, jstring scopetype) {
    const char *str = env->GetStringUTFChars(scopetype, 0);
    core.set_scope_type(std::string(str));
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
    }
    env->ReleaseIntArrayElements(toVertices, arr, JNI_ABORT);
}
    

/*
 * Class:     graphlab_wrapper_GraphLabJNIWrapper
 * Method:    runGraphlab
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_runGraphlab
  (JNIEnv * env, jobject obj) {
    // HACK - where to get number of workers?
    envs = std::vector<JNIEnv * >(512, NULL);
    if (jvm == NULL) env->GetJavaVM(&jvm);

    cachedObj = obj;
    jclass clazz = env->GetObjectClass(obj);
    wrapperMethodID = env->GetMethodID(clazz, "execUpdate", "(II)[I");
    std::cout << "Got method ID " << wrapperMethodID << std::endl;
    std::cout << "Graph has: " << core.graph().num_vertices() << " vertices and " << 
                core.graph().num_edges() << " edges." << std::endl;
    core.get_engine_options().print();
    double runtime = core.start(); 
    std::cout << "Runtime: " << runtime << " seconds." << std::endl;
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