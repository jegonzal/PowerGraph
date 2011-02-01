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


int taskbudget=0;
int maxiter=0;
std::string metrics_type = "none";
 
void detach_from_JVM() {
  int threadid = graphlab::thread::thread_id();
  if (envs[threadid] != NULL) {
    int res = jvm->DetachCurrentThread();
    std::cout << "Detached from JVM: " << res << std::endl;
    assert(res>=0);
  }
}

void jni_update_wrapper(gl_types::iscope &scope,
			gl_types::icallback &scheduler,
			gl_types::ishared_data* shared_data,
			jint functionid);

void jni_update_0(gl_types::iscope &scope,
		  gl_types::icallback &scheduler,
		  gl_types::ishared_data* shared_data) {
  jni_update_wrapper(scope, scheduler, shared_data, 0);                  
}

void jni_update_1(gl_types::iscope &scope,
		  gl_types::icallback &scheduler,
		  gl_types::ishared_data* shared_data) {
  jni_update_wrapper(scope, scheduler, shared_data, 1);                  
}
 
void jni_update_2(gl_types::iscope &scope,
		  gl_types::icallback &scheduler,
		  gl_types::ishared_data* shared_data) {
  jni_update_wrapper(scope, scheduler, shared_data, 2);                  
}
 
void jni_update_3(gl_types::iscope &scope,
		  gl_types::icallback &scheduler,
		  gl_types::ishared_data* shared_data) {
  jni_update_wrapper(scope, scheduler, shared_data, 3);                  
}
 
 
void jni_update_4(gl_types::iscope &scope,
		  gl_types::icallback &scheduler,
		  gl_types::ishared_data* shared_data) {
  jni_update_wrapper(scope, scheduler, shared_data, 4);                  
}

gl_types::update_function functions[5] = {jni_update_0, jni_update_1, jni_update_2, jni_update_3, jni_update_4};

/**
 * The Page rank update function
 */
void jni_update_wrapper(gl_types::iscope &scope,
			gl_types::icallback &scheduler,
			gl_types::ishared_data* shared_data,
			jint functionid) {
  jint vertexid = scope.vertex();        
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
  jsize task_sz = jenv->GetArrayLength(result)/3;
  //  std::cout << "New tasks: " << task_sz << std::endl; 
  jboolean isCopy = false;
  jint * arr = jenv->GetIntArrayElements(result, &isCopy);
        
  // Check for exception
  jthrowable exc = jenv->ExceptionOccurred();
  if (exc) {
    std::cout << "Exception occured!!" << std::endl;
    /* We don't do much with the exception, except that
       we print a debug message for it, clear it, and 
       throw a new exception. */
    jclass newExcCls;
    jenv->ExceptionDescribe();
    jenv->ExceptionClear();
    newExcCls = jenv->FindClass("java/lang/IllegalArgumentException");
    if (newExcCls == NULL) {
      /* Unable to find the exception class, give up. */
      return;
    }
    jenv->ThrowNew(newExcCls, "thrown from C code");
  }
        
  for(int i=0; i<task_sz; i++) {
    // TODO: support multiple tasks, priority
    scheduler.add_task(gl_types::update_task(arr[i], functions[arr[i+task_sz]]), arr[i+2*task_sz]/(1.0e6));
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
    

  JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setVertexColors
  (JNIEnv * env, jobject obj, jintArray colors) {
    jni_graph & graph = core.graph();
    jsize sz = env->GetArrayLength(colors);
    jboolean isCopy = false;
    jint * arr = env->GetIntArrayElements(colors, &isCopy);
    for(int i=0; i<sz; i++) {
      graph.color(i) = (graphlab::vertex_color_type) arr[i];
    }
    env->ReleaseIntArrayElements(colors, arr, JNI_ABORT);
  }

  /*
   * Class:     graphlab_wrapper_GraphLabJNIWrapper
   * Method:    runGraphlab
   * Signature: ()V
   */
  JNIEXPORT void JNICALL 
  Java_graphlab_wrapper_GraphLabJNIWrapper_runGraphlab(JNIEnv * env, jobject obj) {
    // HACK - where to get number of workers?
    envs.clear(); envs.resize(512, NULL);

    if (jvm == NULL) env->GetJavaVM(&jvm);

    cachedObj = obj;
    jclass clazz = env->GetObjectClass(obj);
    wrapperMethodID = env->GetMethodID(clazz, "execUpdate", "(II)[I");
    std::cout << "Got method ID " << wrapperMethodID << std::endl;
    std::cout << "Graph has: " << core.graph().num_vertices() << " vertices and " << 
      core.graph().num_edges() << " edges." << std::endl;
    core.get_engine_options().print();
    if (taskbudget>0)
      core.engine().set_task_budget(taskbudget);
    if (maxiter>0)
      core.sched_options().add_option("max_iterations", maxiter);
 
    // Set thread destroy callback
    graphlab::thread::set_thread_destroy_callback(detach_from_JVM);
   
    double runtime = core.start(); 
    if (metrics_type != "none") {
      core.set_metrics_type(metrics_type);
      core.fill_metrics();
      core.report_metrics();
      // Hack: this prevents core destructor from dumping the metrics.
      // ... which leads to some weird mutex error.
      core.set_metrics_type("none");
    }
    std::cout << "Finished after " << core.engine().last_update_count() << " updates." << std::endl;
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
    jint * funcarr = env->GetIntArrayElements(funcs, &isCopy);
    for(int i=0; i<sz; i++) {
      core.add_task(gl_types::update_task(arr[i], functions[funcarr[i]]), 1.0);
    }
    env->ReleaseIntArrayElements(vertices, arr, JNI_ABORT);
    env->ReleaseIntArrayElements(funcs, funcarr, JNI_ABORT);
  }

  JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setTaskBudget
  (JNIEnv * env, jobject obj, jint budget) {
    std::cout << "Set task budget: " << budget << std::endl;
    taskbudget = budget;
  }

  /*
   * Class:     graphlab_wrapper_GraphLabJNIWrapper
   * Method:    setIterations
   * Signature: (I)V
   */
  JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setIterations
  (JNIEnv * env, jobject obj, jint iter) {
    maxiter = iter;
  }

  JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setMetrics
  (JNIEnv * env, jobject obj, jstring schedulertype) {
    const char *str = env->GetStringUTFChars(schedulertype, 0);
    metrics_type = std::string(str);
  }

  JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setNumCPUs
  (JNIEnv * env, jobject obj, jint ncpus) {
    core.set_ncpus(ncpus);
  }



  JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_computeGraphColoring
  (JNIEnv * env, jobject obj, jint ncpus) {
    core.graph().compute_coloring();
  }



  /*
   * Class:     graphlab_wrapper_GraphLabJNIWrapper
   * Method:    scheduleAll
   * Signature: (I)V
   */
  JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_scheduleAll
  (JNIEnv * env, jobject obj, jint funcid) {
    core.add_task_to_all(functions[0], 1.0);
  }

  JNIEXPORT jint JNICALL Java_graphlab_test_JniTest_dummy
  (JNIEnv * env, jobject obj, jint i) {
    return i+1;
  }

#ifdef __cplusplus
}
#endif
