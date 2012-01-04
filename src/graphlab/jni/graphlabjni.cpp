/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */

/**
 * @file graphlabjni.cpp
 *
 * Contains the JNI interface between Java and C++. In general, applications
 * will keep their graphs in the Java layer and access the engine through the
 * JNI. This wrapper provides a proxy graph for the engine to manipulate and
 * forwards update calls to the Java layer. To learn how to use this interface,
 * refer to the org.graphlab.Core class and to the examples.
 *
 * Most of the methods in this file belong to the org.graphlab.Core class.
 * @author akyrola
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */

#include <string.h>

#include <cctype>
#include <string>

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

#include "jni_core.hpp"
#include "org_graphlab_Core.h"
#include "org_graphlab_Context.h"

using namespace graphlab;

/** Proxy edge */
struct edge_data {};

/** Proxy vertex */
struct vertex_data {
  int app_id;   /**< corresponding application vertex ID */
};

/** Proxy graph */
typedef graph<vertex_data, edge_data> graph_type;

/**
 * Proxy updater
 * Forwards update calls to the Java layer
 */
class proxy_updater : 
  public iupdate_functor<graph_type, proxy_updater> {
  
public:

  static jmethodID java_method_id;  // id of Core#execUpdate
  
  jobject obj;    // java core object
  jint id;        // ID of java updater object
  
  void operator()(icontext_type& context) {

    JNIEnv *jenv = jni_core<graph_type, proxy_updater>::get_JNIEnv ();
    jint app_vertex_id = context.vertex_data().app_id;
    jenv->CallVoidMethod (obj, java_method_id,
                         &context,
                         app_vertex_id,
                         id);
    
    // check for exception
    jthrowable exc = jenv->ExceptionOccurred();
    if (exc) {
      logstream(LOG_ERROR)
        << "Exception occured!!"
        << std::endl;
        
      jclass new_exc;
      jenv->ExceptionDescribe();
      jenv->ExceptionClear();
      new_exc = jenv->FindClass("java/lang/IllegalArgumentException");
      if (new_exc == NULL) return;
      jenv->ThrowNew(new_exc, "thrown from C code");
      
    }

  }
  
};
jmethodID proxy_updater::java_method_id = 0;

typedef jni_core<graph_type, proxy_updater> jni_core_type;
typedef iupdate_functor<graph_type, proxy_updater>::icontext_type icontext_type;

#ifdef __cplusplus
extern "C" {
#endif

  JNIEXPORT jlong JNICALL
  Java_org_graphlab_Core_createCore
  (JNIEnv *env, jobject obj){
  
    // TODO: allow config of log options and command line options
    // TODO: error and exception handling
  
    // configure log level (TODO: allow config)
    global_logger().set_log_level(LOG_DEBUG);
    global_logger().set_log_to_console(true);
 
    // setup the parser
    command_line_options clopts("JNI options. TODO.");

    // set jvm, if we don't have it already
    if (NULL == jni_core_type::get_jvm()){
      JavaVM* jvm = NULL;
      env->GetJavaVM(&jvm);
      jni_core_type::set_jvm(jvm);
    }
    
    // get the method ID, if we don't have it already
    jclass clazz = env->GetObjectClass(obj);
    if (0 == proxy_updater::java_method_id){
      proxy_updater::java_method_id =
        env->GetMethodID(clazz, "execUpdate", "(JII)V");
    }
    
    // allocate and configure core
    jni_core_type *jni_core = new jni_core_type(env, obj);
    jni_core->core().set_options(clopts);
    
    logstream(LOG_DEBUG)
      << "GraphLab core initialized in JNI."
      << std::endl;
    
    // return address of core
    return (long) jni_core;
    
  }
  
  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_destroyCore
  (JNIEnv *env, jobject obj, jlong ptr){
    
    if (NULL == env || 0 == ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
      return;
    }
    
    // cleanup
    jni_core_type *jni_core = (jni_core_type *) ptr;
    delete jni_core;
    
    logstream(LOG_DEBUG)
      << "GraphLab core deleted in JNI."
      << std::endl;
    
  }

  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_resizeGraph
  (JNIEnv *env, jobject obj, jlong ptr, jint count){
    
    if (NULL == env || 0 == ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
      return;
    }
    
    jni_core_type *jni_core = (jni_core_type *) ptr;
    jni_core->core().graph().resize(count);
    
  }
  
  JNIEXPORT jint JNICALL
  Java_org_graphlab_Core_addVertex
  (JNIEnv *env, jobject obj, jlong ptr, jint id){
  
    if (NULL == env || 0 == ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
      return -1;
    }
    
    jni_core_type *jni_core = (jni_core_type *) ptr;
    
    // init vertex
    vertex_data vertex;
    vertex.app_id = id;
    
    // add to graph
    return jni_core->core().graph().add_vertex(vertex);
  
  }
  
  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_addEdge
  (JNIEnv *env, jobject obj, jlong ptr, jint source, jint target){
  
    if (NULL == env || 0 == ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }
    
    jni_core_type *jni_core = (jni_core_type *) ptr;
    
    // add to graph
    jni_core->core().graph().add_edge(source, target, edge_data());
  
  }
  
  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_schedule
  (JNIEnv * env, jobject obj, jlong ptr, jint vertex_id, jint updater_id){
  
    if (NULL == env || 0 == ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }

    jni_core_type *jni_core = (jni_core_type *) ptr;
  
    // initialize proxy updater
    proxy_updater updater;
    updater.obj = jni_core->obj();
    updater.id = updater_id;

    // schedule vertex
    jni_core->core().schedule(vertex_id, updater);
    
  }
  
  JNIEXPORT jdouble JNICALL
  Java_org_graphlab_Core_start
  (JNIEnv *env, jobject obj, jlong ptr){
    
    if (NULL == env || 0 == ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return 0;
    }
    
    jni_core_type *jni_core = (jni_core_type *) ptr;

    logstream(LOG_DEBUG)
      << "Graph has: "
      << jni_core->core().graph().num_vertices() << " vertices and "
      << jni_core->core().graph().num_edges() << " edges."
      << std::endl;
    jni_core->core().engine().get_options().print();
    
//     if (taskbudget>0)
//       core.engine().set_task_budget(taskbudget);
//     if (maxiter>0)
//       core.sched_options().add_option("max_iterations", maxiter);

    // set thread destroy callback
    thread::set_thread_destroy_callback(jni_core_type::detach_from_jvm);
    double runtime = jni_core->core().start(); 

//     if (metrics_type != "none") {
//       core.set_metrics_type(metrics_type);
//       core.fill_metrics();
//       core.report_metrics();
//       // Hack: this prevents core destructor from dumping the metrics.
//       // ... which leads to some weird mutex error.
//       core.set_metrics_type("none");
//     }

    logstream(LOG_INFO)
        << "Finished after " 
	      << jni_core->core().engine().last_update_count() << " updates."
	      << std::endl;
    logstream(LOG_INFO)
        << "Runtime: " << runtime 
	      << " seconds."
	      << std::endl;
    
    return runtime;
    
  }
  
  JNIEXPORT void JNICALL
  Java_org_graphlab_Context_schedule
  (JNIEnv *env, jobject obj, jlong core_ptr, jlong context_ptr, jint vertex_id, jint updater_id){
    
    if (NULL == env || 0 == core_ptr || 0 == context_ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }

    jni_core_type *jni_core = (jni_core_type *) core_ptr;
    icontext_type *context = (icontext_type *) context_ptr;
  
    // initialize proxy updater
    proxy_updater updater;
    updater.obj = jni_core->obj();
    updater.id = updater_id;
    
    context->schedule(vertex_id, updater);
    
  }

  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_setNCpus
  (JNIEnv * env, jobject obj, jlong ptr, jlong ncpus) {
  
    if (NULL == env || 0 == ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }
  
    jni_core_type *jni_core = (jni_core_type *) ptr;
    jni_core->core().set_ncpus(ncpus);
    
  }

  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_setSchedulerType
  (JNIEnv * env, jobject obj, jlong ptr, jstring scheduler_str) {
  
    if (NULL == env || 0 == ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }
  
    const char *str = env->GetStringUTFChars(scheduler_str, NULL);
    if (NULL == str) return;  // OutOfMemoryError already thrown
    
    jni_core_type *jni_core = (jni_core_type *) ptr;
    jni_core->core().set_scheduler_type(std::string(str));
    env->ReleaseStringUTFChars(scheduler_str, str);
    
  }

  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_setScopeType
  (JNIEnv * env, jobject obj, jlong ptr, jstring scope_str) {
  
    if (NULL == env || 0 == ptr){
      jni_core_type::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }
  
    const char *str = env->GetStringUTFChars(scope_str, NULL);
    if (NULL == str) return;  // OutOfMemoryError already thrown
    
    jni_core_type *jni_core = (jni_core_type *) ptr;
    jni_core->core().set_scope_type(std::string(str));
    env->ReleaseStringUTFChars(scope_str, str);
    
  }
 
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setVertexColors
//   (JNIEnv * env, jobject obj, jintArray colors) {
//     jni_graph & graph = core.graph();
//     jsize sz = env->GetArrayLength(colors);
//     jboolean isCopy = false;
//     jint * arr = env->GetIntArrayElements(colors, &isCopy);
//     for(int i=0; i<sz; i++) {
//       graph.color(i) = gl_types::vertex_color(arr[i]);
//     }
//     env->ReleaseIntArrayElements(colors, arr, JNI_ABORT);
//   }
// 
// 
//   /*
//    * Class:     graphlab_wrapper_GraphLabJNIWrapper
//    * Method:    schedule
//    * Signature: ([I[I)V
//    */
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_schedule
//   (JNIEnv * env, jobject obj, jintArray vertices, jintArray funcs) {
//     // TODO: support multiple update functions
//     jboolean isCopy = false;
//     jsize sz = env->GetArrayLength(vertices);
//     jint * arr = env->GetIntArrayElements(vertices, &isCopy);
//     jint * funcarr = env->GetIntArrayElements(funcs, &isCopy);
//     for(int i=0; i<sz; i++) {
//       core.add_task(gl_types::update_task(arr[i], functions[funcarr[i]]), 1.0);
//     }
//     env->ReleaseIntArrayElements(vertices, arr, JNI_ABORT);
//     env->ReleaseIntArrayElements(funcs, funcarr, JNI_ABORT);
//   }
// 
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setTaskBudget
//   (JNIEnv * env, jobject obj, jint budget) {
//     std::cout << "Set task budget: " << budget << std::endl;
//     taskbudget = budget;
//   }
// 
//   /*
//    * Class:     graphlab_wrapper_GraphLabJNIWrapper
//    * Method:    setIterations
//    * Signature: (I)V
//    */
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setIterations
//   (JNIEnv * env, jobject obj, jint iter) {
//     maxiter = iter;
//   }
// 
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setMetrics
//   (JNIEnv * env, jobject obj, jstring schedulertype) {
//     const char *str = env->GetStringUTFChars(schedulertype, 0);
//     metrics_type = std::string(str);
//   }
// 
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_computeGraphColoring
//   (JNIEnv * env, jobject obj, jint ncpus) {
//     core.graph().compute_coloring();
//   }
// 
//   /*
//    * Class:     graphlab_wrapper_GraphLabJNIWrapper
//    * Method:    scheduleAll
//    * Signature: (I)V
//    */
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_scheduleAll
//   (JNIEnv * env, jobject obj, jint funcid) {
//     core.add_task_to_all(functions[0], 1.0);
//   }

  JNIEXPORT jint JNICALL Java_graphlab_test_JniTest_dummy
  (JNIEnv * env, jobject obj, jint i) {
    return i+1;
  }

#ifdef __cplusplus
}
#endif

#include <graphlab/macros_undef.hpp>