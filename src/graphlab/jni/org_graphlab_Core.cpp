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
 * @file org_graphlab_Core.cpp
 *
 * Contains the JNI interface for org.graphlab.Core. In general, applications
 * will keep their graphs in the Java layer and access the engine through the
 * JNI. This wrapper provides a proxy graph for the engine to manipulate and
 * forwards update calls to the Java layer. To learn how to use this interface,
 * refer to the org.graphlab.Core class and to the examples.
 *
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */

#include <wordexp.h>
#include "org_graphlab_Core.hpp"
#include "org_graphlab_Updater.hpp"

using namespace graphlab;

//---------------------------------------------------------------
// jni_core static members
//---------------------------------------------------------------

template<typename G, typename U>
JavaVM* jni_core<G, U>::mjvm = NULL;

template<typename G, typename U>
const size_t jni_core<G, U>::ENV_ID = 1;

#ifdef __cplusplus
extern "C" {
#endif

//---------------------------------------------------------------
// static helper functions
//---------------------------------------------------------------

  static jlong createCore (JNIEnv *env, jobject obj, int argc, char **argv){
  
    // configure log level (TODO: allow config)
    global_logger().set_log_level(LOG_DEBUG);
    global_logger().set_log_to_console(true);

    // set jvm, if we don't have it already
    if (NULL == proxy_updater::core::get_jvm()){
      JavaVM* jvm = NULL;
      env->GetJavaVM(&jvm);
      proxy_updater::core::set_jvm(jvm);
    }
    
    // get the method ID for Updater#execUpdate, if we don't have it already
    if (0 == proxy_updater::java_method_id){
      jclass updater_class = env->FindClass("org/graphlab/Updater");
      proxy_updater::java_method_id =
        env->GetMethodID(updater_class, "execUpdate", "(JI)V");
    }
    
    // store env for this thread
    thread::get_local(proxy_updater::core::ENV_ID) = env;
    
    // allocate and configure core
    proxy_updater::core *jni_core = new proxy_updater::core(env, obj);
    if (NULL != argv){
      (*jni_core)().parse_options(argc, argv);
    }
    
    // return address of jni_core
    return (long) jni_core;
  
  }
  
//---------------------------------------------------------------
// JNI functions
//---------------------------------------------------------------

  JNIEXPORT jlong JNICALL
  Java_org_graphlab_Core_createCore
  (JNIEnv *env, jobject obj){
    return createCore(env, obj, 0, NULL);
  }
  
  JNIEXPORT jlong JNICALL Java_org_graphlab_Core_createCore__Ljava_lang_String_2
  (JNIEnv *env, jobject obj, jstring command_line_args){
  
    // convert jstring to c string
    const char *cstr = NULL;
    cstr = env->GetStringUTFChars(command_line_args, NULL);
    if (NULL == cstr) {
       return 0; /* OutOfMemoryError already thrown */
    }
    
    // prepend with dummy name
    char buffer[1024];
    snprintf(buffer, 1024, "x %s", cstr);
    env->ReleaseStringUTFChars(command_line_args, cstr);
    
    // split string
    wordexp_t we_result;
    if (0 != wordexp(buffer, &we_result, 0)) return 0;
   
    // create core
    jlong ptr = createCore(env, obj, we_result.we_wordc, we_result.we_wordv);
   
    // cleanup
    wordfree(&we_result);
    return ptr;
  
  }
  
  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_destroyCore
  (JNIEnv *env, jobject obj, jlong ptr){
    
    if (NULL == env || 0 == ptr){
      proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
      return;
    }
    
    // cleanup
    proxy_updater::core *jni_core = (proxy_updater::core *) ptr;
    delete jni_core;
    
    logstream(LOG_DEBUG)
      << "GraphLab core deleted in JNI."
      << std::endl;
    
  }

  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_resizeGraph
  (JNIEnv *env, jobject obj, jlong ptr, jint count){
    
    if (NULL == env || 0 == ptr){
      proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
      return;
    }
    
    proxy_updater::core *jni_core = (proxy_updater::core *) ptr;
    (*jni_core)().graph().resize(count);
    
  }
  
  JNIEXPORT jint JNICALL
  Java_org_graphlab_Core_addVertex
  (JNIEnv *env, jobject obj, jlong ptr, jint id){
  
    if (NULL == env || 0 == ptr){
      proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
      return -1;
    }
    
    proxy_updater::core *jni_core = (proxy_updater::core *) ptr;
    
    // init vertex
    proxy_vertex vertex;
    vertex.app_id = id;
    
    // add to graph
    return (*jni_core)().graph().add_vertex(vertex);
  
  }
  
  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_addEdge
  (JNIEnv *env, jobject obj, jlong ptr, jint source, jint target){
  
    if (NULL == env || 0 == ptr){
      proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }
    
    proxy_updater::core *jni_core = (proxy_updater::core *) ptr;
    
    // add to graph
    (*jni_core)().graph().add_edge(source, target, proxy_edge());
  
  }
  
  JNIEXPORT jdouble JNICALL
  Java_org_graphlab_Core_start
  (JNIEnv *env, jobject obj, jlong ptr){
    
    if (NULL == env || 0 == ptr){
      proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return 0;
    }
    
    proxy_updater::core *jni_core = (proxy_updater::core *) ptr;

    logstream(LOG_DEBUG)
      << "Graph has: "
      << (*jni_core)().graph().num_vertices() << " vertices and "
      << (*jni_core)().graph().num_edges() << " edges."
      << std::endl;
    (*jni_core)().engine().get_options().print();
    
    double runtime = (*jni_core)().start(); 

    logstream(LOG_INFO)
        << "Finished after " 
	      << (*jni_core)().engine().last_update_count() << " updates."
	      << std::endl;
    logstream(LOG_INFO)
        << "Runtime: " << runtime
	      << " seconds."
	      << std::endl;
    
    return runtime;
    
  }

  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_setNCpus
  (JNIEnv * env, jobject obj, jlong ptr, jlong ncpus) {
  
    if (NULL == env || 0 == ptr){
      proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }
  
    proxy_updater::core *jni_core = (proxy_updater::core *) ptr;
    (*jni_core)().set_ncpus(ncpus);
    
  }

  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_setSchedulerType
  (JNIEnv * env, jobject obj, jlong ptr, jstring scheduler_str) {
  
    if (NULL == env || 0 == ptr){
      proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }
  
    const char *str = env->GetStringUTFChars(scheduler_str, NULL);
    if (NULL == str) return;  // OutOfMemoryError already thrown
    
    proxy_updater::core *jni_core = (proxy_updater::core *) ptr;
    (*jni_core)().set_scheduler_type(std::string(str));
    env->ReleaseStringUTFChars(scheduler_str, str);
    
  }

  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_setScopeType
  (JNIEnv * env, jobject obj, jlong ptr, jstring scope_str) {
  
    if (NULL == env || 0 == ptr){
      proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "ptr must not be null.");
        return;
    }
  
    const char *str = env->GetStringUTFChars(scope_str, NULL);
    if (NULL == str) return;  // OutOfMemoryError already thrown
    
    proxy_updater::core *jni_core = (proxy_updater::core *) ptr;
    (*jni_core)().set_scope_type(std::string(str));
    env->ReleaseStringUTFChars(scope_str, str);
    
  }
  
  JNIEXPORT void JNICALL
  Java_org_graphlab_Core_schedule
  (JNIEnv * env, jobject obj,
  jlong core_ptr, jobject updater, jint vertex_id){
  
    if (NULL == env || 0 == core_ptr){
      proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "core_ptr must not be null.");
        return;
    }

    // get objects from pointers
    proxy_updater::core *jni_core = (proxy_updater::core *) core_ptr;

    // schedule vertex
    (*jni_core)().schedule(vertex_id, proxy_updater(env, updater));
    
  }
  
  JNIEXPORT void JNICALL 
  Java_org_graphlab_Core_scheduleAll
  (JNIEnv * env, jobject obj,
  jlong core_ptr, jobject updater) {

    if (NULL == env || 0 == core_ptr){
    proxy_updater::core::throw_exception(
        env,
        "java/lang/IllegalArgumentException",
        "core_ptr and updater_ptr must not be null.");
        return;
    }

    // get objects from pointers
    proxy_updater::core *jni_core = (proxy_updater::core *) core_ptr;

    // schedule vertex
    (*jni_core)().schedule_all(proxy_updater(env, updater));

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

#ifdef __cplusplus
}
#endif