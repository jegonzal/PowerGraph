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
 * refer to the GraphLabJNIWrapper class and to the examples.
 *
 * Most of the methods in this file belong to the GraphLabJNIWrapper class.
 * @author akyrola
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */

#include <jni.h>
#include <string.h>

#include <cctype>
#include <string>

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

/** Proxy edge */
struct edge_data {};
/** Proxy vertex */
struct vertex_data {};
/** Proxy graph */
typedef graphlab::graph<vertex_data, edge_data> graph_type;

// static JavaVM *jvm = NULL;
// JNIEnv * cachedEnv;
// jobject cachedObj;
// jmethodID wrapperMethodID;
// std::vector<JNIEnv *> envs;
// 
// 
// int taskbudget=0;
// int maxiter=0;
// std::string metrics_type = "none";
//  
// void detach_from_JVM() {
//   int threadid = graphlab::thread::thread_id();
//   if (envs[threadid] != NULL) {
//     int res = jvm->DetachCurrentThread();
//     std::cout << "Detached from JVM: " << res << std::endl;
//     assert(res>=0);
//   }
// }
// 
// void jni_update_wrapper(gl_types::iscope &scope,
// 			gl_types::icallback &scheduler,
// 			jint functionid);
// 
// void jni_update_0(gl_types::iscope &scope,
// 		  gl_types::icallback &scheduler) {
//   jni_update_wrapper(scope, scheduler, 0);                  
// }
// 
// void jni_update_1(gl_types::iscope &scope,
// 		  gl_types::icallback &scheduler) {
//   jni_update_wrapper(scope, scheduler, 1);                  
// }
//  
// void jni_update_2(gl_types::iscope &scope,
// 		  gl_types::icallback &scheduler) {
//   jni_update_wrapper(scope, scheduler, 2);                  
// }
//  
// void jni_update_3(gl_types::iscope &scope,
// 		  gl_types::icallback &scheduler) {
//   jni_update_wrapper(scope, scheduler,  3);                  
// }
//  
//  
// void jni_update_4(gl_types::iscope &scope,
// 		  gl_types::icallback &scheduler) {
//   jni_update_wrapper(scope, scheduler,  4);                  
// }
// 
// gl_types::update_function functions[5] = 
//   {jni_update_0, jni_update_1, jni_update_2, jni_update_3, jni_update_4};
// 
// /**
//  * The Page rank update function
//  */
// void jni_update_wrapper(gl_types::iscope &scope,
// 			gl_types::icallback &scheduler,
// 			jint functionid) {
//   jint vertexid = scope.vertex();        
//   int threadid = graphlab::thread::thread_id();
//   JNIEnv * jenv;
//   if (envs[threadid] == NULL) {
//     int res = jvm->AttachCurrentThread((void **)&jenv, NULL);
//     std::cout << "Attached to JVM: " << res << std::endl;
//     assert(res>=0);
//     envs[threadid] = jenv;
//   }
//   jenv = envs[threadid];
//         
//   jintArray result = (jintArray) jenv->CallObjectMethod(cachedObj, wrapperMethodID, vertexid, functionid);
//   jsize task_sz = jenv->GetArrayLength(result)/3;
//   //  std::cout << "New tasks: " << task_sz << std::endl; 
//   jboolean isCopy = false;
//   jint * arr = jenv->GetIntArrayElements(result, &isCopy);
//         
//   // Check for exception
//   jthrowable exc = jenv->ExceptionOccurred();
//   if (exc) {
//     std::cout << "Exception occured!!" << std::endl;
//     /* We don't do much with the exception, except that
//        we print a debug message for it, clear it, and 
//        throw a new exception. */
//     jclass newExcCls;
//     jenv->ExceptionDescribe();
//     jenv->ExceptionClear();
//     newExcCls = jenv->FindClass("java/lang/IllegalArgumentException");
//     if (newExcCls == NULL) {
//       /* Unable to find the exception class, give up. */
//       return;
//     }
//     jenv->ThrowNew(newExcCls, "thrown from C code");
//   }
//         
//   for(int i=0; i<task_sz; i++) {
//     // TODO: support multiple tasks, priority
//     scheduler.add_task(gl_types::update_task(arr[i], functions[arr[i+task_sz]]), arr[i+2*task_sz]/(1.0e6));
//   }
//   jenv->ReleaseIntArrayElements(result, arr, JNI_ABORT);
// }  

/**
 * Proxy updater
 * Forwards update calls to the Java layer
 */
class proxy_updater : 
  public graphlab::iupdate_functor<graph_type, proxy_updater> {
};

/**
 * GraphLab core engine
 * For some reason, pthread assertions throw up when this is statically
 * allocated. So this implementation dynamically allocates the core in
 * the initGraphLab function.
 */
static graphlab::core<graph_type, proxy_updater> *core = NULL;

#ifdef __cplusplus
extern "C" {
#endif
  
  /*
   * Class:     graphlab_wrapper_GraphLabJNIWrapper
   * Method:    initGraphLab
   * Signature: ()V
   */
  JNIEXPORT void JNICALL
  Java_graphlab_wrapper_GraphLabJNIWrapper_initGraphLab
  (JNIEnv * env, jobject obj) {
  
    // configure log level
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
 
    // setup the parser
    graphlab::command_line_options clopts("JNI options. TODO.");

    // set the engine options
    if (NULL == core) core = new graphlab::core<graph_type, proxy_updater>();
    core->set_options(clopts);
    
    logger(LOG_DEBUG, "GraphLab core initialized in JNI.\n");
    
  }


//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setScheduler
//   (JNIEnv * env, jobject obj, jstring schedulertype) {
//     const char *str = env->GetStringUTFChars(schedulertype, 0);
//     core.set_scheduler_type(std::string(str));
//   }
// 
//   /*
//    * Class:     graphlab_wrapper_GraphLabJNIWrapper
//    * Method:    setScopeType
//    * Signature: (Ljava/lang/String;)V
//    */
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setScopeType
//   (JNIEnv * env, jobject obj, jstring scopetype) {
//     const char *str = env->GetStringUTFChars(scopetype, 0);
//     core.set_scope_type(std::string(str));
//   }
// 
//   /*
//    * Class:     graphlab_wrapper_GraphLabJNIWrapper
//    * Method:    createGraph
//    * Signature: (I)V
//    */
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_createGraph
//   (JNIEnv * env, jobject obj, jint numvertices) {
//     // Create the dummy vertices
//     jni_graph & graph = core.graph();
//     for(int i=0; i<numvertices; i++) {
//       graph.add_vertex(vertex_data());
//     }
//   }
// 
//   /*
//    * Class:     graphlab_wrapper_GraphLabJNIWrapper
//    * Method:    addEdges
//    * Signature: (I[I)V
//    */
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_addEdges
//   (JNIEnv * env, jobject obj, jint fromVertex, jintArray toVertices) {
//     jni_graph & graph = core.graph();
//     jsize sz = env->GetArrayLength(toVertices);
//     jboolean isCopy = false;
//     jint * arr = env->GetIntArrayElements(toVertices, &isCopy);
//     for(int i=0; i<sz; i++) {
//       graph.add_edge(fromVertex, arr[i]);
//     }
//     env->ReleaseIntArrayElements(toVertices, arr, JNI_ABORT);
//   }
//     
// 
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
//   /*
//    * Class:     graphlab_wrapper_GraphLabJNIWrapper
//    * Method:    runGraphlab
//    * Signature: ()V
//    */
//   JNIEXPORT void JNICALL 
//   Java_graphlab_wrapper_GraphLabJNIWrapper_runGraphlab(JNIEnv * env, jobject obj) {
//     // HACK - where to get number of workers?
//     envs.clear(); envs.resize(512, NULL);
// 
//     if (jvm == NULL) env->GetJavaVM(&jvm);
// 
//     cachedObj = obj;
//     jclass clazz = env->GetObjectClass(obj);
//     wrapperMethodID = env->GetMethodID(clazz, "execUpdate", "(II)[I");
//     std::cout << "Got method ID " << wrapperMethodID << std::endl;
//     std::cout << "Graph has: " << core.graph().num_vertices() << " vertices and " << 
//       core.graph().num_edges() << " edges." << std::endl;
//     core.get_engine_options().print();
//     if (taskbudget>0)
//       core.engine().set_task_budget(taskbudget);
//     if (maxiter>0)
//       core.sched_options().add_option("max_iterations", maxiter);
//  
//     // Set thread destroy callback
//     graphlab::thread::set_thread_destroy_callback(detach_from_JVM);
//    
//     double runtime = core.start(); 
//     if (metrics_type != "none") {
//       core.set_metrics_type(metrics_type);
//       core.fill_metrics();
//       core.report_metrics();
//       // Hack: this prevents core destructor from dumping the metrics.
//       // ... which leads to some weird mutex error.
//       core.set_metrics_type("none");
//     }
//     std::cout << "Finished after " 
// 	      << core.engine().last_update_count() 
// 	      << " updates." << std::endl;
//     std::cout << "Runtime: " << runtime 
// 	      << " seconds." << std::endl;
//   }
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
//   JNIEXPORT void JNICALL Java_graphlab_wrapper_GraphLabJNIWrapper_setNumCPUs
//   (JNIEnv * env, jobject obj, jint ncpus) {
//     core.set_ncpus(ncpus);
//   }
// 
// 
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