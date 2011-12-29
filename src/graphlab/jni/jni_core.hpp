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
 
#ifndef JNI_CORE_HPP
#define JNI_CORE_HPP

#include <jni.h>
#include <graphlab.hpp>

namespace graphlab {

  /**
   * Wrapper for graphlab::core. Contains the core and a reference to
   * the Java core object (so that it doesn't get garbage collected.)
   */
  template <typename Graph, typename UpdateFunctor>
  class jni_core {

  public:
    typedef Graph graph_type;
    typedef UpdateFunctor update_functor_type;
    typedef core<graph_type, update_functor_type> core_type;
    
  private:
    
    core_type *mcore;                   /**< graphlab::core */
    jobject mobj;                       /**< java object reference */
    
    /** Java virtual machine reference */
    static JavaVM *mjvm;
    
    /** Map of thread ids to their respective JNI envs */
    static std::vector<JNIEnv *> menvs;
    
  public:

    /**
     * Creates a new graphlab core and a global reference to the
     * Java object.
     */
    jni_core (JNIEnv *env, jobject &obj){
      // allocate a new core
      this->mcore = new core_type();
      // create a new global reference to the java object
      this->mobj = env->NewGlobalRef(obj);
    }
    
    core_type &core(){
      return *mcore;
    }
    
    jobject &obj(){
      return mobj;
    }
    
    ~jni_core(){
    
      delete mcore;
      
      JNIEnv *env;
      mjvm->AttachCurrentThread((void **)&env, NULL);
      env->DeleteGlobalRef(mobj);
      
    }
    
    /**
     * Saves a reference to the Java Virtual Machine
     */
    static void set_jvm (JavaVM *jvm){
      mjvm = jvm;
    }
    
    /**
     * Retrieves a reference to the Java Virtual Machine
     */
    static JavaVM *get_jvm (){
      return mjvm;
    }
    
    /**
     * Detaches the current thread from the JVM
     */
    static void detach_from_jvm() {
    
      logstream(LOG_INFO)
          << "detach_from_JVM invoked."
          << std::endl;
      
      int thread_id = thread::thread_id();
      
      if (NULL != menvs[thread_id]) {
        int res = mjvm->DetachCurrentThread();
        logstream(LOG_INFO)
          << "Detached from JVM: " << res
          << std::endl;
        assert(res>=0);
      }
      
    }
    
    /**
     * Retrieves the JNI environment for the current thread.
     */
    static JNIEnv *get_JNIEnv (){
    
      JNIEnv *jenv = NULL;
      int thread_id = thread::thread_id();
      
      logstream(LOG_DEBUG)
        << "Thread ID: " << thread_id
        << std::endl;
    
      // if current thread is not already on the JVM, attach it    
      if (NULL == menvs[thread_id]) {
        logstream(LOG_INFO)
          << "Attaching thread to JVM ... ";
        int res = mjvm->AttachCurrentThread((void **)&jenv, NULL);
        logstream(LOG_INFO)
          << "Done: " << res
          << std::endl;
        assert(res >= 0);
        menvs[thread_id] = jenv;
      }
      
      // now we have the environment associated with the current thread
      return menvs[thread_id];
      
    }
    
  };
  
  // class initializers
  template<typename G, typename U>
  JavaVM* jni_core<G, U>::mjvm = NULL;
  
  template<typename G, typename U>
  std::vector<JNIEnv *> jni_core<G, U>::menvs(thread::cpu_count());

}

#endif