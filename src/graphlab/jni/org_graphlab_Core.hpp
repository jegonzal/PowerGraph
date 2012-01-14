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
 * @file org_graphlab_Core.hpp
 * \c javah will generate \c org_graphlab_Core.h from the native methods
 * defined in \c org.graphlab.Context (and so will overwrite the file every time).
 * Define any additional classes/structs/typedefs in this hpp file.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */

#ifndef ORG_GRAPHLAB_CORE_HPP
#define ORG_GRAPHLAB_CORE_HPP

#include <graphlab.hpp>
#include "org_graphlab_Core.h"

namespace graphlab {

  /**
   * Wrapper for graphlab::core.
   * Contains the core, a reference to the Java core object (so that it
   * doesn't get garbage collected), and other utility functions for dealing
   * with the JVM.
   */
  template <typename Graph, typename UpdateFunctor>
  class jni_core {

  public:
  
    typedef core<Graph, UpdateFunctor> core_type;
    
  private:
    
    /** graphlab::core object - the soul that this body wraps around */
    core_type *mcore;
    /** associated org.graphlab.Core object */
    jobject mobj;
    
    /** Java virtual machine reference - set only once for each process */
    static JavaVM *mjvm;
    
    /** Map of thread ids to their respective JNI envs - BAD CODE */
    static std::vector<JNIEnv *> menvs;
    
  public:

    /**
     * Creates a new graphlab core and a new reference to the associated
     * Java org.graphlab.Core object (so that it doesn't get garbage collected.)
     * @param[in] env   JNI environment, which will be used to create the
     *                  reference to the Java object.
     * @param[in] obj   associated org.graphlab.Core object.
     */
    jni_core (JNIEnv *env, jobject &obj){
      this->mcore = new core_type();
      this->mobj = env->NewGlobalRef(obj);
    }
    
    /**
     * Gets the real graphlab core that this method wraps around
     * @return graphlab::core
     */
    core_type &operator()(){
      return *mcore;
    }
    
    /**
     * Gets the associated org.graphlab.Core object
     * @return org.graphlab.Core
     */
    jobject &obj(){
      return mobj;
    }
    
    /**
     * Deallocates the graphlab core and deletes the reference to the
     * Java object (so that it may be garbage collected.)
     */
    ~jni_core(){
    
      delete mcore;
      
      // -- BEGIN BAD CODE --
      JNIEnv *env;
      /* TODO: fix - this has to be done because I can't add this env back to
       * the vector.
       */
      mjvm->AttachCurrentThread((void **)&env, NULL);
      // -- END BAD CODE --
      
      // allow associated org.graphlab.Core to be gc'ed
      env->DeleteGlobalRef(mobj);
      
    }
    
    /**
     * Saves a reference to the Java Virtual Machine.
     * @param[in] jvm   pointer to the Java Virtual Machine
     */
    static void set_jvm (JavaVM *jvm){
      mjvm = jvm;
    }
    
    /**
     * Gets a reference to the Java Virtual Machine.
     * @return pointer to the Java Virtual Machine
     */
    static JavaVM *get_jvm (){
      return mjvm;
    }
    
    /**
     * Detaches the current thread from the JVM. If the current
     * thread is associated with an environment that we kept in
     * the vector, remove it.
     * Note: this is unsafe - update after Yucheng provides thread-local
     * storage.
     */
    static void detach_from_jvm() {
      
      int thread_id = thread::thread_id();
      
      if (NULL != menvs[thread_id]) {
        int res = mjvm->DetachCurrentThread();
        logstream(LOG_INFO)
          << "Detached from JVM: " << res
          << std::endl;
        assert(res >= 0);
        menvs[thread_id] = NULL;
      }
      
    }
    
    /** Convenience method for throwing Java exceptions */
    static void throw_exception(JNIEnv* env,
                                const char *exception,
                                const char *message){
      jclass exc = env->FindClass(exception);
      if (NULL == exc) return;
      env->ThrowNew(exc, message);
    }
    
    /**
     * Retrieves the JNI environment for the current thread. This caches
     * the reference to the JNIEnv in a vector.
     * Note: this is unsafe - update after Yucheng provides thread-local
     * storage.
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

}

#endif