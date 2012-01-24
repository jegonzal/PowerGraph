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

#include "org_graphlab_Updater.hpp"

using namespace graphlab;

//---------------------------------------------------------------
// proxy_updater static members
//---------------------------------------------------------------

jmethodID proxy_updater::java_exec_update = 0;
jmethodID proxy_updater::java_add = 0;
jmethodID proxy_updater::java_priority= 0;

void proxy_updater::init(JNIEnv *env){

  jclass updater_class = env->FindClass("org/graphlab/Updater");

  // get the method ID for Updater#execUpdate, if we don't have it already
  if (0 == proxy_updater::java_exec_update){
    proxy_updater::java_exec_update =
      env->GetMethodID(updater_class, "execUpdate", "(JI)V");
  }
  
  // get the method ID for Updater#add, if we don't have it already
  if (0 == proxy_updater::java_add){
    proxy_updater::java_add =
      env->GetMethodID(updater_class, "add", "(Lorg/graphlab/Updater;)V");
  }
  
  // get the method ID for Updater#priority, if we don't have it already
  if (0 == proxy_updater::java_priority){
    proxy_updater::java_priority = 
      env->GetMethodID(updater_class, "priority", "()D");
  }
  
}

//---------------------------------------------------------------
// proxy_updater instance members
//---------------------------------------------------------------

proxy_updater::proxy_updater(JNIEnv *env, jobject &java_updater){
  mjava_updater = env->NewGlobalRef(java_updater);
}

proxy_updater::proxy_updater() : mjava_updater(NULL){}

proxy_updater::proxy_updater(const proxy_updater& other){
    
  if (NULL == other.mjava_updater){
    this->mjava_updater = NULL;
    return;
  }
  
  // create another reference
  JNIEnv *env = core::get_jni_env();
  this->mjava_updater = env->NewGlobalRef(other.mjava_updater);
  
}

proxy_updater &proxy_updater::operator=(const proxy_updater& other){
    
  if (this == &other) return *this;
  
  JNIEnv *env = core::get_jni_env();
  jobject java_updater = NULL;
  
  // if other has a java object, create a new ref
  if (NULL != other.mjava_updater)
    java_updater = env->NewGlobalRef(other.mjava_updater);
  
  // if this has a java object, delete ref
  if (NULL != this->mjava_updater)
    env->DeleteGlobalRef(this->mjava_updater);
    
  // assign!
  this->mjava_updater = java_updater;
  
  return *this;
  
}

proxy_updater::~proxy_updater(){
  if (NULL == mjava_updater) return;
  // delete reference to allow garbage collection
  JNIEnv *env = core::get_jni_env();
  env->DeleteGlobalRef(mjava_updater);
  mjava_updater = NULL;
}

//---------------------------------------------------------------
// proxy_updater instance members - the update function
//---------------------------------------------------------------

void proxy_updater::operator()(icontext_type& context){
  
  JNIEnv *env = core::get_jni_env();
  
  // forward call to org.graphlab.Updater
  jint app_vertex_id = context.vertex_data().app_id;
  env->CallVoidMethod (mjava_updater, java_exec_update,
                       &context,
                       app_vertex_id);
  
  // check for exception
  jthrowable exc = env->ExceptionOccurred();
  if (exc) {
  
    // TODO: better error handling
  
    logstream(LOG_ERROR)
      << "Exception occured!!"
      << std::endl;
      
    jclass new_exc;
    env->ExceptionDescribe();
    env->ExceptionClear();
    new_exc = env->FindClass("java/lang/IllegalArgumentException");
    if (new_exc == NULL) return;
    env->ThrowNew(new_exc, "thrown from C code");
    
  }

}

//---------------------------------------------------------------
// proxy_updater instance members - the add function
//---------------------------------------------------------------

void proxy_updater::operator+=(const proxy_updater& other) const {
  
  JNIEnv *env = core::get_jni_env();
  
  // forward call to org.graphlab.Updater
  env->CallVoidMethod (mjava_updater, java_add, other.mjava_updater);
  
  // check for exception
  jthrowable exc = env->ExceptionOccurred();
  if (exc) {
  
    // TODO: better error handling
  
    logstream(LOG_ERROR)
      << "Exception occured!!"
      << std::endl;
      
    jclass new_exc;
    env->ExceptionDescribe();
    env->ExceptionClear();
    new_exc = env->FindClass("java/lang/IllegalArgumentException");
    if (new_exc == NULL) return;
    env->ThrowNew(new_exc, "thrown from C code");
    
  }

}