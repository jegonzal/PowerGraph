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
jmethodID proxy_updater::java_clone    = 0;

// initialize JNI method IDs
void proxy_updater::init(JNIEnv *env){

  jclass updater_class = env->FindClass("org/graphlab/Updater");

  // get the method ID for Updater#execUpdate, if we don't have it already
  if (0 == java_exec_update){
    java_exec_update =
      env->GetMethodID(updater_class, "execUpdate", "(JLorg/graphlab/data/Vertex;)V");
  }
  
  // get the method ID for Updater#add, if we don't have it already
  if (0 == java_add){
    java_add =
      env->GetMethodID(updater_class, "add", "(Lorg/graphlab/Updater;)V");
  }
  
  // get the method ID for Updater#priority, if we don't have it already
  if (0 == java_priority){
    java_priority = 
      env->GetMethodID(updater_class, "priority", "()D");
  }
  
  // get the method ID for Updater#clone, if we don't have it already
  if (0 == java_clone){
    java_clone = 
      env->GetMethodID(updater_class, "clone", "()Lorg/graphlab/Updater;");
  }
  
}

//---------------------------------------------------------------
// proxy_updater instance members
//---------------------------------------------------------------

proxy_updater::
  proxy_updater(JNIEnv *env, jobject &obj)
  : java_any(env, obj){}

proxy_updater::proxy_updater(){}

proxy_updater::
  proxy_updater(const proxy_updater& other){

  // other doesn't have an existing ref
  if (NULL == other.obj()){
    set_obj(NULL);
    return;
  }
  
  // clone the java object
  JNIEnv *env = core::get_jni_env();
  set_obj(env->CallObjectMethod(other.obj(), java_clone));  
  
}

proxy_updater &proxy_updater::operator=(const proxy_updater& other){
    
  if (this == &other) return *this;
  java_any::operator=(other);
  return *this;
  
}

proxy_updater::~proxy_updater(){}

//---------------------------------------------------------------
// proxy_updater instance members - the update function
//---------------------------------------------------------------

void proxy_updater::operator()(icontext_type& context){
  
  JNIEnv *env = core::get_jni_env();
  
  // forward call to org.graphlab.Updater
  env->CallVoidMethod (obj(), java_exec_update,
                       &context,
                       context.vertex_data().app_vertex);
  
  // check for exception
  jthrowable exc = env->ExceptionOccurred();
  if (exc) {
  
    // TODO: better error handling
      
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

void proxy_updater::operator+=(const proxy_updater& other){
  
  JNIEnv *env = core::get_jni_env();
  
  // forward call to org.graphlab.Updater
  env->CallVoidMethod (obj(), java_add, other.obj());
  
  // check for exception
  jthrowable exc = env->ExceptionOccurred();
  if (exc) {
  
    // TODO: better error handling
      
    jclass new_exc;
    env->ExceptionDescribe();
    env->ExceptionClear();
    new_exc = env->FindClass("java/lang/IllegalArgumentException");
    if (new_exc == NULL) return;
    env->ThrowNew(new_exc, "thrown from C code");
    
  }

}