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

#include "org_graphlab_Aggregator.hpp"

using namespace graphlab;

//---------------------------------------------------------------
// proxy_aggregator static members
//---------------------------------------------------------------

jmethodID proxy_aggregator::java_exec     = 0;
jmethodID proxy_aggregator::java_add      = 0;
jmethodID proxy_aggregator::java_finalize = 0;

void proxy_aggregator::init(JNIEnv *env){

  jclass aggregator_class = env->FindClass("org/graphlab/Aggregator");

  // get the method ID for Aggregator#exec, if we don't have it already
  if (0 == java_exec){
    java_exec =
      env->GetMethodID(aggregator_class, "exec", "(JLorg/graphlab/data/Vertex;)V");
  }
  
  // get the method ID for Aggregator#add, if we don't have it already
  if (0 == java_add){
    java_add =
      env->GetMethodID(aggregator_class, "add", "(Lorg/graphlab/Aggregator;)V");
  }
  
  // get the method ID for Aggregator#finalize, if we don't have it already
  if (0 == java_finalize){
    java_finalize = 
      env->GetMethodID(aggregator_class, "finalize", "(J)V");
  }
  
  if (0 == java_clone){
    java_clone = 
      env->GetMethodID(aggregator_class, "clone", "()Lorg/graphlab/Aggregator;");
  }
  
}

//---------------------------------------------------------------
// proxy_aggregator instance members
//---------------------------------------------------------------

proxy_aggregator::
  proxy_aggregator(JNIEnv *env, jobject &obj)
  : java_any(env, obj){}

proxy_aggregator::proxy_aggregator(){}

proxy_aggregator::
  proxy_aggregator(const proxy_aggregator& other)
  : java_any(other) {}

proxy_aggregator &proxy_aggregator::operator=(const proxy_aggregator& other){
    
  if (this == &other) return *this;
  java_any::operator=(other);
  return *this;
  
}

proxy_aggregator::~proxy_aggregator(){}

//---------------------------------------------------------------
// proxy_aggregator instance members - the update function
//---------------------------------------------------------------

void proxy_aggregator::operator()(icontext_type& context){
  
  JNIEnv *env = core::get_jni_env();
  
  // forward call to org.graphlab.Aggregator
  env->CallVoidMethod (obj(), java_exec,
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
// proxy_aggregator instance members - the add function
//---------------------------------------------------------------

void proxy_aggregator::operator+=(const proxy_aggregator& other) {
  
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

//---------------------------------------------------------------
// proxy_aggregator instance members - the finalize function
//---------------------------------------------------------------

void proxy_aggregator::finalize(iglobal_context& context){

  JNIEnv *env = core::get_jni_env();
  
  // forward call to org.graphlab.Aggregator
  env->CallVoidMethod (obj(), java_finalize, &context);
  
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