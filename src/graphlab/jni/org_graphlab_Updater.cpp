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

jmethodID proxy_updater::java_method_id = 0;

#ifdef __cplusplus
extern "C" {
#endif

  JNIEXPORT jlong JNICALL
  Java_org_graphlab_Updater_createUpdater
  (JNIEnv *env, jobject updater){
    proxy_updater *proxy = new proxy_updater(env, updater);
    return (long) proxy;
  }

#ifdef __cplusplus
}
#endif