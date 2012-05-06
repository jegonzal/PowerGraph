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
 * @file message.cpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
 
#include "json_message.hpp"
#include <graphlab/macros_def.hpp>

using namespace graphlab;
namespace json = rapidjson;

/////////////////////////////// INSTANCE MEMBERS ///////////////////////////////

json_message::
  json_message(const std::string method, const std::string state) : mdocument() {

  // make object root
  mdocument.SetObject();
  
  // add method if exists
  if (method.length() > 0){
  	json::Value methodv(method.c_str(), mdocument.GetAllocator());
  	mdocument.AddMember("method", methodv, mdocument.GetAllocator());
  }
  
  // add state if exists
  if (state.length() > 0){
  	json::Value statev (state.c_str(), mdocument.GetAllocator());
	  mdocument.AddMember("state", statev, mdocument.GetAllocator());
	}
  
}

json_message::
  ~json_message(){}

///////////////////////////////// CLASS MEMBERS ////////////////////////////////

#include <graphlab/macros_undef.hpp>
