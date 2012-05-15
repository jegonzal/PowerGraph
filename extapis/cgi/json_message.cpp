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

///////////////////////////////// C CALLBACKS //////////////////////////////////


/////////////////////////////// INSTANCE MEMBERS ///////////////////////////////

json_message::
  json_message(const std::string method, const std::string state) :
  mdocument() {

  // make object root
  mdocument.SetObject();
  
  // add method if exists
  if (0 < method.length()){
  	json::Value methodv(method.c_str(), mdocument.GetAllocator());
  	mdocument.AddMember("method", methodv, mdocument.GetAllocator());
  }
  
  // add state if exists
  if (0 < state.length()){
  	json::Value statev (state.c_str(), mdocument.GetAllocator());
	  mdocument.AddMember("state", statev, mdocument.GetAllocator());
	}
  
}

bool json_message::feed(byte *data, std::size_t nread){
  return false;
}

json_message::
  ~json_message(){
}

// important lesson: must use graphlab::
std::ostream& graphlab::operator<< (std::ostream &out, json_message &message){
  rapidjson::StringBuffer buffer;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
  message.mdocument.Accept(writer);
  return out << buffer.GetString();
}

///////////////////////////////// CLASS MEMBERS ////////////////////////////////

#include <graphlab/macros_undef.hpp>
