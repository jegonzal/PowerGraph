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
typedef json_message jm;
typedef dispatcher_update dp;

/////////////////////////////// INSTANCE MEMBERS ///////////////////////////////

jm::
  json_message(const std::string method, const std::string state) :
  mdocument() {

  // make object root
  mdocument.SetObject();
  
  // add method if exists
  json::Value methodv(method.c_str(), mdocument.GetAllocator());
  mdocument.AddMember("method", methodv, mdocument.GetAllocator());
  
  // add state if exists
  json::Value statev (state.c_str(), mdocument.GetAllocator());
	mdocument.AddMember("state", statev, mdocument.GetAllocator());
	  
}

void jm::parse(byte *data, std::size_t bytes){
  
  CHECK(NULL != data);
  
  // assume null-terminated string for now
  if (mdocument.Parse<0>(data).HasParseError()){/* TODO: error handling */}
  
}

jm::
  ~json_message(){
}

// important lesson: must use graphlab::
std::ostream& graphlab::operator<< (std::ostream &out, jm &message){
  rapidjson::StringBuffer buffer;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
  message.mdocument.Accept(writer);
  return out << buffer.GetString();
}

const char * jm::updater(){
  if (mdocument.HasMember("updater"))
    return mdocument["updater"].GetString();
  return NULL;
}

const char * jm::vertex(){
  if (mdocument.HasMember("vertex"))
    return mdocument["vertex"].GetString();
  return NULL;
}

void jm::add_context(dp::icontext_type& context, byte flags){

  json::Value contextv;
  contextv.SetObject();
  
  // add params if it does not exists
  if (!mdocument.HasMember("params")){
    json::Value paramsv;
    paramsv.SetObject();
    mdocument.AddMember("params", paramsv, mdocument.GetAllocator());
  }

  // add vertex state to document
  if (flags & VERTEX > 0){
    json::Value vertexv;
    vertexv.SetObject();
    vertexv.AddMember("id",
      context.vertex_id(),
      mdocument.GetAllocator());
    vertexv.AddMember("state",
      context.const_vertex_data().state.c_str(),
      mdocument.GetAllocator());
    contextv.AddMember("vertex", vertexv, mdocument.GetAllocator());
  }
  
//   if (flags & EDGES > 0){
//   }
  
  mdocument["params"].AddMember("context", contextv, mdocument.GetAllocator());

}

///////////////////////////////// CLASS MEMBERS ////////////////////////////////

#include <graphlab/macros_undef.hpp>
