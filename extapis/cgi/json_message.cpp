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
typedef json_invocation ji;
typedef json_return jr;

/////////////////////////////// JSON MESSAGE ///////////////////////////////////

jm::json_message() : mdocument(), mallocator(mdocument.GetAllocator()) {
  // make object root
  mdocument.SetObject();
}

void jm::parse(const byte *data, std::size_t bytes){
  
  if (NULL == data)
    throw "json must not be empty.";
  
  // TODO: assume null-terminated string for now (may switch to binary)
  if (mdocument.Parse<0>(data).HasParseError())
    throw "json syntax error.";
  
}

// important lesson: must use graphlab::
std::ostream& graphlab::operator<< (std::ostream &out, jm &message){
  rapidjson::StringBuffer buffer;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
  message.mdocument.Accept(writer);
  return out << buffer.GetString();
}

/////////////////////////////// JSON INVOCATION ////////////////////////////////

ji::json_invocation(const std::string& method, const std::string& state){
  // add method
  json::Value methodv(method.c_str(), mallocator);
  mdocument.AddMember("method", methodv, mallocator);
  // add state
  json::Value statev (state.c_str(), mallocator);
	mdocument.AddMember("state", statev, mallocator);
}

ji::~json_invocation(){}

json::Value&
ji::create_vertex(json::Value& vertexv, const dispatcher::graph_type::vertex_type& vertex){
  // construct JSON vertex object
  vertexv.SetObject();
  vertexv.AddMember("id", vertex.id(), mallocator);
  vertexv.AddMember("state", vertex.data().c_str(), mallocator);
  return vertexv;
}

json::Value&
ji::create_edge(json::Value& edgev, dispatcher::graph_type::edge_type& edge){
  
  // construct JSON edge object
  edgev.SetObject();
  edgev.AddMember("state", edge.data().c_str(), mallocator);
  
  // add source and target vertices
  json::Value vertexv;
  edgev.AddMember("source", create_vertex(vertexv, edge.source()), mallocator);
  edgev.AddMember("target", create_vertex(vertexv, edge.target()), mallocator);
  
  return edgev;
  
}

void ji::ensure_params_exist(){
  // add params if it does not exists
  if (!mdocument.HasMember("params")){
    json::Value paramsv;
    paramsv.SetObject();
    mdocument.AddMember("params", paramsv, mallocator);
  }
}

void ji::add_vertex(const dispatcher::vertex_type& vertex){
  ensure_params_exist();
  json::Value vertexv; create_vertex(vertexv, vertex);
  mdocument["params"].AddMember("vertex", vertexv, mallocator);
}
    
void ji::add_gather(const dispatcher::gather_type& gather_total){
  ensure_params_exist();
  mdocument["params"].AddMember("gather", serialize_to_string(gather_total).c_str(), mallocator);
}

void ji::add_edge(dispatcher::graph_type::edge_type& edge){
  ensure_params_exist();
  json::Value edgev; create_edge(edgev, edge);
  mdocument["params"].AddMember("edge", edgev, mallocator);
}

void ji::add_other(const std::string& other){
  ensure_params_exist();
  mdocument["params"].AddMember("other", other.c_str(), mallocator);
}

/////////////////////////////// JSON RETURN ////////////////////////////////////

jr::~json_return(){}

const char * jr::program() const {
  if (mdocument.HasMember("program"))
    return mdocument["program"].GetString();
  return NULL;
}

const char * jr::vertex() const {
  if (mdocument.HasMember("vertex"))
    return mdocument["vertex"].GetString();
  return NULL;
}

const edge_dir_type jr::edge_dir() const {

  if (!mdocument.HasMember("edges"))
    throw "json missing 'edges' field.";
  
  std::string edge_set = std::string(mdocument["edges"].GetString());
  if ("IN_EDGES" == edge_set) return IN_EDGES;
  if ("OUT_EDGES" == edge_set) return OUT_EDGES;
  if ("ALL_EDGES" == edge_set) return ALL_EDGES;
  if ("NO_EDGES" == edge_set) return NO_EDGES;
  
  logstream(LOG_ERROR) << "Unrecognized value: " << edge_set << ". Using NO_EDGES." << std::endl;
  return graphlab::NO_EDGES;

}

const char *jr::result() const {
  if (mdocument.HasMember("result"))
    return mdocument["result"].GetString();
  return NULL;
}

const char *jr::signal() const {
  if (mdocument.HasMember("signal"))
    return mdocument["signal"].GetString();
  return NULL;
}

#include <graphlab/macros_undef.hpp>
