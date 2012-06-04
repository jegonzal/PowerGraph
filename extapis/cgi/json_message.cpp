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

typedef json_schedule js;
typedef json_message jm;
typedef json_invocation ji;
typedef json_return jr;
typedef dispatcher_update dp;

/////////////////////////////// JSON SCHEDULE //////////////////////////////////

js::json_schedule(js::neighbors_enum targets, const std::string &updater) :
  mupdater(updater),
  mtargets(targets),
  mvertices(0){
}

js::json_schedule
  (js::neighbors_enum targets,
    const std::string &updater,
    const json::Value &vertices) :
  mupdater(updater),
  mtargets(targets),
  mvertices(vertices.Size()){
  for (json::SizeType i = 0; i < vertices.Size(); i++)
    mvertices[i] = vertices[i].GetInt();
}

js::neighbors_enum js::
  targets() const {
  return mtargets;
}

const std::string& js::
  updater() const {
  return mupdater;
}

const std::vector<unsigned>& js::
  vertices() const {
  return mvertices;
}

/////////////////////////////// JSON MESSAGE ///////////////////////////////////

jm::json_message() : mdocument(), mallocator(mdocument.GetAllocator()) {
  // make object root
  mdocument.SetObject();
}

void jm::parse(const byte *data, std::size_t bytes){
  
  CHECK(NULL != data);
  
  // assume null-terminated string for now
  if (mdocument.Parse<0>(data).HasParseError()){/* TODO: error handling */}
  
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
  // add method if exists
  json::Value methodv(method.c_str(), mallocator);
  mdocument.AddMember("method", methodv, mallocator);
  
  // add state if exists
  json::Value statev (state.c_str(), mallocator);
	mdocument.AddMember("state", statev, mallocator);
}

ji::~json_invocation(){}

json::Value&
ji::create_vertex(
  json::Value& vertexv,
  const dp::graph_type::vertex_id_type vertex_id,
  const dp::graph_type::vertex_data_type& vertex_data){
  
  // construct JSON vertex object
  vertexv.SetObject();
  vertexv.AddMember("id", vertex_id, mallocator);
  vertexv.AddMember("state", vertex_data.state.c_str(), mallocator);
  return vertexv;
  
}

json::Value&
ji::create_edge(
  json::Value& edgev,
  const dp::icontext_type& context,
  const dp::graph_type::edge_type& edge){
  
  // construct JSON edge object
  edgev.SetObject();
  edgev.AddMember("state", context.const_edge_data(edge).state.c_str(), mallocator);
  
  // add source and target vertices
  json::Value vertexv;
  edgev.AddMember("source",
    create_vertex(vertexv, edge.source(), context.const_vertex_data(edge.source())),
    mallocator);
  edgev.AddMember("target",
    create_vertex(vertexv, edge.target(), context.const_vertex_data(edge.target())),
    mallocator);
  
  return edgev;
  
}

void ji::add_in_edges(dp::icontext_type& context, json::Value& parent){
  
  // construct JSON in_edges array
  json::Value in_edgesv;
  in_edgesv.SetArray();
  
  // TODO: optimize by removing target vertex
  
  json::Value edgev; // reusable JSON edge object
  foreach(dp::graph_type::edge_type edge, context.in_edges()){
    in_edgesv.PushBack(create_edge(edgev, context, edge), mallocator);
  }
  
  parent.AddMember("in_edges", in_edgesv, mallocator);
  
}

void ji::add_out_edges(dp::icontext_type& context, json::Value& parent){

  json::Value out_edgesv;
  out_edgesv.SetArray();
  
  // TODO: optimize by removing source vertex
  
  json::Value edgev; // reusable JSON edge object
  foreach(dp::graph_type::edge_type edge, context.out_edges()){
    out_edgesv.PushBack(create_edge(edgev, context, edge), mallocator);
  }
  
  parent.AddMember("out_edges", out_edgesv, mallocator);

}

void ji::add_context(dp::icontext_type& context, byte flags){

  json::Value contextv;
  contextv.SetObject();
  
  // add params if it does not exists
  if (!mdocument.HasMember("params")){
    json::Value paramsv;
    paramsv.SetObject();
    mdocument.AddMember("params", paramsv, mallocator);
  }

  // add vertex state to document
  if (flags & VERTEX > 0){
    json::Value vertexv;
    contextv.AddMember("vertex",
      create_vertex(vertexv, context.vertex_id(), context.const_vertex_data()),
      mallocator);
  }
  
  // add edge states to document
  if (flags & EDGES > 0){
    add_in_edges(context, contextv);   // in_edges : [<edge>, ... , <edge>]
    add_out_edges(context, contextv);  // out_edges: [<edge>, ..., <edge>]
  }
  
  mdocument["params"].AddMember("context", contextv, mallocator);

}

void ji::add_edge(dp::icontext_type& context, const dp::graph_type::edge_type& edge){

  // add params if it does not exists
  if (!mdocument.HasMember("params")){
    json::Value paramsv;
    paramsv.SetObject();
    mdocument.AddMember("params", paramsv, mallocator);
  }

  // JSON edge object
  json::Value edgev;
  create_edge(edgev, context, edge);
  
  mdocument["params"].AddMember("edge", edgev, mallocator);

}

void ji::add_other(const std::string& other){

  // add params if it does not exists
  if (!mdocument.HasMember("params")){
    json::Value paramsv;
    paramsv.SetObject();
    mdocument.AddMember("params", paramsv, mallocator);
  }
  
  mdocument["params"].AddMember("other", other.c_str(), mallocator);

}

/////////////////////////////// JSON RETURN ////////////////////////////////////

jr::~json_return(){}

const char * jr::updater() const {
  if (mdocument.HasMember("updater"))
    return mdocument["updater"].GetString();
  return NULL;
}

const char * jr::vertex() const {
  if (mdocument.HasMember("vertex"))
    return mdocument["vertex"].GetString();
  return NULL;
}

const char * jr::scatter_schedule() const {
  if (mdocument.HasMember("schedule") && mdocument["schedule"].IsString())
    return mdocument["schedule"].GetString();
  return NULL;
}

const js jr::schedule() const {
  
  if (!mdocument.HasMember("schedule") || !mdocument["schedule"].IsObject()){
    return json_schedule(js::NONE);
  }
  
  // TODO: throw error if necessary parts are missing
  std::string updater = std::string(mdocument["schedule"]["updater"].GetString());
  if (mdocument["schedule"]["vertices"].IsArray()){
    return json_schedule(js::SOME, updater, mdocument["schedule"]["vertices"]);
  }else {
    const char *vertices = mdocument["schedule"]["vertices"].GetString();
    js::neighbors_enum neighbors;
    if (0 == strcmp(vertices, "ALL"))
      neighbors = js::ALL;
    else if (0 == strcmp(vertices, "IN_NEIGHBORS"))
      neighbors = js::IN_NEIGHBORS;
    else if (0 == strcmp(vertices, "OUT_NEIGHBORS"))
      neighbors = js::OUT_NEIGHBORS;
    else
      neighbors = js::NONE;
    return js(neighbors, updater);
  }
  
}

const consistency_model jr::consistency() const {

  if (!mdocument.HasMember("consistency")){
    // TODO: throw error
  }
  
  std::string consistency = std::string(mdocument["consistency"].GetString());
  if ("VERTEX" == consistency) return VERTEX_CONSISTENCY;
  if ("EDGE" == consistency) return EDGE_CONSISTENCY;
  if ("FULL" == consistency) return FULL_CONSISTENCY;
  if ("NULL" == consistency) return NULL_CONSISTENCY;
  
  // TODO: throw error
  return NULL_CONSISTENCY;

}

const edge_dir_type jr::edge_dir() const {

  if (!mdocument.HasMember("edges")){
    // TODO: throw error
  }
  
  std::string edge_set = std::string(mdocument["edges"].GetString());
  if ("IN_EDGES" == edge_set) return IN_EDGES;
  if ("OUT_EDGES" == edge_set) return OUT_EDGES;
  if ("ALL_EDGES" == edge_set) return ALL_EDGES;
  if ("NO_EDGES" == edge_set) return NO_EDGES;
  
  // TODO: throw error
  return graphlab::NO_EDGES;

}

#include <graphlab/macros_undef.hpp>
