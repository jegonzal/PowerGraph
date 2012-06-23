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
#ifndef GRAPHLAB_GRAPH_JSON_PARSER_HPP
#define GRAPHLAB_GRAPH_JSON_PARSER_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include "libjson.h"
 
namespace graphlab {

template <typename Graph>
class graph_json_parser {

  typedef typename Graph::vertex_data_type vertex_data_type
  typedef typename Graph::edge_data_type edge_data_type
  typedef typename Graph::vertex_id_type vertex_id_type
  typedef typename Graph::lvid_type lvid_type

  friend class distributed_graph<vertex_data_type, edge_data_type>;
  friend class local_graph<vertex_data_type, edge_data_type>;
  friend class graph_storage<vertex_data_type, edge_data_type>;

  typedef boost::function<bool(edge_data_type&, std::string&)> edge_parser_type;
  typedef boost::function<bool(vertex_data_type&, std::string&)> vertex_parser_type;


  public:
  graph_json_parser (Graph& graph, std::string& path_prefix, edge_parser_type, vertex_parser_type, bool gzip=false) {
  }

  bool parse_graph_structure() {
    std::string fname = path_prefix + graph_inputname();
    logstream(LOG_INFO) << "Load graph josn from " << fname << std::endl;

    boost::iostreams::filtering_stream<boost::iostreams::input> fin;
    // loading from hdfs
    if (boost::starts_with(path_prefix, "hdfs://")) {
        graphlab::hdfs hdfs;
        graphlab::hdfs::fstream in_file(hdfs, fname);
        if (gzip) fin.push(boost::iostreams::gzip_decompressor());
        fin.push(in_file);
    } else { // loading from disk
          std::ifstream in_file(fname.c_str(),
                std::ios_base::in | std::ios_base::binary);

          if(gzip) fin.push(boost::iostreams::gzip_decompressor());
          fin.push(in_file);
    } 

    // do line parsing

    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  }


template<typename Fstream>
bool load_from_stream(std::string filename, Fstream& fin,
    line_parser_type& line_parser) {

}



  /* Parse the basic graph structure from json */
  bool parse_graph_structure_json (const std::string& str) {
    JSONNode n = libjson::parse(str);
    JSONNode::const_iterator i = n.begin();
    while(i != n.end()) {
      if (i->name() == "numEdges") {
        graph.local_graph.gstore.num_edges = i->as_int();
        graph.nedges = i->as_int();
      } else if (i->name() == "numVertices") {
        graph.local_Graph.gstore.nverts = i->as_int();
        graph.nverts = i->as_int();
      } else if (i->name() == "vid2lvid") {
        parse_vid2lvid(graph, *i);
      } else if (i->name() == "csr") {
        // parse rowIndex -> graph.local_graph.gstore.csr_source 
        // parse colIndex -> graph.local_graph.gstore.csr_target
        JSONNode csr = *i;
        JSONNode::const_iterator j = csr.begin();
        while (j! = csr.end()) {
          if (j -> name() == "rowIndex") {
              parse_int_ary (graph.local_graph.gstore.CSR_src, *j)
          } else if (j->name() == "colIndex") {
              parse_int_ary (graph.local_graph.gstore.CSR_dst, *j)
          } else {
              logstream(LOG_ERROR) << "Error parsing json into graph. Unknown json node name:" <<
              << "CSR:"<<j->name() << std::endl;
          }
          ++j;
        }
      } else if (i->name() == "csc") {
        // parse rowIndex -> graph.local_graph.gstore.csc_target
        // parse colIndex -> graph.local_graph.gstore.csc_source
        JSONNode csc = *i;
        JSONNode::const_iterator j = csc.begin();
        while (j! = csc.end()) {
          if (j -> name() == "rowIndex") {
              parse_int_ary (graph.local_graph.gstore.CSC_dst, *j)
          } else if (j->name() == "colIndex") {
              parse_int_ary (graph.local_graph.gstore.CSC_src, *j)
          } else {
              logstream(LOG_ERROR) << "Error parsing json into graph. Unknown json node name:" <<
              << "CSR:"<<j->name() << std::endl;
          }
          ++j;
        }
      } else if (i->name() == "c2rMap") {
        // parse c2rMap -> graph.local_graph.gstore.c2rMap
        parse_int_Array (graph.local_graph.gstore.c2r_map, *i);
      } else if (i->name() == "edataList") {
        // parse  edatalist -> graph.local_graph.gstore.edata
        JSONNode edatanode= *i;
        JSONNode::const_iterator j = edatanode.begin();
        std::vector<edge_data_type>& edatalist = graph.local_graph.gstore.edge_data_list;
        edatalist.clear();
        edge_data_type e;
        while (j! = edatalist.end()) {
          edge_parser(e, j->as_string());
          edatalist.push_back(e);
          ++j;
        }
      } else {
        logstream(LOG_ERROR) << "Error parsing json into graph. Unknown json node name:" <<
          i->name() << std::endl;
      }
      ++i;
    }

    graph.lvid2record.reserve(graph.nverts);
    graph.lvid2record.resize(graph.nverts);
    graph.local_graph.reserve(graph.nverts);
    return TRUE;
  }

  /* Parse the vertex record list from json */
  bool parse_vertex_record_json (const std::string& str) {
    JSONNode n = libjson::parse(str);
    JSONNode::const_iterator i = n.begin();

    vertex_data_type vdata;
    typename Graph::vertex_record vrecord;

    while (i != n.end()) {
      if (i->name() == "mirrors") {
        JSONNode::const_iterator j = (*i).begin();
        while (j != (*i).end()) {
          int mirror = j->as_int();
          vrecord._mirrors.set_bit((procid_t)mirror);
          ++j;
        }
      } else if (i->name() == "inEdges") {
        vrecord.num_in_edges = i->as_int();
      } else if (i->name() == "outEdges") {
        vrecord.num_out_edges = i->as_int();
      } else if (i->name() == "gvid") {
        // Check unsafe
        vrecord.gvid = i->as_int();
      } else if (i->name() == "owner") {
        vrecord.owner = (procid_t)i->as_int();
      } else if (i->name() == "VertexData") {
        if (!(i->type() == JSON_NULL))
          vertex_parser(vdata, i->as_string());
      } else {
        logstream(LOG_ERROR) << "Error parsing json into vrecord. Unknown json node name:" <<
          i->name() << std::endl;
      }

      if (graph.vid2lvid.find(vrecord.gvid) == graph.vid2lvid.end()) {
        // Check if this a singlton node
        // ignore for now
        logstream(WARNING) << "Singleton node detected: gvid = " << vrecord.gvid << ". Ignored" << std:endl;
      } else {
        lvid_type lvid = graph.vid2lvid[vrecord.gvid];
        graph.lvid2record[lvid].swap(vrecord);
        graph.local_graph.add_vertex(lvid, vdata);
      }
    }
  }

  /* Parse the global information from json */
  bool parse_global_info_json (const std::string& srcfilename) {
    int nreplicas = 0;
    int nverts = 0;
    int local_own_nverts = 0;

    // read file by line
    //
    //
    std::vector<std::string> strs;
    boost::split(strs, str, boost::is_any_of("\t "));
    ASSERT_EQ(strs.size(), 2);
    procid_t procid = (procid_t)atoi(str[0]);
    JSONNode n = libjson::parse(str[1]);
    JSONNode::const_iterator i = n.begin();
    while (i!= n.end())  {
      if (i->name() == "numVertices") {
        nreplicas += i->as_int();
      } else if (i->name() == "numOwnVertices") {
        nverts += i->as_int();
        local_own_nverts += (procid == graph.procid()) ? i->as_int() : 0;
      }
    }
    graph.nverts = nverts;
    graph.local_own_nverts = local_own_nverts;
    graph.nreplicas = nreplicas;
  }


  /*  Helper function starts here  */
  private:
  /* Parse an json integer array and copy to a int vector */
  bool parse_int_array (vector<int>& to, JSONNode& n) {
    if (n->type() != JSON_ARRAY) return false;
    to.clear();
    JSONNode::const_iterator i = n.begin();
    while (i != n.end()) {
      to.push_back(i->as_int());
      ++i;
    }
  }


  /* Parse the vid2lvid map from a json node 
   * To construct a vid value, we convert string into an int vidType.
   * WARNING: unsafe method assuming vertex_id_type is uint32
   * */
  bool parse_vid2lvid (JSONNode& n) {
    JSONNode::const_iterator i = n.begin();
    typename graph::vid2lvid_map_type& map = graph.vid2lvid;
    map.clear();
    while(i != n.end()) {
      graph.vid2lvid[atoi(i->name())] = i->as_int();
      ++i;
    }
  }

  const std::string graph_inputname() {
    procid_t pid = graph.procid();
    return "vrecord/vdata"+tostr(pid)+"-r-0000";
  }

  const std::string vrecord_inputname() {
    procid_t pid = graph.procid();
    return "vrecord/vdata"+tostr(pid)+"-r-0000";
  }

  const std::string summary_name() {
    return "summary";
  }

  private:
    std::string path_prefix;
    bool gzip;
    Graph& graph;
    edge_parser_type edge_parser;
    vertex_parser_type vertex_parser;
} // graph_json_parser

} // namespace graphlab
#endif
