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
 * @file dispatcher.hpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
 
#ifndef GRAPHLAB_DISPATCHER_HPP
#define GRAPHLAB_DISPATCHER_HPP

#include <graphlab.hpp>

////////////////////////////// CGI VERTEX PROGRAM //////////////////////////////

namespace graphlab {
  
  class json_return; // forward declaration
  
  /**
   * The gather type for CGI programs is just a string. To merge two gathers,
   * a merge is invoked on the child process.
   */
  class cgi_gather {
  private:
    std::string mstate;
  public:
    cgi_gather(const std::string& state="");
    void operator+=(const cgi_gather& other);
    const char *c_str() const;
    void save(oarchive& oarc) const;
    void load(iarchive& iarc);
  }; // end of cgi_gather_type
  
  /**
   * The dispatcher vertex program extends ivertex_program and specifies the
   * vertex type, edge type, and gather type. It receives calls from the engine
   * and forwards them to child processes using JSON via UNIX pipes.
   * @internal
   */  
  class dispatcher :
    public ivertex_program<distributed_graph<std::string, std::string>, cgi_gather> {
  public:
    /** vertex state is just a string - user decides format */
    typedef std::string cgi_vertex;
    /** edge state is just a string - user decides format */
    typedef std::string cgi_edge;
    /** distributed graph with cgi_vertex and cgi_edge */
    typedef distributed_graph<cgi_vertex, cgi_edge> graph_type;
    typedef cgi_gather gather_type;
  private:
    /** vertex program state - */
    std::string mstate;
  public:
    dispatcher(const std::string& state="");
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const;
    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const;
    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const;
    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total);
    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const;
    void save(oarchive& oarc) const;
    void load(iarchive& iarc);
    static void init_edge(graph_type::edge_type edge);
    static void init_vertex(graph_type::vertex_type vertex);
    static bool read_graph(graph_type& graph,
                          const std::string& filename,
                          const std::string& textline);
  }; // end of dispatcher
  
};

////////////////////////////////////////////////////////////////////////////////

#endif /* #ifndef GRAPHLAB_DISPATCHER_HPP */

