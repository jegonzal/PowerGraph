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
 * The original graph.hpp contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Inteface Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 * Implementation by Danny Bickson, CMU
 *
 */


#ifndef GRAPHLAB_MULTIGRAPH_HPP
#define GRAPHLAB_MULTIGRAPH_HPP

#include <omp.h>
#include <cmath>
#include <stdio.h>
#include <string>
#include <list>
#include <vector>

#include <fstream>

#include <boost/bind.hpp>
#include <boost/unordered_set.hpp>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/util/random.hpp>
#include <graphlab/graph/graph_storage.hpp>
#include <graphlab/graph/graph3.hpp>
#include "../toolkits/shared/io.hpp"
#include <graphlab/macros_def.hpp>


namespace graphlab { 

  /**
 * CSR/CSC implementation of graph.
 * Assumptions: number of nodes and number of edges are below MAX_UINT32
 */

  template<typename VertexData, typename EdgeData>
  class multigraph {

    std::vector<graph3<VertexData,EdgeData> * > graphs;
    std::vector<VertexData> node_vdata_array;
    uint num_nodes;
    EdgeData _edge;
    char _color; //placeholder 
    std::vector<std::string> in_files;

  public:
    /// The type of a vertex is a simple size_t
    //typedef graphlab::vertex_id_type vertex_id_type;
    typedef uint vertex_id_type;    

    /// The type of an edge id 
    //typedef size_t edge_id_type;
    
    /// Type for vertex colors 
    typedef char vertex_color_type;

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;
 
    /** Represents an edge with source() and target()*/
    /** The type of the edge list */
    //typedef typename gstore_type::edge_list edge_list;
    typedef edge_type_impl edge_type;
    typedef edge_type_impl edge_id_type;
 
    /** Interface for iupdate functor.*/
    //typedef typename gstore_type::edge_list edge_list_type;
    typedef edge_list edge_list_type;

    int num_graphs(){ return in_files.size(); };

    graph3<VertexData,EdgeData>* graph(int i){ return graphs[i]; }

    multigraph(){
      num_nodes = 0;
      undirected = false;
    }


    const std::vector<VertexData> *get_node_vdata(){ return &node_vdata_array; }

    edge_list_type in_edges(const vertex_id_type v) const{
      assert(false);
      return in_edges(v,0);
    }
    edge_list_type out_edges(const vertex_id_type v) const{
      assert(false);
      return out_edges(v,0);
    }


    /**
     * Create a graph with nverts vertices.
     */
    //graph3(size_t nverts) { }

    //multigraph(const graph<VertexData, EdgeData>& g) { (*this) = g; }

    // METHODS =================================================================>

    /**
     * Finalize a graph by sorting its edges to maximize the
     * efficiency of graphlab.  
     * This function takes O(|V|log(degree)) time and will 
     * fail if there are any duplicate edges.
     * This is also automatically invoked by the engine at
     * start.
     */
    void finalize() {   
      finalized = true;
    } // End of finalize
            
    /** \brief Get the number of vertices */
    size_t num_vertices() const {
       return num_nodes;
    } // end of num vertices

    /** \brief Get the number of vertices local to this machine */
    size_t local_vertices() const {
      return num_nodes;
    } 

    /** \brief Get the number of edges */
    size_t num_edges() const {
      size_t total = 0;
      for (int i=0; i< (int)graphs.size(); i++)
         total += graphs[i]->num_edges();
      return total;
    } 

    /** \brief Finds an edge.
        The value of the first element of the pair will be true if an 
        edge from src to target is found and false otherwise. If the 
        edge is found, the edge ID is returned in the second element of the pair. */
    edge_type find(const vertex_id_type source,
                   const vertex_id_type target) const {
       assert(false); //not implemented yet
       return edge_type(-1,-1);
   
    } // end of find

    edge_type reverse_edge(const edge_type& edge) const {
      //return gstore.find(edge.target(), edge.source());
       assert(false); //not implemented yet
      return edge_type(-1,-1);
    }


    /** 
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     */
    vertex_id_type add_vertex(const VertexData& vdata = VertexData() ) {
       assert(false); //not implemented yet
      if (finalized)
      {
        logstream(LOG_FATAL)
          << "Attempting add vertex"
          << "to a finalized graph." << std::endl;
        ASSERT_MSG(false, "Add vertex to a finalized graph.");
      }
    } // End of add vertex;


    /** 
     * \brief Add additional vertices up to provided num_vertices.  This will
     * fail if resizing down.
     */
    void resize(size_t num_vertices ) {
      ASSERT_GE(num_vertices, num_nodes());
       assert(false); //not implemented yet
      //TODO
    } // End of resize
    
    
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    edge_id_type add_edge(vertex_id_type source, vertex_id_type target, 
                          const EdgeData& edata = EdgeData()) {
       assert(false); //not implemented yet
      return edge_id_type(source, target); //TODO
    } // End of add edge
        
    
    /** \brief Returns a reference to the data stored on the vertex v. */
    VertexData& vertex_data(vertex_id_type v) {
      ASSERT_LT(v, num_nodes);
      return node_vdata_array[v];
    } // end of data(v)
    
    /** \brief Returns a constant reference to the data stored on the vertex v */
    const VertexData& vertex_data(vertex_id_type v) const {
      ASSERT_LT(v, num_nodes);
      return node_vdata_array[v];
    } // end of data(v)

    /** \brief Returns a reference to the data stored on the edge source->target. */
    EdgeData& edge_data(vertex_id_type source, vertex_id_type target) {
     ASSERT_TRUE(finalized);
       assert(false); //not implemented yet
     return _edge;
    } // end of edge_data(u,v)
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_type source, vertex_id_type target) const {
     ASSERT_TRUE(finalized);
       assert(false); //not implemented yet
     return _edge;
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(edge_type edge) { 
      ASSERT_TRUE(finalized);
       assert(false); //not implemented yet
      return _edge;
    }
    const EdgeData& edge_data(edge_type edge) const {
       assert(false); //not implemented yet
      //return 
    }

    size_t num_in_edges(const vertex_id_type v, int id = 0) const {
      return graphs[id]->num_in_edges(v);
    }

    size_t num_out_edges(const vertex_id_type v, int id = 0) const {
      return graphs[id]->num_out_edges(v);
    }


    edge_list_type in_edges(vertex_id_type v, int id) {
      return graphs[id]->in_edges(v);
    }

    edge_list_type out_edges(vertex_id_type v, int id) {
      return graphs[id]->out_edges(v);
    }

    const edge_list_type in_edges(vertex_id_type v, int id) const {
      return graphs[id]->in_edges(v);
    }

    const edge_list_type out_edges(vertex_id_type v, int id) const {
       return graphs[id]->out_edges(v);
     }



   
    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    const vertex_color_type& color(vertex_id_type vertex) const {
       assert(false); //not implemented yet
      return _color;
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    vertex_color_type& color(vertex_id_type vertex) {
       assert(false); //not implemented yet
      return _color;
    }

    vertex_color_type get_color(vertex_id_type vid) const{
       assert(false); //not implemented yet
      return char(0);
    }
    
    void set_color(vertex_id_type vid, vertex_color_type col) {
       assert(false); //not implemented yet
    }
    
    /** \brief This function constructs a heuristic coloring for the 
        graph and returns the number of colors */
    size_t compute_coloring() {
       assert(false); //not implemented yet
      return -1;
    } // end of compute coloring


    /**
     * \brief Check that the colors satisfy a valid coloring of the graph.
     * return true is coloring is valid;
     */
    bool valid_coloring() const {
       assert(false); //not implemented yet
      return true;
    }
    
    /** \brief count the number of times the graph was cleared and rebuilt */
    size_t get_changeid() const {
       assert(false); //not implemented yet
      return 0;
    }

    size_t estimate_sizeof() const {
      /*size_t vlist_size = sizeof(vertices) + sizeof(VertexData) * vertices.capacity();
      size_t vcolor_size = sizeof(vcolors) + sizeof(vertex_color_type) * vcolors.capacity();
      size_t elist_size = edges_tmp.estimate_sizeof(); 
      size_t store_size = gstore.estimate_sizeof();*/

//      printf("graph3: tmplist size: %u, gstoreage size: %u \n", elist_size, store_size);
       return num_nodes*sizeof(uint)+num_edges*sizeof(int);
    }

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
       assert(false); //not implemented yet
    } // end of load

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
       assert(false); //not implemented yet
    } // end of save

    std::string reference_graph_name(int ref){
       assert(in_files.size() > 0);
       return in_files[ref];
    }  
 
    void doload(int i){
       graphlab::timer mt; mt.start();
       graph3<VertexData,EdgeData> *graph = new graph3<VertexData,EdgeData>();
       graph->load_directed(in_files[i], true, false);
       logstream(LOG_INFO)<<"Time taken to load: " << mt.current_time() << std::endl;
       num_nodes = graph->num_vertices();
       if (node_vdata_array.size() == 0)
         node_vdata_array.resize(num_nodes);
       graph->set_node_vdata_array(&node_vdata_array);
       graphs.push_back(graph);
        }

    void unload_all(){
       for (int i=0; i< num_graphs(); i++){
           graphs[i]->clear();
           delete graphs[i];
       }
       graphs.clear();
    } 

    void unload(int i){
      graphs[i]->clear();
      delete graphs[i];
      graphs.erase(graphs.begin()+i);
    }

    /** \brief Load the graph from a file */
    void load(const std::string & listdir, const std::string& dirname, const std::string & prefix, bool delayed) {
     in_files = list_all_files_in_dir(listdir, prefix);
     for (int i=0; i< (int)in_files.size(); i++)
       in_files[i] = dirname + in_files[i];

     if (!delayed){
     for (int i=0; i< (int)in_files.size(); i++){
        graph3<VertexData, EdgeData> graph;
       logstream(LOG_INFO)<<"loading graph " << i <<"/" << in_files.size() << " " << dirname << in_files[i] << std::endl;
        doload(i);
     }
     logstream(LOG_INFO)<<"Total edges read: " << num_edges() << std::endl;
     num_nodes = graph(0)->num_vertices();
     }
     else logstream(LOG_INFO)<<"preparing to load " << in_files.size() << " input graphs" << std::endl;
    } // end of load

    /**
     * \brief save the adjacency structure to a text file.
     *
     * Save the adjacency structure as a text file in:
     *    src_Id, dest_Id \n
     *    src_Id, dest_Id \n
     * format.
     */
    void save_adjacency(const std::string& filename) const {
       assert(false); //not implemented yet
    }



    /**
     * builds a topological_sort of the graph returning it in topsort. 
     * 
     * \param[out] topsort Resultant topological sort of the graph vertices.
     *
     * function will return false if graph is not acyclic.
     */
    bool topological_sort(std::vector<vertex_id_type>& topsort) const {
       assert(false); //not implemented yet
      return true;
    } // end of topological sort

    void set_undirected(){ undirected = true; }

    
  private:    
    /** Internal edge class  */   
    bool finalized;
    bool undirected;
    
  }; // End of class graph3

} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>
#endif

