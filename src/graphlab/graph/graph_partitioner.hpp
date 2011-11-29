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
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */


/**
 * \file partition_method.hpp 
 *
 * This file contains a basic struct enumerating supported
 * partitioning methods.
 *
 */

#ifndef GRAPHLAB_GRAPH_PARTITIONER_HPP
#define GRAPHLAB_GRAPH_PARTITIONER_HPP




#include <omp.h>
#include <cmath>

#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>


#include <queue>
#include <algorithm>
#include <functional>
#include <fstream>


#include <boost/bind.hpp>
#include <boost/unordered_set.hpp>



#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/extern/metis/metis.hpp>

#include <graphlab/util/random.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab { 


  struct graph_partitioner {

    typedef uint32_t part_id_type;
      
    /**
       \brief the partition methods
    */
    enum partition_method {
      PARTITION_RANDOM, /**< Vertices are randomly assigned to partitions. Every partition
                           will have roughly the same number of vertices. (between
                           #vertices / nparts and (#vertices / nparts + 1)) */
      PARTITION_METIS,  /**< Partitions the graph using the
                           <a href="http://glaros.dtc.umn.edu/gkhome/views/metis"> METIS</a>
                           graph partitioning package. <b>A modified version of METIS is
                           included with GraphLab, so contact us if you encounter any
                           issues.</b>*/
      PARTITION_BFS, /**< Partitions the graph using Breadth First Searches. A random
                        root vertex is selected and inserted into the first
                        partition. Additional vertices are inserted into the partition
                        by starting a BFS from the root vertex.  When the partition
                        reaches capacity: (#vertices / nparts), the procedure restarts
                        with a next partition with a new random root vertex. */
      PARTITION_EDGE_NUM, /**< Partitions the vertices such that every partition 
                             has roughly the same number of edges. */
    };
    
    /// Converts a partition_method_enum to a string
    inline static std::string enum_to_string(partition_method val) {
      switch(val) {
      case PARTITION_RANDOM:
        return "random";
      case PARTITION_METIS:
        return "metis";
      case PARTITION_BFS:
        return "bfs";
      case PARTITION_EDGE_NUM:
        return "edge_num";
      default:
        return "";
      }
    }
    
    /// Converts a string to a partition_method_enum. Returns true on success
    inline static bool string_to_enum(std::string s, partition_method &val) {
      if (s == "random") {
        val = PARTITION_RANDOM;
        return true;
      }
      else if (s == "metis") {
        val = PARTITION_METIS;
        return true;
      }
      else if (s == "bfs") {
        val = PARTITION_BFS;
        return true;
      }
      else if (s == "edge_num") {
        val = PARTITION_EDGE_NUM;
        return true;
      }
      return false;
    }
  

  

    /**
     * \brief Randomly assign vertices to partitions.  This will assign
     * vertices evenly to each partition.
     * Equivalent to calling partition() with the 
     * partition_method::PARTITION_RANDOM parameter
     *
     * \param nparts The number of parts to partition into
     * \param[out] vertex2part A vector providing a vertex_id -> partition_id mapping
     */
    template <typename Graph>
    inline static void random_partition(const Graph& graph,
                                        const size_t nparts, 
                                        std::vector<part_id_type>& vertex2part) {
      vertex2part.resize(graph.num_vertices());
      for (size_t i = 0; i < graph.num_vertices(); ++i) 
        vertex2part[i] = part_id_type(i % nparts);
      random::shuffle(vertex2part.begin(), vertex2part.end());  
    }


    /**
     * \brief partition the graph with roughly the same number of edges for
     * each part. Equivalent to calling partition() with the 
     * partition_method::PARTITION_EDGE_NUM parameter.
     *
     * \param nparts The number of parts to partition into
     * \param[out] vertex2part A vector providing a vertex_id -> partition_id mapping
     */
    template <typename Graph>
    inline static void edge_num_partition(const Graph& graph,
                            const size_t nparts, 
                            std::vector<part_id_type>& vertex2part){
      vertex2part.resize(graph.num_vertices());
      const size_t num_edges = 2 * graph.num_edges();
      const size_t edge_per_part = num_edges / nparts;
      part_id_type curpart = 0;
      std::vector<part_id_type> parts;
      parts.resize(nparts, 0);      
      for (size_t i = 0; i< graph.num_vertices(); i++){
        const part_id_type ne = 
          part_id_type(graph.out_edges(i).size() + graph.in_edges(i).size());
        vertex2part[i] = curpart;
        parts[curpart]+= ne;           
        if (parts[curpart] >= edge_per_part  && curpart < nparts-1){
          curpart++;
        }
      }
    }


   


    /**
     * \brief Use a modified version of the METIS library to partition the
     * graph. Equivalent to calling partition() with the 
     * partition_method::PARTITION_METIS paramter
     *
     * \param numparts The number of parts to partition into
     * \param[out] ret_part A vector providing a vertex_id -> partition_id mapping
     *
     * The metis library is described in:
     *
     *   "A Fast and Highly Quality Multilevel Scheme for Partitioning
     *   Irregular Graphs”. George Karypis and Vipin Kumar. SIAM
     *   Journal on Scientific Computing, Vol. 20, No. 1, pp. 359—392,
     *   1999.
     *
     * We have modified an alpha version (5.0) to work with the
     * GraphLab framework.  Therefore users having trouble with this
     * function or the included Metis source should direct concerns to
     * the contact information provided at:
     *
     *   http://www.select.cs.cmu.edu/code
     */
    template <typename Graph>
    inline static void metis_partition(const Graph& graph,
                                       size_t numparts, 
                                       std::vector<part_id_type>& ret_part) {
      typedef typename Graph::vertex_id_type vertex_id_type;
      typedef typename Graph::edge_type   edge_type;

      if (numparts == 1) {
        ret_part.assign(graph.num_vertices(), 0);
        return;
      }
      // Determine parameters needed to construct the partitioning
      metis::idxtype numverts(graph.num_vertices());
      ASSERT_GT(numverts, 0);
      // Compute the number of edges 
      metis::idxtype numedges(graph.num_edges());

      // allocate metis data structures
      metis::idxtype* vweight = new metis::idxtype[numverts];
      ASSERT_NE(vweight, NULL);    
      metis::idxtype* xadj = new metis::idxtype[numverts + 1];
      ASSERT_NE(xadj, NULL);
      metis::idxtype* adjacency = new metis::idxtype[2 * numedges];
      ASSERT_NE(adjacency, NULL);
      metis::idxtype* eweight = NULL;
      //       if(weighted) {
      //         eweight = new idxtype[numedges];
      //         assert(eweigth != NULL);
      //       }
      metis::idxtype* res = new metis::idxtype[numverts];   
      ASSERT_NE(res, NULL);

      // Pass through vertices filling in the metis data structures
      size_t offset = 0;
      for(size_t u = 0; u < graph.num_vertices(); ++u) {
        // Update vertex weight
        // Set weight
        vweight[u] = 1;
        // Update the offset
        xadj[u] = offset;
        // Fill the the adjacency data
      
        std::set<vertex_id_type> neighbors;
        foreach(const edge_type& edge, graph.out_edges(u)) {
          neighbors.insert(edge.target());
        }
        foreach(const edge_type& edge, graph.in_edges(u)) {
          neighbors.insert(edge.source());
        }
        foreach(vertex_id_type vid, neighbors) {
          if (vid == u) continue;
          adjacency[offset] = vid;
          ASSERT_GE(adjacency[offset], 0);
          offset++;
          ASSERT_GE(offset, 0);
        }
    
      } // end of data structure creation
      
      // Set the last entry in xadj to the end of the adjacency array
      xadj[numverts] = offset;
    
      // Set additional metis flags
      /**
       * 0 No weights (vwgts and adjwgt are NULL) 
       * 1 Weights on the edges only (vwgts = NULL) 
       * 2 Weights on the vertices only (adjwgt = NULL) 
       * 3 Weights both on vertices and edges. 
       */
      metis::idxtype weightflag = 2;
      // 0 for C-style numbering starting at 0 (1 for fortran style)
      metis::idxtype numflag = 0;
      // the number of parts to cut into
      metis::idxtype nparts = numparts;     
      // Options array (only care about first element if first element
      // is zero
      metis::idxtype options[5] = {0}; 
      // output argument number of edges cut
      metis::idxtype edgecut = 0;
    
      // Call kmetis
      metis::METIS_PartGraphKway(&(numverts), 
                                 xadj,
                                 adjacency,
                                 vweight,
                                 eweight,
                                 &(weightflag),
                                 &(numflag),
                                 &(nparts),
                                 options,
                                 &(edgecut),
                                 res);
    
      //     // Call pmetis
      //     metis::METIS_PartGraphRecursive(&(numverts), 
      //                                     xadj,
      //                                     adjacency,
      //                                     vweight,
      //                                     eweight,
      //                                     &(weightflag),
      //                                     &(numflag),
      //                                     &(nparts),
      //                                     options,
      //                                     &(edgecut),
      //                                     res);
    
      // destroy all unecessary data structures except res
      if(xadj != NULL) delete [] xadj;
      if(adjacency != NULL) delete [] adjacency;
      if(vweight != NULL) delete [] vweight;
      if(eweight != NULL) delete [] eweight;

      // Resize the partition
      ret_part.resize(graph.num_vertices());
      // process the final results
      ASSERT_NE(res, NULL);
      for(vertex_id_type v = 0; v < graph.num_vertices(); ++v) {
        ret_part[v] = (vertex_id_type)res[v];
      }    
      // Delete the result array
      if(res != NULL) delete [] res;
    } // end of metis partition




    /**
     * This function computes a weighted graph partition using METIS.
     * \param numparts The number of parts to partition into
     * \param[out] ret_part A vector providing a vertex_id -> partition_id mapping
     * \param vfunction A function of the type size_t (*)(const VertexData &v)
     *                  This function takes in the data on a vertex, and 
     *                  returns the weight of the vertex
     * \param wfunction A function of the type size_t (*)(const EdgeData &e)
     *                  This function takes in the data on an edge, and 
     *                  returns the weight of the edge
     * \param usemetisdefaults If set to true, uses Metis default parameter set.
     *                         defaults to false.
     *                         
     *
     * Use a modified version of the METIS library to partition the
     * graph using user provided edge and vertex weight functions.
     * The methis library is described in:
     *
     *   "A Fast and Highly Quality Multilevel Scheme for Partitioning
     *   Irregular Graphs”. George Karypis and Vipin Kumar. SIAM
     *   Journal on Scientific Computing, Vol. 20, No. 1, pp. 359—392,
     *   1999.
     *
     * We have modified an alpha version (5.0) to work with the
     * GraphLab framework.  Therefore users having trouble with this
     * function or the included Metis source should direct concerns to
     * the contact information provided at:
     *
     *   http://www.select.cs.cmu.edu/code
     *
     */
    template <typename Graph,
              typename EdgeWeightFunction, 
              typename VertexWeightFunction>
    inline static void metis_weighted_partition(const Graph& graph,
                                                const size_t numparts,
                                                std::vector<part_id_type>& ret_part,
                                                VertexWeightFunction vfunction,
                                                EdgeWeightFunction wfunction,
                                                bool usemetisdefaults = false) {
      typedef typename Graph::vertex_id_type vertex_id_type;
      typedef typename Graph::edge_id_type edge_id_type;
      typedef typename Graph::edge_type edge_type;
      if (numparts == 1) {
        ret_part.assign(graph.num_vertices(), 0);
        return;
      }
      // Determine parameters needed to construct the partitioning
      metis::idxtype numverts(graph.num_vertices());
      ASSERT_GT(numverts, 0);
      // Compute the number of edges 
      metis::idxtype numedges (graph.num_edges());

      // allocate metis data structures
      metis::idxtype* vweight = new metis::idxtype[numverts];
      ASSERT_NE(vweight, NULL);
      metis::idxtype* xadj = new metis::idxtype[numverts + 1];
      ASSERT_NE(xadj, NULL);
      metis::idxtype* adjacency = new metis::idxtype[2 * numedges];
      ASSERT_NE(adjacency, NULL);
      metis::idxtype* eweight = NULL;
      eweight = new metis::idxtype[ 2 * numedges];
      ASSERT_NE(eweight, NULL);
      
      metis::idxtype* res = new metis::idxtype[numverts];   
      ASSERT_NE(res, NULL);

      // Pass through vertices filling in the metis data structures
      size_t offset = 0;
      for(size_t u = 0; u < graph.num_vertices(); ++u) {
        // Update vertex weight
        // Set weight
        vweight[u] = double(vfunction(graph.vertex_data(u)));
        // Update the offset 
        xadj[u] = offset;
        // Fill the the adjacency data
      
        std::set<size_t> neighbors;
        std::map<size_t, double> nbrtoweight;
        foreach(const edge_type& edge, graph.out_edges(u)) {
          neighbors.insert(edge.target());
          nbrtoweight[edge.target()] = double(wfunction(graph.edge_data(edge)));
        }
        foreach(const edge_type& edge, graph.in_edges(u)) {
          neighbors.insert(edge.source());
          nbrtoweight[edge.source()] = double(wfunction(graph.edge_data(edge)));
        }
        foreach(vertex_id_type vid, neighbors) {
          if (vid == u) continue;
          adjacency[offset] = vid;
          eweight[offset] = nbrtoweight[vid];
          ASSERT_GE(adjacency[offset], 0);
          offset++;
          ASSERT_GE(offset, 0);
        }
    
      } // end of data structure creation
      
      // Set the last entry in xadj to the end of the adjacency array
      xadj[numverts] = offset;
    
      // Set additional metis flags
      /**
       * 0 No weights (vwgts and adjwgt are NULL) 
       * 1 Weights on the edges only (vwgts = NULL) 
       * 2 Weights on the vertices only (adjwgt = NULL) 
       * 3 Weights both on vertices and edges. 
       */
      metis::idxtype weightflag = 3;
      // 0 for C-style numbering starting at 0 (1 for fortran style)
      metis::idxtype numflag = 0;
      // the number of parts to cut into
      metis::idxtype nparts = numparts;     
      // Options array (only care about first element if first element
      // is zero
      metis::idxtype options[5] = {0};
      
      options[0] = 1;
      options[1]=3;
      options[2]=1;
      options[3]=2;
      options[4]=0;
      if (usemetisdefaults) options[0] = 0;
      // output argument number of edges cut
      metis::idxtype edgecut = 0;
    
      // Call kmetis
      metis::METIS_PartGraphKway(&(numverts), 
                                 xadj,
                                 adjacency,
                                 vweight,
                                 eweight,
                                 &(weightflag),
                                 &(numflag),
                                 &(nparts),
                                 options,
                                 &(edgecut),
                                 res);
    
      // Call pmetis
      /*   metis::METIS_PartGraphRecursive(&(numverts), 
           xadj,
           adjacency,
           vweight,
           eweight,
           &(weightflag),
           &(numflag),
           &(nparts),
           options,
           &(edgecut),
           res);*/
    
      // destroy all unecessary data structures except res
      if(xadj != NULL) delete [] xadj;
      if(adjacency != NULL) delete [] adjacency;
      if(vweight != NULL) delete [] vweight;
      if(eweight != NULL) delete [] eweight;

      // Resize the partition
      ret_part.resize(graph.num_vertices());
      // process the final results
      ASSERT_NE(res, NULL);
      for(vertex_id_type v = 0; v < graph.num_vertices(); ++v) {
        ret_part[v] = res[v];
      }    
      // Delete the result array
      if(res != NULL) delete [] res;
    } // end of metis partition




    /**
     * \brief Performs a breadth first search partitioning of the graph.
     * Equivalent to calling partition() with the 
     * partition_method::PARTITION_EDGE_NUM paramter
     *
     * \param nparts The number of parts to partition into
     * \param[out] vertex2part A vector providing a vertex_id -> partition_id mapping
     *
     * The algorithm works by picking up a random vertex and performing a breadth
     * first search until the number of vertices touched is |V|/nparts.
     * This will then be assigned as the first partition. This procedure repeats
     * until all partitions are filled.
     */
    template<typename Graph>
    inline static void bfs_partition(const Graph& graph,
                                     const size_t nparts, 
                                     std::vector<part_id_type>& vertex2part) {
      typedef typename Graph::vertex_id_type vertex_id_type;
      typedef typename Graph::edge_type   edge_type;
      // create a list of all unassigned variables
      std::set<vertex_id_type> unassigned;
      vertex2part.resize(graph.num_vertices());
      // initialize the unassigned vertices
      for(vertex_id_type v = 0; v < graph.num_vertices(); ++v) {
        unassigned.insert(v);
      }
      // Compute the partition size
      size_t maxpartsize = (size_t)(std::ceil(double(unassigned.size()) / (double)nparts));
      size_t partid = 0;
      while(!unassigned.empty()) {  
        std::list<vertex_id_type> queue;    // Breadth first queue 
        std::set<vertex_id_type>  visited;  // Set of visited vertices
        // While the task is still too small and their remains
        // unassigned vertices
        size_t curpartsize = 0;
        while(curpartsize < maxpartsize 
              && !unassigned.empty()) {
          if(queue.empty()) { 
            queue.push_front(*unassigned.begin());
            visited.insert(*unassigned.begin());
          }
          ASSERT_FALSE(queue.empty());
          // Pop the first element off the queue 
          vertex_id_type v = queue.front(); queue.pop_front();
          ASSERT_LT(partid, nparts);
          // Add the element to the task
          vertex2part[v] = (uint32_t)partid;
          ++curpartsize;
          // Remove the vertex from the set of unassigned vertices
          unassigned.erase(v); 
          // Add all its unassigned and unvisited neighbors to the queue
          foreach(const edge_type& edge, graph.out_edges(v)) {
            vertex_id_type u = edge.target();
            if(unassigned.find(u) != unassigned.end() &&
               visited.find(u) == visited.end()) {
              queue.push_back(u);
              visited.insert(u);
            }
          }
          foreach(const edge_type& edge, graph.in_edges(v)) {
            vertex_id_type u = edge.source();
            if(unassigned.find(u) != unassigned.end() &&
               visited.find(u) == visited.end()) {
              queue.push_back(u);
              visited.insert(u);
            }
          }
        } // End of block build foor loop
        // move to the next part
        partid++;
      }// end of outer while loop
    } // end of bfs partition


    /**
     * Partition the graph using one of the available partitioning
     * methods.
     * \param partmethod A partition method. \ref partition_method
     *
     * \param nparts The number of parts to partition into
     * \param[out] vertex2part A vector providing a vertex_id -> partition_id mapping
     * 
     */
    template<typename Graph>
    inline static void partition(const partition_method partmethod,
                                 const Graph& graph,
                                 const size_t nparts, 
                                 std::vector<part_id_type>& vertex2part) {
      switch (partmethod) {
      case PARTITION_METIS:
        return metis_partition(graph, nparts, vertex2part);
      case PARTITION_BFS:
        return bfs_partition(graph, nparts, vertex2part);
      case PARTITION_RANDOM:
        return random_partition(graph, nparts, vertex2part);
      case PARTITION_EDGE_NUM:
        return edge_num_partition(graph, nparts, vertex2part);
      default:
        ASSERT_TRUE(false); //shoud never ever happen
      }
    }



    /**
     * Partition the graph using one of the available partitioning
     * methods.
     * \param partmethod A partition method. \ref partition_method
     *
     * \param nparts The number of parts to partition into
     * \param[out] vertex2part A vector providing a vertex_id -> partition_id mapping
     * 
     */
    template<typename Graph>
    inline static void partition(const std::string& partition_method_str,
                                 const Graph& graph,
                                 const size_t nparts, 
                                 std::vector<part_id_type>& vertex2part) {

      // Parse the partitioning method string
      partition_method partmethod(PARTITION_RANDOM);
      const bool successful_parse =
        string_to_enum(partition_method_str, partmethod);
      if(!successful_parse) {
        logstream(LOG_FATAL) << "Invalid partitioning method string: "
                             << partition_method_str << std::endl;
      }

      switch (partmethod) {
      case PARTITION_METIS:
        return metis_partition(graph, nparts, vertex2part);
      case PARTITION_BFS:
        return bfs_partition(graph, nparts, vertex2part);
      case PARTITION_RANDOM:
        return random_partition(graph, nparts, vertex2part);
      case PARTITION_EDGE_NUM:
        return edge_num_partition(graph, nparts, vertex2part);
      default:
        ASSERT_TRUE(false); //shoud never ever happen
      }
    }



  }; // end of namespace graph partitioner
  
 

} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

