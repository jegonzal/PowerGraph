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



#ifndef GRAPHLAB_RW_LOCK_MANAGER_HPP
#define GRAPHLAB_RW_LOCK_MANAGER_HPP

#include <vector>

#include <graphlab/parallel/pthread_tools.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {


  template<typename Graph>
  class rw_lock_manager {

  public:
    typedef Graph graph_type;
    typedef typename graph_type::vertex_id_type  vertex_id_type;
    typedef typename graph_type::edge_list_type  edge_list_type;
    typedef typename edge_list_type::const_iterator iterator;
    typedef typename graph_type::edge_type       edge_type;

  private:
    const graph_type& graph;
    std::vector<rwlock> locks;

  public:
    rw_lock_manager(const graph_type& graph) : 
      graph(graph), locks(graph.num_vertices()) { }

    inline void writelock_vertex(const vertex_id_type& vid) {
      // ASSERT_LT(vid, locks.size()); 
      locks[vid].writelock();
    }

    inline void readlock_vertex(const vertex_id_type& vid) {
      // ASSERT_LT(vid, locks.size()); 
      locks[vid].readlock();
    }

    inline bool try_readlock_vertex(const vertex_id_type& vid) {
      // ASSERT_LT(vid, locks.size()); 
      return locks[vid].try_readlock();
    }

    inline void release_vertex(const vertex_id_type& vid) {
      // ASSERT_LT(vid, locks.size()); 
      locks[vid].unlock();
    }

    inline void lock_edge_context(const vertex_id_type& vid) {
      const edge_list_type in_edges  = graph.in_edges(vid);
      const edge_list_type out_edges = graph.out_edges(vid);
      // Loop over edges
      bool processed_center = false;
      iterator in = in_edges.begin(), out = out_edges.begin();
      while(in != in_edges.end() && out != out_edges.end()) {
        const vertex_id_type in_vid  = in->source();
        const vertex_id_type out_vid = out->target();
        if(__builtin_expect(!processed_center && 
                            vid < in_vid && vid < out_vid, false)) {
          processed_center = true; writelock_vertex(vid);
        } else if (in_vid == out_vid) {
          readlock_vertex(in_vid); ++in; ++out;
        } else if(in_vid < out_vid) {
          readlock_vertex(in_vid); ++in;
        } else {
          readlock_vertex(out_vid); ++out;
        }
      } // loop over both in and out
      // Finish locking any remaining in edges
      while(in != in_edges.end()) {
        const vertex_id_type in_vid  = in->source();        
        if(__builtin_expect(!processed_center && vid < in_vid, false)) {
          processed_center = true; writelock_vertex(vid);
        } else { readlock_vertex(in_vid); ++in; }
      } // loop over in edges
      // Finish locking any remaining out edges
      while(out != out_edges.end()) {
        const vertex_id_type out_vid  = out->target();
        if(__builtin_expect(!processed_center && vid < out_vid, false)) {
          processed_center = true; writelock_vertex(vid);
        } else { readlock_vertex(out_vid); ++out; }
      } // loop over out edges
      // If we still have not locked the center do so now
      if(__builtin_expect(!processed_center, false)) { writelock_vertex(vid); }
    } // end of lock edge context

    inline void lock_full_context(const vertex_id_type& vid) {
      const edge_list_type in_edges  = graph.in_edges(vid);
      const edge_list_type out_edges = graph.out_edges(vid);
      // Loop over edges
      bool processed_center = false;
      iterator in = in_edges.begin(), out = out_edges.begin();     
      while(in != in_edges.end() && out != out_edges.end()) {
        const vertex_id_type in_vid  = in->source();
        const vertex_id_type out_vid = out->target();
        if(__builtin_expect(!processed_center && 
                            vid < in_vid && vid < out_vid, false)) {
          processed_center = true; writelock_vertex(vid);
        } else if (in_vid == out_vid) {
          writelock_vertex(in_vid); ++in; ++out;
        } else if(in_vid < out_vid) {
          writelock_vertex(in_vid); ++in;
        } else { 
          writelock_vertex(out_vid); ++out;
        }
      } // loop over both in and out
      // Finish locking any remaining in edges
      while(in != in_edges.end()) {
        const vertex_id_type in_vid  = in->source();        
        if(__builtin_expect(!processed_center && vid < in_vid, false)) {
          processed_center = true; writelock_vertex(vid);
        } else { writelock_vertex(in_vid); ++in; }
      } // loop over in edges
      // Finish locking any remaining out edges
      while(out != out_edges.end()) {
        const vertex_id_type out_vid  = out->target();
        if(__builtin_expect(!processed_center && vid < out_vid, false)) {
          processed_center = true; writelock_vertex(vid);
        } else { writelock_vertex(out_vid); ++out; }
      } // loop over out edges
      // If we still have not locked the center do so now
      if(__builtin_expect(!processed_center, false)) { writelock_vertex(vid); }
    } // end of lock full context


    inline void release_context(const vertex_id_type& vid) {
      const edge_list_type in_edges  = graph.in_edges(vid);
      const edge_list_type out_edges = graph.out_edges(vid);
      // Loop over edges
      bool processed_center = false;
      iterator in = in_edges.begin(), out = out_edges.begin();     
      while(in != in_edges.end() && out != out_edges.end()) {
        const vertex_id_type in_vid  = in->source();
        const vertex_id_type out_vid = out->target();
        if(__builtin_expect(!processed_center && 
                            vid < in_vid && vid < out_vid, false)) {
          processed_center = true; release_vertex(vid);
        } else if (in_vid == out_vid) {
          release_vertex(in_vid); ++in; ++out;
        } else if(in_vid < out_vid) {
          release_vertex(in_vid); ++in;
        } else { 
          release_vertex(out_vid); ++out;
        }
      } // loop over both in and out
      // Finish releasing any remaining in edges
      while(in != in_edges.end()) {
        const vertex_id_type in_vid  = in->source();        
        if(__builtin_expect(!processed_center && vid < in_vid, false)) {
          processed_center = true; release_vertex(vid);
        } else { release_vertex(in_vid); ++in; }
      } // loop over in edges
      // Finish releasing any remaining out edges
      while(out != out_edges.end()) {
        const vertex_id_type out_vid  = out->target();
        if(__builtin_expect(!processed_center && vid < out_vid, false)) {
          processed_center = true; release_vertex(vid);
        } else { release_vertex(out_vid); ++out; }
      } // loop over out edges
      // If we still have not released the center do so now
      if(__builtin_expect(!processed_center, false)) { release_vertex(vid); }
    } // end of release context



    inline void lock_single_edge(vertex_id_type center, 
                                 vertex_id_type neighbor) {
      if(center < neighbor) { writelock_vertex(center); readlock_vertex(neighbor); }
      else {  readlock_vertex(neighbor); writelock_vertex(center); }
    }

    inline void release_single_edge(vertex_id_type center, 
                                    vertex_id_type neighbor) {
      release_vertex(center); release_vertex(neighbor);
    }

    inline void swap_single_edge(vertex_id_type center, 
                                 vertex_id_type old_neighbor,
                                 vertex_id_type new_neighbor) {
      release_vertex(old_neighbor);
      if(!try_readlock_vertex(new_neighbor)) {
        release_vertex(center);
        lock_single_edge(center, new_neighbor);
      }      
    } // end of swap_single_edge

    
  }; // end of rw_lock_manager


}; // end of graphlab namespace
#include <graphlab/macros_undef.hpp>



#endif

