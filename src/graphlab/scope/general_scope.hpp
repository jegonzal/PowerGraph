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


#ifndef GRAPHLAB_GENERAL_SCOPE_HPP
#define GRAPHLAB_GENERAL_SCOPE_HPP

//#include <boost/bind.hpp>

#include <graphlab/scope/iscope.hpp>

#include <graphlab/scope/idiffable.hpp>


namespace graphlab {



  /**
   * This defines a general scope type 
   */
  template<typename Graph>
  class general_scope :
    public iscope<Graph> {
  public:
    typedef iscope<Graph> base;
    typedef typename base::vertex_id_type   vertex_id_type;
    typedef typename base::edge_id_type     edge_id_type;
    typedef typename base::vertex_data_type vertex_data_type;
    typedef typename base::edge_data_type   edge_data_type;
    typedef typename base::edge_list_type   edge_list_type;
    

    using base::_vertex;
    using base::_graph_ptr;

   


    // Cache related members --------------------------------------------------
    struct cache_entry {
      vertex_data_type old;
      vertex_data_type current;
      uint16_t reads, writes;
      cache_entry() : reads(0), writes(0) { }
    };
    typedef std::map<vertex_id_type, cache_entry> cache_map_type;


    edge_id_type eid;
    cache_map_type cache;

    


  public:

    general_scope() :
      base(NULL,NULL) { }

    general_scope(Graph* graph_ptr, 
                  vertex_id_type vertex,
                  consistency_model::model_enum stype 
                  = consistency_model::EDGE_CONSISTENCY) :
      base(graph_ptr, vertex, stype) { }

    

    void init(Graph* graph, vertex_id_type vertex, 
              consistency_model::model_enum consistency) {
      base::_graph_ptr = graph;
      base::_vertex = vertex;
      base::_consistency = consistency;
    } // end of init

    edge_id_type& edge_id() { return eid; }

    vertex_data_type& vertex_data(const vertex_id_type vid) {
      typedef typename cache_map_type::iterator iterator_type;
      iterator_type iter = cache.find(vid);
      if(iter != cache.end()) {
        return iter->second.current;
      } else {
        return _graph_ptr->vertex_data(vid);
      }
    }

    const vertex_data_type& vertex_data(const vertex_id_type vid) const {
      typedef typename cache_map_type::const_iterator iterator_type;
      iterator_type iter = cache.find(vid);
      if(iter != cache.end()) {
        return iter->second.current;
      } else {
        return _graph_ptr->vertex_data(vid);
      }
    }



    vertex_data_type& vertex_data() {
      return vertex_data(_vertex);
    }

    const vertex_data_type& vertex_data() const {
      return vertex_data(_vertex);
    }

    const vertex_data_type& const_vertex_data() const {
      return vertex_data(_vertex);
    }

    vertex_data_type& neighbor_vertex_data(vertex_id_type vid) {
      return vertex_data(vid);
    }

    const vertex_data_type& neighbor_vertex_data(vertex_id_type vid) const {
      return vertex_data(vid);
    }

    const vertex_data_type& const_neighbor_vertex_data(vertex_id_type vid) const {
      return vertex_data(vid);
    }

    
    // bool experimental_scope_upgrade(consistency_model::model_enum newrange) { 
    //   assert(manager != NULL);
    //   manager->release_scope(this);

    //   ASSERT_TRUE(manager->get_scope(thread::thread_id(), 
    //                                  base::_vertex,newrange) == this);
    //   stype = newrange;
    //   return true;
    // }

    //friend class scope_manager<Graph>;
  }; // end of ext_locked_scope

} // end of graphlab namespace


#endif

