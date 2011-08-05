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
//#include <graphlab/scope/iscope_manager.hpp>



namespace graphlab {
  template<typename Graph> class general_scope_factory;

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

    consistency_model::model_enum stype;
    //    iscope_manager<Graph>* manager;

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

