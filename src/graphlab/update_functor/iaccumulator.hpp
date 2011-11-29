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
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved. 
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



#ifndef GRAPHLAB_IACCUMULATOR_HPP
#define GRAPHLAB_IACCUMULATOR_HPP



#include <graphlab/context/icontext.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {





  template<typename Graph, typename UpdateFunctor, typename Accumulator>
  class iaccumulator {    
  public:
    typedef Graph graph_type;
    typedef UpdateFunctor update_functor_type;
    typedef Accumulator   accumulator_type;

    typedef typename graph_type::vertex_data_type  vertex_data_type;
    typedef typename graph_type::vertex_id_type    vertex_id_type;   
    typedef typename graph_type::vertex_color_type vertex_color_type;

    typedef typename graph_type::edge_data_type    edge_data_type;
    typedef typename graph_type::edge_type         edge_type;
    typedef typename graph_type::edge_list_type    edge_list_type;

    typedef icontext<graph_type, update_functor_type> icontext_type;
    typedef iglobal_context iglobal_context_type;

    enum edge_set {IN_EDGES, OUT_EDGES, ALL_EDGES, NO_EDGES};

    
    virtual ~iaccumulator() { }

    /**
     * When multiple update functors are scheduled to be run on the
     * same function they are added. The default behavior is to simply
     * ignore the later update functors.
     */
    virtual void operator+=(const accumulator_type& other) = 0;

    /**
     * The main part of an update functor
     */
    virtual void operator()(icontext_type& context) = 0;


    /**
     * The finalize routine is called on the final accumulator after
     * accumulating over the entire graph
     */
    virtual void finalize(iglobal_context& context) = 0;

    /**
     * Gets the context range required by this update functor.  If not
     * implemented by the derived class then the default context range
     * is returned.
     */
    virtual consistency_model::model_enum consistency() const {
      return consistency_model::USE_DEFAULT;
    }


    virtual bool is_factorizable() const { return false; }

    virtual bool writable_gather() const { return false; }
    virtual bool writable_scatter() const { return true; }
    
    virtual edge_set gather_edges() const { return IN_EDGES; }
    virtual edge_set scatter_edges() const { return OUT_EDGES; }
    
    virtual void gather(icontext_type& context, const edge_type& edge) { };
    virtual void apply(icontext_type& context) { };
    virtual void scatter(icontext_type& context, const edge_type& edge) { };
  
  };  // end of iaccumulator
 
}; //end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
