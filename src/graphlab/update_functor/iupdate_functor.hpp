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



#ifndef GRAPHLAB_IUPDATE_FUNCTOR_HPP
#define GRAPHLAB_IUPDATE_FUNCTOR_HPP


// #include <boost/type_traits.hpp>
#include <graphlab/context/icontext.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {



  
  /**
   * This interface should be used as the base class of all user
   * defined update functors.
   *
   *  \todo: Provide detailed explanation of standard usage pattern.
   */
  template<typename Graph, typename UpdateFunctor> 
  class iupdate_functor {    
  public:
    typedef Graph graph_type;
    typedef UpdateFunctor update_functor_type;

    typedef typename graph_type::vertex_data_type  vertex_data_type;
    typedef typename graph_type::vertex_id_type    vertex_id_type;   
    typedef typename graph_type::vertex_color_type vertex_color_type;

    typedef typename graph_type::edge_data_type    edge_data_type;
    typedef typename graph_type::edge_type         edge_type;
    typedef typename graph_type::edge_list_type    edge_list_type;

   
    typedef icontext<graph_type, update_functor_type> icontext_type;
    typedef iglobal_context iglobal_context_type;

    /**
     * The set of edges that are operated on during gather and scatter
     * operations.
     */
    enum edge_set {IN_EDGES, OUT_EDGES, ALL_EDGES, NO_EDGES};
    
    virtual ~iupdate_functor() { }

    /**
     * Gets the context range required by this update functor.  If not
     * implemented by the derived class then the default context range
     * is returned.
     */
    virtual consistency_model::model_enum consistency() const {
      return consistency_model::USE_DEFAULT;
    }

    /**
     * When multiple update functors are scheduled to be run on the
     * same function they are added. The default behavior is to simply
     * ignore the later update functors.
     */
    virtual void operator+=(const update_functor_type& other) const { }

    /**
     * Get the priority of the update functor
     */
    virtual double priority() const { return double(0); }        

    /**
     * The main part of an update functor
     */
    virtual void operator()(icontext_type& context) { 
      logstream(LOG_FATAL) 
        << "Operator() not implemented!" << std::endl;
    } 

    /**
     * Returns true if the factorized (gather, apply, scatter) version
     * of the update functor is to be used.
     */
    virtual bool is_factorizable() const { return false; }

    /**
     * Returns true of the adjacent edge and vertex are modified
     * during gather.
     */
    virtual bool writable_gather() const { return false; }

    /**
     * Returns true of the adjacent edge and vertex are modified
     * during the gather.
     */
    virtual bool writable_scatter() const { return true; }
    
    /**
     * Returns the set of edges to gather 
     */
    virtual edge_set gather_edges() const { return IN_EDGES; }

    /**
     * Returns the set of edges to scatter
     */
    virtual edge_set scatter_edges() const { return OUT_EDGES; }
    
    /**
     * Init gather is called before gathering
     */
    virtual void init_gather(iglobal_context_type& context) { };

    /**
     * Gather is called on all gather_edges() and may be called in
     * parallel.  The merge() operation is used to join update
     * functors.
     */
    virtual void gather(icontext_type& context, const edge_type& edge) { 
      logstream(LOG_FATAL) << "Gather not implemented!" << std::endl;
    };

    /**
     * Merges update functors during the gather process.
     */
    virtual void merge(const update_functor_type& other) {
      logstream(LOG_FATAL) << "Gather not implemented!" << std::endl;
    }

    /**
     * Apply is called within the vertex consistency model on the
     * center vertex after all gathers have completed.
     */
    virtual void apply(icontext_type& context) { 
      logstream(LOG_FATAL) << "Apply not implemented!" << std::endl;
    };
    
    
    /**
     * Scatter is invoked on all scatter_edges() after calling
     * init_scatter() and may be called in parallel.
     */
    virtual void scatter(icontext_type& context, const edge_type& edge) { 
      logstream(LOG_FATAL) << "Scatter not implemented!" << std::endl;
    }
  };  // end of iupdate_functor
 
}; //end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
