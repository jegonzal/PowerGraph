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
#include <graphlab/scope/iscope.hpp>
#include <graphlab/engine/callback/icallback.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  // /**
  //  * type deduction
  //  */
  // template<typename G, typename U, template T<G, U> > 
  // struct is_factorizable : public boost::false_type { }; 
  // template<typename G, typename U, template T< typename G, typename U> >
  // struct is_factorizable<G, U, iupdate_functor<G,U>::factorizble > : 
  //   public boost::true_type { }; 



  template<typename Graph, typename UpdateFunctor> 
  class iupdate_functor {    
  public:
    typedef Graph graph_type;
    typedef UpdateFunctor update_functor_type;

    typedef typename graph_type::vertex_data_type  vertex_data_type;
    typedef typename graph_type::vertex_id_type    vertex_id_type;   
    typedef typename graph_type::vertex_color_type vertex_color_type;

    typedef typename graph_type::edge_data_type    edge_data_type;
    typedef typename graph_type::edge_id_type      edge_id_type;
    typedef typename graph_type::edge_list_type    edge_list_type;

   
    typedef iscope<graph_type> iscope_type;
    typedef icallback<graph_type, update_functor_type> icallback_type;

    enum edge_set {IN_EDGES, OUT_EDGES, ALL_EDGES, NO_EDGES};
    
    virtual ~iupdate_functor() { }

    

    /**
     * The main part of an update functor
     */
    virtual void operator()(iscope_type& scope, 
                            icallback_type& callback) {
      // Gather
      if(gather_edges() == IN_EDGES || gather_edges() == ALL_EDGES) {
        foreach(const edge_id_type eid, scope.in_edge_ids()) 
          gather(scope, callback, eid);
      }
      if(gather_edges() == OUT_EDGES || gather_edges() == ALL_EDGES) {
        foreach(const edge_id_type eid, scope.out_edge_ids()) 
          gather(scope, callback, eid);
      }
      // Apply
      apply(scope, callback);
      // scatter
      if(scatter_edges() == IN_EDGES || scatter_edges() == ALL_EDGES) {
        foreach(const edge_id_type eid, scope.in_edge_ids()) 
          scatter(scope, callback, eid);
      }
      if(scatter_edges() == OUT_EDGES || scatter_edges() == ALL_EDGES) {
        foreach(const edge_id_type eid, scope.out_edge_ids()) 
          scatter(scope, callback, eid);
      }
    } // end of operator()

    virtual bool is_factorizable() const { return false; }
    
    virtual bool writable_gather() const { return false; }
    virtual bool writable_scatter() const { return true; }
    
    virtual edge_set gather_edges() const { return IN_EDGES; }
    virtual edge_set scatter_edges() const { return OUT_EDGES; }
    
    virtual void gather(iscope_type& scope, icallback_type& callback,
                        edge_id_type eid) { };
    virtual void apply(iscope_type& scope, icallback_type& callback) { };
    virtual void scatter(iscope_type& scope, icallback_type& callback,
                         edge_id_type eid) { };

    
        
    //    virtual void operator()(iscope_type& scope, icallback_type& callback) = 0;

    /**
     * Gets the scope range required by this update functor.  If not
     * implemented by the derived class then the default scope range
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
    virtual double priority() const { return double(1.0); }        
  }; 

}; //end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
