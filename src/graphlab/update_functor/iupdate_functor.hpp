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


#include <graphlab/scope/iscope.hpp>
#include <graphlab/engine/callback/icallback.hpp>

namespace graphlab {


  template<typename Graph> 
  class iupdate_functor {    
  public:
    typedef Graph graph_type;
    typedef typename graph_type::vertex_data_type vertex_data_type;
    typedef typename graph_Type::edge_data_type   edge_data_type;
    typedef typename iscope<graph_type> iscope_type;
    typedef typename icallback<Graph, iupdate_functor> icallback_type;

    
    /**
     * The main part of an update functor
     */
    virtual void operator()(iscope_type& scope, icallback_type& callback) = 0;


    /**
     * Gets the scope range required by this update functor.  If not
     * implemented by the derived class then the default scope range
     * is returned.
     */
    virtual scope_range::scope_range_enum scope() const {
      return scope_range::USE_DEFAULT;
    }

    /**
     * When multiple update functors are scheduled to be run on the
     * same function they are added. The default behavior is to simply
     * ignore the later update functors.
     */
    virtual void operator+=(const iupdate_functor& other) const { }

    /**
     * Get the priority of the update functor
     */
    virtual double priority() const { return double(1.0); }        
  }; 

}; //end of namespace graphlab

#endif
