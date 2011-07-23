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

#ifndef GRAPHLAB_TYPES_HPP
#define GRAPHLAB_TYPES_HPP

#include <graphlab.hpp>


/**
   \namespace graphlab
   \brief The namespace containing all graphlab objects and functions.

   All objects and functions used in graphlab are contained within the
   graphlab namespace.  Forexample to access the graph type a user
   must therefore either include the graphlab namespace or use:

   <code>
     graphlab::graph<VertexType, EdgeType>
   </code>

   Because most of the graphlab types depend on the graph type we have
   created a templated struct called graphlab::types. \todo finish
   explanation.
                           
*/
namespace graphlab {
 
  //! Predecleration of the core type
  template <typename Graph, typename UpdateFunctor> class core;

  /**
  \brief A types datastructure which provides convenient
  specializations of all user-facing GraphLab types.
  
  GraphLab is heavily templatized. The graphlab::types object provides
  a convenient way to access the GraphLab classes without requiring
  excessive angle brackets (< , >). The GraphLab types object is
  located in <graphlab.hpp>.  To define a graphlab type object:
  
  \code
  typedef graphlab::graph<vertex_data, edge_data> graph_type;
  typedef graphlab::types<graph_type> gl;
  \endcode
  
  Now we can use gl::... to access all the available graphlab types. 
  */
  template<typename Graph, typename UpdateFunctor>
  struct types {
    ///  \brief The type of the Graph. 
    typedef Graph graph;
    typedef UpdateFunctor update_functor;

    /// \brief The type of the data stored on each vertex of the Graph. 
    typedef typename graph::vertex_data_type vertex_data;
    
    /// \brief The type of the data stored on each edge of the Graph.   
    typedef typename graph::edge_data_type   edge_data;

    
    typedef graphlab::disk_graph<vertex_data, edge_data> disk_graph;

    /** \brief A convenient wrapper object around the commonly used
    portions of GraphLab.  This is useful for most GraphLab
    applications. See the \ref graphlab::core object for more details.
    */
    typedef graphlab::core<graph, update_functor> core;



    typedef graphlab::command_line_options command_line_options;
    typedef graphlab::graphlab_options graphlab_options;
    
    
  
    typedef iupdate_functor<graph> iupdate_functor;
    
    typedef graphlab::iengine<graph, update_functor>    iengine;
    typedef graphlab::iscope<graph>                     iscope;
    typedef graphlab::ischeduler<iengine>               ischeduler;
    typedef graphlab::icallback<graph, update_functor>  icallback;
    
    //    typedef graphlab::imonitor<graph>            imonitor;

    typedef graphlab::glshared_sync_ops<Graph> glshared_sync_ops;
    typedef graphlab::glshared_apply_ops glshared_apply_ops;
    typedef graphlab::glshared_merge_ops glshared_merge_ops;


    
    typedef shared_memory_engine<graph, update_functor> shared_memory_engine;



    typedef graphlab::fifo_scheduler<graph> fifo_scheduler;

    // typedef graphlab::priority_scheduler<graph> priority_scheduler;
    // typedef graphlab::sampling_scheduler<graph> sampling_scheduler;
    // typedef graphlab::sweep_scheduler<graph> sweep_scheduler;
    // typedef graphlab::multiqueue_fifo_scheduler<graph> multiqueue_fifo_scheduler;
    // typedef graphlab::multiqueue_priority_scheduler<graph> 
    // multiqueue_priority_scheduler;
    // typedef graphlab::clustered_priority_scheduler<graph> clustered_priority_scheduler;
    // typedef graphlab::round_robin_scheduler<graph> round_robin_scheduler;
    // typedef graphlab::chromatic_scheduler<graph> chromatic_scheduler;
    
    
    
    
    

    /**
     * \brief The type of id assigned to each vertex. Equivalent to
     * graphlab::vertex_id_type 
     */
    typedef typename graph::vertex_id_type vertex_id;    


    
    /**
     * \brief The color type associated with each vertex.
     */
    typedef typename graph::vertex_color_type vertex_color;    


    /**
     * \brief The type of id assigned to each vertex. Equivalent to
     * graphlab::edge_id_type
     */
    typedef typename graph::edge_id_type edge_id;


    /**
     * \brief The type of id assigned to each vertex. Equivalent to
     * graphlab::edge_id_t
     */
    typedef typename graph::edge_list_type edge_list;

    
    typedef graphlab::options_map          scheduler_options;
    typedef graphlab::sched_status         sched_status;
    typedef graphlab::consistency_model    consistency_model;

    template <typename T>
    class glshared : public graphlab::glshared<T> { };

    template <typename T>
    class glshared_const : public graphlab::glshared_const<T> { };
  }; // end of types

}; // end of graphlab namespace



#endif

