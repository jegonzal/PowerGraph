/**
   \mainpage 
   
   \section intro_sec Introduction
   
   GraphLab is a powerful new system for designing and implementing
   parallel algorithms in machine learning.  While the current targets
   multi-core shared memory parallel systems we are in the process of
   implementing a distributed version and plan to provide support for
   alternative parallel architectures include GPUs in the near future.
 
   For a more user friendly tour of GraphLab and its features visit
   the site: <a href="http://www.graphlab.ml.cmu.edu/details.html">
   http://www.graphlab.ml.cmu.edu/details.html </a>
   
*/


#ifndef GRAPHLAB_MASTER_INCLUDES
#define GRAPHLAB_MASTER_INCLUDES


//
#include <graphlab/engine/engine_includes.hpp>
#include <graphlab/factors/factor_includes.hpp>
#include <graphlab/graph/graph_includes.hpp>
#include <graphlab/logger/logger_includes.hpp>
#include <graphlab/monitoring/monitoring_includes.hpp>
#include <graphlab/parallel/parallel_includes.hpp>
#include <graphlab/schedulers/scheduler_includes.hpp>
#include <graphlab/scope/scope_includes.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/shared_data/shared_data_includes.hpp>
#include <graphlab/tasks/task_includes.hpp>
#include <graphlab/util/util_includes.hpp>

#include <graphlab/core.hpp>

#ifdef GLDISTRIBUTED
#include <graphlab/distributed/distributed_includes.hpp>
#include <graphlab/distributed/distributed_engine.hpp>
#endif 



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
 
  
  template<typename Graph>
  struct types {
    typedef Graph graph;

    typedef graphlab::core<typename graph::vertex_data_type,
                           typename graph::edge_data_type> core;

    typedef graphlab::command_line_options command_line_options;
    typedef graphlab::engine_options engine_options;
    
    typedef typename graph::vertex_data_type vertex_data_type;
    typedef typename graph::edge_data_type   edge_data_type;
    
    typedef graphlab::update_task<graph>        update_task;
    typedef typename update_task::update_function_type update_function;
    
    typedef graphlab::iscope<graph>              iscope;
    typedef graphlab::ischeduler<graph>          ischeduler;
    typedef graphlab::icallback<graph>           icallback;
    typedef graphlab::iengine<graph>             iengine;
    typedef graphlab::imonitor<graph>            imonitor;

    typedef graphlab::ishared_data<graph>        ishared_data;
    typedef graphlab::ishared_data_manager<graph> ishared_data_manager;
    typedef graphlab::sync_ops<Graph> sync_ops;
    typedef graphlab::apply_ops<Graph> apply_ops;

    typedef graphlab::thread_shared_data<graph>  thread_shared_data;
    
    typedef graphlab::ivertex_set<graph>         ivertex_set;
    typedef graphlab::vertex_set<graph>          vset;
    typedef graphlab::restricted_vertex_set<graph> rvset;
    typedef typename rvset::selector_function_type  selector_function;
    typedef typename graphlab::execution_plan<graph> execution_plan;
    


    template<typename Scheduler, typename ScopeFactory>
    struct engines {
      typedef graphlab::synchronous_engine<graph> synchronous;
      typedef graphlab::
      asynchronous_engine<graph, Scheduler, ScopeFactory> asynchronous;
      #ifdef GLDISTRIBUTED
      typedef graphlab::distributed_engine<graph, Scheduler> distributed;
      #endif
    };
    

    typedef graphlab::fifo_scheduler<graph> fifo_scheduler;
    typedef graphlab::priority_scheduler<graph> priority_scheduler;
    typedef graphlab::sampling_scheduler<graph> sampling_scheduler;
    typedef graphlab::splash_scheduler<graph> splash_scheduler;
    typedef graphlab::sweep_scheduler<graph> sweep_scheduler;
    typedef graphlab::multiqueue_fifo_scheduler<graph> multiqueue_fifo_scheduler;
    typedef graphlab::multiqueue_priority_scheduler<graph> multiqueue_priority_scheduler;
    typedef graphlab::set_scheduler<graph> set_scheduler;
    typedef graphlab::clustered_priority_scheduler<graph> clustered_priority_scheduler;
    typedef graphlab::round_robin_scheduler<graph> round_robin_scheduler;
    typedef graphlab::colored_scheduler<graph> colored_scheduler;
    
    
    
    
    

    
    typedef graphlab::vertex_id_t vertex_id_t;
    typedef graphlab::edge_id_t edge_id_t;
    typedef graphlab::edge_list edge_list;
    
    typedef graphlab::scheduler_options          scheduler_options;
    typedef graphlab::sched_status               sched_status;
    typedef graphlab::partition_method           partition_method;
    typedef graphlab::scope_range scope_range;

    typedef graphlab::random  random;
  };
  
  // typedef types<blob_graph> blob_types;

}


#endif
