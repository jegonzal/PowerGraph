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


#ifndef GRAPHLAB_DISTRIBUTED_AGGREGATOR
#define GRAPHLAB_DISTRIBUTED_AGGREGATOR

#ifndef __NO_OPENMP__
#include <omp.h>
#endif

#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/vertex_program/icontext.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/util/generics/conditional_addition_wrapper.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * \internal
   * Implements a distributed aggregator interface which can be plugged
   * into the engine. This class includes management of periodic aggregators.
   * 
   * Essentially, the engine should ideally pass-through all calls to
   *  - add_vertex_aggregator()
   *  - add_edge_aggregator()
   *  - aggregate_now()
   *  - aggregate_periodic()
   * 
   * On engine start(), the engine should call aggregate_all_periodic() 
   * to ensure all periodic aggregators are called once prior to vertex program
   * execution. After which, the start() function should be called to prepare
   * the state of the schedule. At termination of the engine, the stop()
   * function should be called to reset the state of the aggregator.
   * 
   * During engine execution, two modes of operations are permitted: 
   * synchronous, and asynchronous. In a synchronous mode of execution,
   * the tick_synchronous() function should be called periodically by 
   * exactly one thread on each machine, at the same time. In an asynchronous
   * mode of execution, tick_asynchronous() should be called periodically
   * on each machine by some arbitrary thread. This polls the state of the 
   * schedule and activates aggregation jobs which are ready. 
   * 
   * tick_synchronous() and tick_asynchronous() should not be used 
   * simultaneously within the same engine execution . For details on their 
   * usage, see their respective documentation.
   * 
   */
  template<typename Graph>
  class distributed_aggregator {
  public:
    typedef Graph graph_type;
    typedef typename graph_type::local_edge_list_type local_edge_list_type;
    typedef typename graph_type::local_edge_type local_edge_type;
    typedef typename graph_type::edge_type edge_type;
    typedef typename graph_type::local_vertex_type local_vertex_type;
    typedef typename graph_type::vertex_type vertex_type ;

    dc_dist_object<distributed_aggregator> rmi;
    graph_type& graph;
    icontext_type* context;
    
  private:
    
    /**
     * \internal
     * A base class which contains a "type-free" specification of the reduction
     * operation, thus allowing the aggregation to be performs at runtime
     * with no other type information whatsoever.
     */
    struct imap_reduce_base {
      /** \brief makes a copy of the current map reduce spec without copying 
       *         accumulator data     */
      virtual imap_reduce_base* clone_empty() const = 0;
      
      /** \brief Performs a map operation on the given vertex adding to the
       *         internal accumulator */
      virtual void perform_map_vertex(icontext_type&,
                                      vertex_type&) = 0;
                                      
      /** \brief Performs a map operation on the given edge adding to the
       *         internal accumulator */
      virtual void perform_map_edge(icontext_type&,
                                    edge_type&) = 0;
                                    
      /** \brief Returns true if the accumulation is over vertices. 
                 Returns false if it is over edges.*/
      virtual bool is_vertex_map() const = 0;      
      
      /** \brief Returns the accumulator stored in an any. 
                 (by some magic, any's can be serialized) */
      virtual any get_accumulator() const = 0;
      
      /** \brief Combines accumulators using a second accumulator 
                 stored in an any (as returned by get_accumulator) */
      virtual void add_accumulator_any(any& other) = 0;
      
      /** \brief Combines accumulators using a second accumulator 
                 stored in a second imap_reduce_base class) */
      virtual void add_accumulator(imap_reduce_base* other) = 0;
      
      /** \brief Resets the accumulator */
      virtual void clear_accumulator();
      
      /** \brief Calls the finalize operation on internal accumulator */
      virtual void finalize(icontext_type&) = 0;
    };
    
    /**
     * \internal
     * A templated implementation of the imap_reduce_base above.
     * \tparam ReductionType The reduction type. (The type the map function returns)
     */
    template <typename ReductionType>
    struct map_reduce_type : public imap_reduce_base {
      conditional_addition_wrapper<ReductionType> acc;
      boost::function<ReductionType(icontext_type&, vertex_type&)> map_vtx_function;
      boost::function<ReductionType(icontext_type&, edge_type&)> map_edge_function;
      boost::function<icontext_type&, const ReductionType&> finalize_function;
      
      bool vertex_map;
      
      /**
       * \brief Constructor which constructs a vertex reduction
       */
      map_reduce_type(
       boost::function<ReductionType(icontext_type&,vertex_type&)> map_vtx_function,
       boost::function<icontext_type&,  const ReductionType&> finalize_function)
          : map_vtx_function(map_vtx_function), finalize_function(finalize_function), vertex_map(true) { }

      /**
       * \brief Constructor which constructs an edge reduction. The last bool
       * is unused and allows for disambiguation between the two constructors
       */
      map_reduce_type(
       boost::function<ReductionType(icontext_type&, edge_type&)> map_edge_function,
       boost::function<icontext_type&,  const ReductionType&> finalize_function,
        bool)
          : map_edge_function(map_edge_function), finalize_function(finalize_function), vertex_map(false) { }


      void perform_map_vertex(icontext_type& context, vertex_type& vertex) {
        acc += map_vtx_function(context, vertex);
      }
      
      void perform_map_edge(icontext_type& context, edge_type& vertex) {
        acc += map_edge_function(context, vertex);
      }
      
      bool is_vertex_map() const {
        return vertex_map;
      }
      
      any get_accumulator() const {
        return acc;
      }
      
      void add_accumulator_any(any& other) {
        acc += other.as<conditional_addition_wrapper<ReductionType> >();
      }

      void add_accumulator(imap_reduce_base* other) {
        acc += dynamic_cast<map_reduce_type<Reduction_type>*>(other)->acc;
      }

      void clear_accumulator() {
        acc.clear();
      }

      void finalize(icontext_type& context) {
        finalize_function(acc.value, other);
      }
      
      imap_reduce_base* clone_empty() const {
        map_reduce_type<ReductionType>* copy;
        if (is_vertex_map) {
          copy = new map_reduce_type<Reduction_type(map_vtx_function, 
                                                    finalize_function);
        }
        else {
          copy = new map_reduce_type<Reduction_type(map_edge_function, 
                                                    finalize_function,
                                                    true);
        }
        return copy;
      }
    };
    

    std::map<std::string, imap_reduce_base*> aggregators;
    std::map<std::string, float> aggregate_period;
    
    float start_time;
    
    /* annoyingly the mutable queue is a max heap when I need a min-heap
     * to track the next thing to activate. So we need to keep 
     *  negative priorities... */
    mutable_queue<std::string, float> schedule;
    size_t ncpus;
  public:

    
    distributed_aggregator(distributed_control& dc, 
                           graph_type& graph, 
                           icontext_type* context):
                            rmi(dc, this), graph(graph), 
                            context(context), ncpus(0) { }

    /** \brief Creates a vertex aggregator. Returns true on success.
               Returns false if an aggregator of the same name already exists */
    template <typename ReductionType>
    bool add_vertex_aggregator(std::string key,
      boost::function<ReductionType(icontext_type&, vertex_type&)> map_function,
      boost::function<icontext_type&, const ReductionType&> finalize_function) {
      
      if (aggregators.count(key) == 0) {
        aggregators[key] = new map_reduce_type<ReductionType>(map_function, 
                                                              inalize_function);
        return true;
      }
      else {
        // aggregator already exists. fail 
        return false;
      }
    }
    
    
    /** \brief Creates a edge aggregator. Returns true on success.
               Returns false if an aggregator of the same name already exists */
    template <typename ReductionType>
    bool add_edge_aggregator(std::string key,
      boost::function<ReductionType(icontext_type&, edge_type&)> map_function,
      boost::function<icontext_type&, const ReductionType&> finalize_function) {
      
      if (aggregators.count(key) == 0) {
        aggregators[key] = new map_reduce_type<ReductionType>(map_function, 
                                                              finalize_function, 
                                                              true);
        return true;
      }
      else {
        // aggregator already exists. fail 
        return false;
      }
    }
    
    
    /**
     * Performs an immediate aggregation on a key. All machines must
     * call this simultaneously. An assertion failure is triggered if
     * the key is not found. This function uses OpenMP parallelism
     * if compiled for it.
     */
    void aggregate_now(std::string& s) {
      if (aggregators.count(key) == 0) {
        ASSERT_MSG(false, "Requested aggregator %s not found", s.c_str());
        return;
      }
      
      imap_reduce_base* mr = aggregators[s];
      mr->clear();
      // ok. now we perform reduction on local data in parallel
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        imap_reduce_base* localmr = mr->clone_empty();
        if (localmr->is_vertex_map()) {
#ifdef _OPENMP
        #pragma omp for
#endif
          for (int i = 0;i < (int)graph.num_local_vertices(); ++i) {
            local_vertex_type lvertex = graph.l_vertex(i);
            if (lvertex.owner() == rmi.procid()) {
              vertex_type vertex(lvertex);
              localmr->perform_map_vertex(*context, vertex);
            }
          }
        }
        else {
          for (int i = 0;i < (int)graph.num_local_vertices(); ++i) {
            foreach(local_edge_type& e, l_vertex(i).in_edges()) {
              edge_type edge(e);
              localmr->perform_map_edge(*context, edge);
            }
          }
        }
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
          mr->add_accumulator(localmr);
        }
        delete localmr;
      }
      
      std::vector<any> gathervec(rmi.numprocs());
      gathervec[rmi.procid()] = mr->get_accumulator();
      
      rmi.gather(gathervec, 0);
      
      if (rmi.procid() == 0) {
        for (procid_t i = 1; i < rmi.numprocs(); ++i) {
          mr->add_accumulator_any(gathervec[i]);
        }
        mr->finalize(*context);
      }
      mr->clear();
      gathervec.clear();
    }
    
    
    /**
     * Requests that the aggregator with a given key be aggregated
     * every certain number of seconds when the engine is running.
     * Note that the period is prescriptive: in practice the actual
     * period will be larger than the requested period. 
     * Seconds must be > 0;
     * 
     * All machines must call simultaneously.
     * Returns true if key is found, and false otherwise.
     */
    bool aggregate_periodic(std::string& key, float seconds) {
      rmi.barrier();
      if (seconds <= 0) return;
      if (aggregators.count(key) == 0) return false;
      else aggregate_period[key] = seconds;
      return true;
    }
    
    /**
     * Performs aggregation on all keys registered with a period.
     * May be used on engine start() to ensure all periodic 
     * aggregators are executed before engine execution.
     */
    void aggregate_all_periodic() {
      std::map<std::string, size_t>::iterator iter = aggregate_period.begin();
      while (iter != aggregate_period.end()) { 
        aggregate_now(iter->first);
        ++iter;
      }
    }
    
    
    /**
     * Must be called on engine start. Initializes the internal scheduler.
     * Must be called on all machines simultaneously.
     * ncpus is really only important for the asynchronous implementation.
     * It must be equal to the number of engine threads.
     */
    void start(size_t ncpus) {
      rmi.barrier();
      schedule.clear();
      start_time = lowres_time_seconds();
      std::map<std::string, size_t>::iterator iter = aggregate_period.begin();
      while (iter != aggregate_period.end()) {
        // schedule is a max heap. To treat it like a min heap
        // I need to insert negative keys
        schedule.push(iter->first, -(start_time + iter->second));
        ++iter;
      }
      this->ncpus = ncpus;
      ti.start();
    }
    
    
    /**
     * If asynchronous aggregation is desired, this function is
     * to be called periodically on each machine. This polls the schedule to see 
     * if there is an aggregator which needs to be activated. If there is an 
     * aggregator to be started, this function will return a non-negative value. 
     * This function is thread 
     * reentrant and each activated aggregator will only return a non-negative 
     * value call to one call to tick_asynchronous().
     * 
     * If a non-negative value is returned, the asynchronous engine
     * must ensure 
     */ 
    int tick_asynchronous() {
      
    }
    
    /**
     * If synchronous aggregation is desired, this function is
     * To be called simultaneously by one thread on each machine. 
     * This polls the schedule to see if there
     * is an aggregator which needs to be activated. If there is an aggregator 
     * to be started, this function will perform aggregation.
     */ 
    void tick_synchronous() {
      // if timer has exceeded our top key
      float curtime = lowres_time_seconds();
      // note that we do not call lowres_time_seconds everytime
      // this ensures that each key will only be run at most once.
      // each time tick_synchronous is called.
      while(-schedule.top().second < curtime) {
        std::string key = schedule.top().first;
        aggregate_now(key);
        schedule.pop();
        schedule.push(key, -(lowres_time_seconds() + aggregate_period[key]));
      }
      
    }

    /**
     * Must be called on engine stop. Clears the internal scheduler
     * And resets all incomplete states.
     */
    void stop() {
      schedule.clear();
      std::map<std::string, imap_reduce_base*>::iterator iter = aggregators.begin();
      while (iter != aggregators.end()) {
        iter->second->clear();
        ++iter;
      }
    }
  }; 


}; // end of graphlab namespace
#include <graphlab/macros_undef.hpp>

#endif
