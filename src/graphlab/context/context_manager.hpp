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



#ifndef GRAPHLAB_CONTEXT_MANAGER_HPP
#define GRAPHLAB_CONTEXT_MANAGER_HPP

#include <vector>

#include <graphlab/parallel/pthread_tools.hpp>

#include <graphlab/context/icontext.hpp>
#include <graphlab/context/context.hpp>
#include <graphlab/context/consistency_model.hpp>






#include <graphlab/macros_def.hpp>
namespace graphlab {


  template<typename Engine>
  class context_manager {

  public:
    typedef Engine engine_type;
    typedef typename engine_type::graph_type graph_type;
    typedef typename engine_type::update_functor_type update_functor_type;
    
    typedef icontext<graph_type, update_functor_type> icontext_type;
    typedef context<engine_type> context_type;

    typedef typename engine_type::ischeduler_type ischeduler_type;
    
    typedef typename graph_type::vertex_id_type  vertex_id_type;
    typedef typename graph_type::edge_list_type  edge_list_type;
    typedef typename graph_type::edge_type       edge_type;

    typedef typename graph_type::vertex_data_type vertex_data_type;
 
    typedef typename context_type::cache_map_type cache_map_type;
    typedef typename context_type::cache_entry cache_entry_type;

    
  private:
    graph_type* graph_ptr;
    std::vector<context_type> contexts;
    std::vector<rwlock> locks;
    consistency_model default_consistency;

  public:

    context_manager(engine_type* engine_ptr,
                    ischeduler_type* ischeduler_ptr,
                    graph_type* graph_ptr,
                    size_t ncpus,
                    consistency_model default_consistency_range 
                    = EDGE_CONSISTENCY) :
      graph_ptr(graph_ptr), contexts(ncpus), locks(graph_ptr->num_vertices()),
      default_consistency(default_consistency_range) {
      ASSERT_NE(engine_ptr, NULL);
      ASSERT_NE(ischeduler_ptr, NULL);
      ASSERT_NE(graph_ptr, NULL);    

      if (default_consistency == DEFAULT_CONSISTENCY)
        default_consistency = EDGE_CONSISTENCY;
      // Initialize all the contexts
      for(size_t i = 0; i < contexts.size(); ++i) {
        contexts[i] = context_type(engine_ptr, graph_ptr, ischeduler_ptr, i);
      }
    } // end of context manager

    /**
     * Start is called before the threads are launch after the engine
     * is started
     */
    void start() {
      // Initialize the start time for all the contexts
      const float start_time = lowres_time_seconds();
      for(size_t i = 0; i < contexts.size(); ++i) 
        contexts[i].set_start_time(start_time);
    } // end of start

    void set_default_consistency(consistency_model 
                                 default_consistency_range) {
      default_consistency = default_consistency_range;
      if (default_consistency == DEFAULT_CONSISTENCY) 
        default_consistency = EDGE_CONSISTENCY;
    }


    // -----------------ACQUIRE CONTEXT-----------------------------
    context_type& get_context(size_t cpuid,
                              vertex_id_type v,
                              consistency_model context = 
                              DEFAULT_CONSISTENCY) {
      // Verify that the cpuid and vertex id are valid
      ASSERT_LT(cpuid, contexts.size());
      ASSERT_LT(v, locks.size());
      if (context == DEFAULT_CONSISTENCY) context = default_consistency;
      switch(context){
      case VERTEX_CONSISTENCY:
        return get_vertex_context(cpuid, v);
      case EDGE_CONSISTENCY:
        return get_edge_context(cpuid, v);
      case FULL_CONSISTENCY:
        return get_full_context(cpuid, v);
      default:
        logstream(LOG_FATAL) << "UNREACHABLE STATE!" << std::endl;
        return get_edge_context(cpuid,v);
      }
    } // end of get context
    


    void flush_cache(size_t cpuid) { }

    void acquire_writelock(const size_t cpuid, const vertex_id_type vid, 
                           const bool is_center = false) {
      locks[vid].writelock();
    }




    void acquire_readlock(const size_t cpuid, const vertex_id_type vid) {      
      locks[vid].readlock();
    }

    
    void release_lock(size_t cpuid, vertex_id_type vid) {
      locks[vid].unlock();
    }

    
    context_type& get_full_context(size_t cpuid, vertex_id_type v) {
      // grab the context
      context_type& context(contexts[cpuid]);
      
      context.init(v, FULL_CONSISTENCY);

      const edge_list_type inedges =  graph_ptr->in_edges(v);
      const edge_list_type outedges = graph_ptr->out_edges(v);

      size_t inidx = 0;
      size_t outidx = 0;

      bool curlocked = false;
      const vertex_id_type numv = 
        vertex_id_type(graph_ptr->num_vertices());
      vertex_id_type curv = context.vertex_id();
      vertex_id_type inv  = 
        (inedges.size() > 0) ? inedges[0].source() : numv;
      vertex_id_type outv  = 
        (outedges.size() > 0) ? outedges[0].target() : numv;

      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curlocked && curv < inv  && curv < outv) {
          acquire_writelock(cpuid, curv, true); // locks[curv].writelock();
          curlocked = true;
          curv = numv;
        } else if (inv < outv) {
          acquire_writelock(cpuid, inv); // locks[inv].writelock(); 
          ++inidx;
          inv = (inedges.size() > inidx) ? 
            inedges[inidx].source() : numv;
        } else if (outv < inv) {
          acquire_writelock(cpuid, outv); // locks[outv].writelock(); 
          ++outidx;
          outv = (outedges.size() > outidx) ? 
            outedges[outidx].target() : numv;
        } else if (inv == outv) {
          acquire_writelock(cpuid, inv); // locks[inv].writelock();
          ++inidx; ++outidx;
          inv = (inedges.size() > inidx) ? 
            inedges[inidx].source() : numv;
          outv = (outedges.size() > outidx) ? 
            outedges[outidx].target() : numv;
        }
      }
      // just in case we never got around to locking it
      if (!curlocked) {
        acquire_writelock(cpuid, curv, true); // locks[curv].writelock();
      }
      return context;
    } // end of get_full_context




    context_type& get_edge_context(size_t cpuid, vertex_id_type v) {
      context_type& context = contexts[cpuid];
      
      context.init(v, EDGE_CONSISTENCY);

      const edge_list_type inedges =  graph_ptr->in_edges(v);
      const edge_list_type outedges = graph_ptr->out_edges(v);

      size_t inidx = 0;
      size_t outidx = 0;

      bool curlocked = false;
      const vertex_id_type numv = vertex_id_type(graph_ptr->num_vertices());
      vertex_id_type curv = context.vertex_id();
      vertex_id_type inv  = 
        (inedges.size() > 0) ? inedges[0].source() : numv;
      vertex_id_type outv  = 
        (outedges.size() > 0) ? outedges[0].target() : numv;
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curlocked && curv < inv  && curv < outv) {
          acquire_writelock(cpuid, curv, true);
          curlocked = true;
          curv = numv;
        } else if (inv < outv) {
          acquire_readlock(cpuid, inv); 
          ++inidx;
          inv = (inedges.size() > inidx) ? 
            inedges[inidx].source() : numv;
        } else if (outv < inv) {
          acquire_readlock(cpuid, outv);
          ++outidx;
          outv = (outedges.size() > outidx) ? 
            outedges[outidx].target() : numv;
        } else if (inv == outv){
          acquire_readlock(cpuid, inv); 
          ++inidx; ++outidx;
          inv = (inedges.size() > inidx) ? 
            inedges[inidx].source() : numv;
          outv = (outedges.size() > outidx) ? 
            outedges[outidx].target() : numv;
        }
      }
      // just in case we never got around to locking it
      if (!curlocked) {
        acquire_writelock(cpuid, curv, true); // locks[curv].writelock();
      }
      return context;
    } // end of get edge context

    

    context_type& get_single_edge_context(size_t cpuid,
                                          vertex_id_type center_vid,
                                          const edge_type& edge,
                                          consistency_model consistency) {
      if (consistency == DEFAULT_CONSISTENCY) consistency = default_consistency;
      context_type& context = contexts[cpuid];
      context.init(center_vid, consistency);      
      const bool src_lt_trg = edge.source() < edge.target();
      const vertex_id_type v1 = src_lt_trg? edge.source() : edge.target();
      const vertex_id_type v2 = src_lt_trg? edge.target() : edge.source();
      ASSERT_NE(v1, v2);
      ASSERT_TRUE(v1 == center_vid || v2 == center_vid);
      switch(consistency) {
      case FULL_CONSISTENCY:
        if(v1 == center_vid) {
          acquire_readlock(cpuid, v1);
          acquire_writelock(cpuid, v2);
        } else {                   
          acquire_writelock(cpuid, v1);
          acquire_readlock(cpuid, v2);
        }
        break;
      case DEFAULT_CONSISTENCY:
      case EDGE_CONSISTENCY:
        acquire_readlock(cpuid, v1);
        acquire_readlock(cpuid, v2);
        break;
      case VERTEX_CONSISTENCY:
        if(v1 == center_vid) acquire_readlock(cpuid, v1);
        else acquire_readlock(cpuid, v2);        
        break;
      case NULL_CONSISTENCY: // NOP
        break;
      }
      return context;
    } // end of get single edge context
    


    context_type& get_vertex_context(size_t cpuid, vertex_id_type v) {
      context_type& context = contexts[cpuid];      
      context.init(v, VERTEX_CONSISTENCY);
      const vertex_id_type curv = context.vertex_id();
      acquire_writelock(cpuid, curv, true); // locks[curv].writelock();     
      return context;
    }
    
    context_type& get_null_context(size_t cpuid, vertex_id_type v) {
      context_type& context = contexts[cpuid];      
      context.init(v, NULL_CONSISTENCY);
      return context;
    }

    /**
     * Get a global context 
     */
    iglobal_context& get_global_context(size_t cpuid) { return contexts[cpuid]; }




    // -----------------RELEASE CONTEXT-----------------------------
    void release_context(size_t cpuid, context_type& context) {
      switch(context.consistency()) {
      case VERTEX_CONSISTENCY:
        release_vertex_context(cpuid, context);
        break;
      case EDGE_CONSISTENCY:
      case FULL_CONSISTENCY:
        release_full_edge_context(cpuid, context);
        break;
      case NULL_CONSISTENCY:
        release_null_context(cpuid, context);
        break;
      default:
        ASSERT_TRUE(false);
      }
    } // end of release context


    void release_single_edge_context(size_t cpuid, 
                                     context_type& context,
                                     const edge_type& edge) {
      const bool src_lt_trg = edge.source() < edge.target();
      const vertex_id_type v1 = src_lt_trg? edge.source() : edge.target();
      const vertex_id_type v2 = src_lt_trg? edge.target() : edge.source();
      ASSERT_NE(v1, v2);
      ASSERT_TRUE(v1 == context.vertex_id() || v2 == context.vertex_id());
      switch(context.consistency()) {
      case FULL_CONSISTENCY:
      case DEFAULT_CONSISTENCY:
      case EDGE_CONSISTENCY:
        release_lock(cpuid, v1);
        release_lock(cpuid, v2);
        break;
      case VERTEX_CONSISTENCY:
        release_lock(cpuid, context.vertex_id());
        break;
      case NULL_CONSISTENCY: // NOP
        break;
      }
    } // end of release single edge context


    void release_full_edge_context(size_t cpuid,
                                   context_type& context) {
      const vertex_id_type v = context.vertex_id();
      const edge_list_type inedges =  graph_ptr->in_edges(v);
      const edge_list_type outedges = graph_ptr->out_edges(v);
      size_t inidx = inedges.size() - 1;
      size_t outidx = outedges.size() - 1;

      bool curvunlocked = false;

      vertex_id_type curv = context.vertex_id();
      vertex_id_type inv  = (inedges.size() > inidx) ?
        (inedges[inidx].source()) : vertex_id_type(-1);
      vertex_id_type outv  = (outedges.size() > outidx) ?
        (outedges[outidx].target()) : vertex_id_type(-1);
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curvunlocked && (curv > inv || inv == vertex_id_type(-1)) &&
            (curv > outv || outv == vertex_id_type(-1))) {
          release_lock(cpuid, curv);
          curvunlocked = true;
        } else if ((inv+1) > (outv+1)) {
          release_lock(cpuid, inv); 
          --inidx;
          inv  = (inedges.size() > inidx) ?
            inedges[inidx].source() : vertex_id_type(-1);
        } else if ((outv+1) > (inv+1)) {
          release_lock(cpuid, outv); 
          --outidx;
          outv  = (outedges.size() > outidx) ?
            outedges[outidx].target() : vertex_id_type(-1);
        } else if (inv == outv){
          release_lock(cpuid, inv);
          --inidx; --outidx;
          inv  = (inedges.size() > inidx) ?
            inedges[inidx].source() : vertex_id_type(-1);
          outv  = (outedges.size() > outidx) ?
            outedges[outidx].target() : vertex_id_type(-1);
        }
      }

      if (!curvunlocked) {
        release_lock(cpuid, curv); 
      }
    }

    void release_vertex_context(size_t cpuid, 
                                context_type& context) {
      vertex_id_type curv = context.vertex_id();
      release_lock(cpuid, curv); 
    }

    void release_null_context(size_t cpuid, 
                              context_type& context) { }

    size_t num_vertices() const { return graph_ptr->num_vertices(); }

    
  }; // end of context_manager


}; // end of graphlab namespace
#include <graphlab/macros_undef.hpp>



#endif

