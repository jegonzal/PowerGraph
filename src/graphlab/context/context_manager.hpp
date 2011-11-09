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

#include <graphlab/context/icontext.hpp>
#include <graphlab/context/context.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/graph/graph.hpp>

#include <graphlab/context/idiffable.hpp>





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
    typedef typename graph_type::edge_id_type    edge_id_type;
    typedef typename graph_type::edge_list_type  edge_list_type;
    typedef typename graph_type::vertex_data_type vertex_data_type;
 
    typedef typename context_type::cache_map_type cache_map_type;
    typedef typename context_type::cache_entry cache_entry_type;
    typedef boost::is_base_of< idiffable<vertex_data_type>, vertex_data_type> is_diffable;
    
  private:
    graph_type* graph_ptr;
    std::vector<context_type> contexts;
    std::vector<rwlock> locks;
    consistency_model::model_enum default_consistency;

  public:

    context_manager(engine_type* engine_ptr,
                    ischeduler_type* ischeduler_ptr,
                    graph_type* graph_ptr,
                    size_t ncpus,
                    consistency_model::model_enum default_consistency_range 
                    = consistency_model::EDGE_CONSISTENCY) :
      graph_ptr(graph_ptr), contexts(ncpus), locks(graph_ptr->num_vertices()),
      default_consistency(default_consistency_range) {
      ASSERT_NE(engine_ptr, NULL);
      ASSERT_NE(ischeduler_ptr, NULL);
      ASSERT_NE(graph_ptr, NULL);    

      if (default_consistency == consistency_model::USE_DEFAULT)
        default_consistency = consistency_model::EDGE_CONSISTENCY;
      // Initialize all the contexts
      for(size_t i = 0; i < contexts.size(); ++i) {
        contexts[i] = context_type(engine_ptr, graph_ptr, ischeduler_ptr, i);
      }
    } // end of context manager


    void set_default_consistency(consistency_model::model_enum 
                                 default_consistency_range) {
      default_consistency = default_consistency_range;
      if (default_consistency == consistency_model::USE_DEFAULT) 
        default_consistency = consistency_model::EDGE_CONSISTENCY;
    }


    // -----------------ACQUIRE CONTEXT-----------------------------
    context_type& get_context(size_t cpuid,
                              vertex_id_type v,
                              consistency_model::model_enum context = 
                              consistency_model::USE_DEFAULT) {
      // Verify that the cpuid and vertex id are valid
      ASSERT_LT(cpuid, contexts.size());
      ASSERT_LT(v, locks.size());

      if (context == consistency_model::USE_DEFAULT) 
        context = default_consistency;      
      switch(context){
      case consistency_model::VERTEX_CONSISTENCY:
        return get_vertex_context(cpuid, v);
      case consistency_model::VERTEX_READ_CONSISTENCY:
        return get_vertex_read_context(cpuid, v);
      case consistency_model::READ_CONSISTENCY:
        return get_read_context(cpuid, v);
      case consistency_model::EDGE_CONSISTENCY:
        return get_edge_context(cpuid, v);
      case consistency_model::FULL_CONSISTENCY:
        return get_full_context(cpuid, v);
      case consistency_model::NULL_CONSISTENCY:
        return get_null_context(cpuid, v);
      default:
        logstream(LOG_FATAL) << "UNREACHABLE STATE!" << std::endl;
        return get_edge_context(cpuid,v);
      }
    } // end of get context
    


    void flush_cache(size_t cpuid) {  flush_cache(cpuid, is_diffable()); }

    //! Nothing to flush if it is not diffable
    void flush_cache(size_t cpuid, const boost::false_type&) { 
      // Check that all the caches are empty (they must be)
      foreach(const context_type& context, contexts)
        ASSERT_EQ(context.cache.size(), 0);
    } // end of flush_cache

    void flush_cache(size_t cpuid, const boost::true_type&) {
      cache_map_type& cache = contexts[cpuid].cache;
      typedef typename cache_map_type::value_type cache_pair_type;
      foreach(const cache_pair_type& pair, cache) {
        const vertex_id_type vid = pair.first;
        const cache_entry_type& entry = pair.second;
        locks[vid].writelock();
        graph_ptr->vertex_data(vid).apply_diff(entry.current, entry.old);
        locks[vid].unlock();        
      }
      // Empty the cache
      cache.clear();
    } // end of flush cache


    void acquire_writelock(const size_t cpuid, const vertex_id_type vid, 
                           const bool is_center = false) {
      acquire_writelock(cpuid, vid, is_center, is_diffable());
    }

    void acquire_writelock(const size_t cpuid, const vertex_id_type vid,
                           const bool is_center,
                           const boost::false_type&) {
      locks[vid].writelock();
    }

    void acquire_writelock(const size_t cpuid, const vertex_id_type vid,
                           const bool is_center,
                           const boost::true_type&) {
      const size_t MAX_WRITES = graph_ptr->vertex_data(vid).lag();
      // First check the cache
      context_type& context = contexts[cpuid];
      typedef typename cache_map_type::iterator iterator_type;
      iterator_type iter = context.cache.find(vid);
      const bool is_cached = iter != context.cache.end();
      if(is_cached) {
        cache_entry_type& cache_entry = iter->second;        
        if(++cache_entry.writes > MAX_WRITES || is_center) {
          //  std::cout << "Flushing" << std::endl;
          locks[vid].writelock();
          vertex_data_type& vdata = graph_ptr->vertex_data(vid);
          vdata.apply_diff(cache_entry.current, cache_entry.old);
          if(is_center) { // if it is the center vertex we evict it
            // from the cache and retain the write lock.
            context.cache.erase(vid);
            return;
          } else {
            cache_entry.current = vdata;
            locks[vid].unlock();
            cache_entry.old = cache_entry.current;
            cache_entry.writes = cache_entry.reads = 0;
          }
        } 
      } else {       
        locks[vid].readlock();
        // create a cache entry
        cache_entry_type& cache_entry = context.cache[vid];
        cache_entry.current = graph_ptr->vertex_data(vid);
        locks[vid].unlock();
        cache_entry.old = cache_entry.current;
      }
    } // end of acquire write lock



    void acquire_readlock(const size_t cpuid, const vertex_id_type vid) {      
      acquire_readlock(cpuid, vid, is_diffable());
    }

    void acquire_readlock(const size_t cpuid, const vertex_id_type vid,
                          const boost::false_type&) {
      locks[vid].readlock();
    }

    void acquire_readlock(const size_t cpuid, const vertex_id_type vid,
                          const boost::true_type&) {
      const size_t MAX_READS = graph_ptr->vertex_data(vid).lag();
      // First check the cache
      context_type& context = contexts[cpuid];
      typedef typename cache_map_type::iterator iterator_type;
      iterator_type iter = context.cache.find(vid);
      const bool is_cached = iter != context.cache.end();
      if(is_cached) {
        cache_entry_type& cache_entry = iter->second;        
        if(++cache_entry.reads > MAX_READS) {        
          locks[vid].readlock();
          const vertex_data_type& vdata = graph_ptr->vertex_data(vid);
          cache_entry.current.apply_diff(vdata, cache_entry.old);
          cache_entry.old = vdata;
          locks[vid].unlock(); 
        }
      } else {
        // Try to get the write lock
        locks[vid].readlock();
        // create a cache entry
        cache_entry_type& cache_entry = context.cache[vid];
        cache_entry.current = graph_ptr->vertex_data(vid);
        locks[vid].unlock();
        cache_entry.old = cache_entry.current;        
      }
    } // end of acquire readlock






    
    void release_lock(size_t cpuid, vertex_id_type vid) {
      release_lock(cpuid, vid, is_diffable());
    }

    void release_lock(size_t cpuid, vertex_id_type vid,
                      const boost::false_type&) {
      locks[vid].unlock();
    }

    void release_lock(size_t cpuid, vertex_id_type vid,
                      const boost::true_type&) {
      // First check the cache
      context_type& context = contexts[cpuid];
      typedef typename context_type::cache_entry cache_entry_type;
      typedef typename context_type::cache_map_type::iterator
        iterator_type;
      iterator_type iter = context.cache.find(vid);
      const bool is_cached = iter != context.cache.end();

      if(is_cached) return;
      else locks[vid].unlock();
    }







    
    context_type& get_full_context(size_t cpuid, vertex_id_type v) {
      // grab the context
      context_type& context(contexts[cpuid]);
      
      context.init(v, consistency_model::FULL_CONSISTENCY);

      const edge_list_type inedges =  graph_ptr->in_edge_ids(v);
      const edge_list_type outedges = graph_ptr->out_edge_ids(v);

      size_t inidx = 0;
      size_t outidx = 0;

      bool curlocked = false;
      const vertex_id_type numv = 
        vertex_id_type(graph_ptr->num_vertices());
      vertex_id_type curv = context.vertex_id();
      vertex_id_type inv  = 
        (inedges.size() > 0) ? graph_ptr->source(inedges[0]) : numv;
      vertex_id_type outv  = 
        (outedges.size() > 0) ? graph_ptr->target(outedges[0]) : numv;

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
            graph_ptr->source(inedges[inidx]) : numv;
        } else if (outv < inv) {
          acquire_writelock(cpuid, outv); // locks[outv].writelock(); 
          ++outidx;
          outv = (outedges.size() > outidx) ? 
            graph_ptr->target(outedges[outidx]) : numv;
        } else if (inv == outv) {
          acquire_writelock(cpuid, inv); // locks[inv].writelock();
          ++inidx; ++outidx;
          inv = (inedges.size() > inidx) ? 
            graph_ptr->source(inedges[inidx]) : numv;
          outv = (outedges.size() > outidx) ? 
            graph_ptr->target(outedges[outidx]) : numv;
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
      
      context.init(v, consistency_model::EDGE_CONSISTENCY);

      const edge_list_type inedges =  graph_ptr->in_edge_ids(v);
      const edge_list_type outedges = graph_ptr->out_edge_ids(v);

      size_t inidx = 0;
      size_t outidx = 0;

      bool curlocked = false;
      const vertex_id_type numv = vertex_id_type(graph_ptr->num_vertices());
      vertex_id_type curv = context.vertex_id();
      vertex_id_type inv  = 
        (inedges.size() > 0) ? graph_ptr->source(inedges[0]) : numv;
      vertex_id_type outv  = 
        (outedges.size() > 0) ? graph_ptr->target(outedges[0]) : numv;
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curlocked && curv < inv  && curv < outv) {
          acquire_writelock(cpuid, curv, true); // locks[curv].writelock();
          curlocked = true;
          curv = numv;
        } else if (inv < outv) {
          acquire_readlock(cpuid, inv); // locks[inv].readlock(); 
          ++inidx;
          inv = (inedges.size() > inidx) ? 
            graph_ptr->source(inedges[inidx]) : numv;
        } else if (outv < inv) {
          acquire_readlock(cpuid, outv); // locks[outv].readlock(); 
          ++outidx;
          outv = (outedges.size() > outidx) ? 
            graph_ptr->target(outedges[outidx]) : numv;
        } else if (inv == outv){
          acquire_readlock(cpuid, inv); // locks[inv].readlock();
          ++inidx; ++outidx;
          inv = (inedges.size() > inidx) ? 
            graph_ptr->source(inedges[inidx]) : numv;
          outv = (outedges.size() > outidx) ? 
            graph_ptr->target(outedges[outidx]) : numv;
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
                                          edge_id_type eid,
                                          bool writable) {
      context_type& context = contexts[cpuid];
      const consistency_model::model_enum consistency = writable? 
        consistency_model::SINGLE_EDGE_WRITE_CONSISTENCY :
        consistency_model::SINGLE_EDGE_READ_CONSISTENCY;
      context.init(center_vid, consistency);
      vertex_id_type source = graph_ptr->source(eid);
      vertex_id_type target = graph_ptr->target(eid);
      ASSERT_NE(source, target);
      ASSERT_TRUE(source == center_vid || target == center_vid);
      if(source < target) std::swap(source, target);
      if(source != center_vid && writable) {
        acquire_writelock(cpuid, source); // locks[source].writelock(); 
      } else { 
        acquire_readlock(cpuid, source); // locks[source].readlock(); 
      }
      if(target != center_vid && writable) { 
        acquire_writelock(cpuid, target); // locks[target].writelock();
      } else { 
        acquire_readlock(cpuid, target); // locks[target].readlock(); 
      }
      return context;
    } // end of get single edge context


    context_type& get_vertex_context(size_t cpuid, vertex_id_type v) {
      context_type& context = contexts[cpuid];      
      context.init(v, consistency_model::VERTEX_CONSISTENCY);
      const vertex_id_type curv = context.vertex_id();
      acquire_writelock(cpuid, curv, true); // locks[curv].writelock();     
      return context;
    }

    context_type& get_vertex_read_context(size_t cpuid, vertex_id_type v) {
      context_type& context = contexts[cpuid];      
      context.init(v, consistency_model::READ_CONSISTENCY);
      const vertex_id_type curv = context.vertex_id();
      acquire_readlock(cpuid, curv); // locks[curv].readlock();      
      return context;
    }

    context_type& get_read_context(size_t cpuid, vertex_id_type v) {
      context_type& context = contexts[cpuid];
      context.init(v, consistency_model::READ_CONSISTENCY);

      const edge_list_type inedges =  graph_ptr->in_edge_ids(v);
      const edge_list_type outedges = graph_ptr->out_edge_ids(v);

      size_t inidx = 0;
      size_t outidx = 0;

      bool curlocked = false;
      const vertex_id_type numv = (vertex_id_type)(graph_ptr->num_vertices());
      vertex_id_type curv = context.vertex_id();
      vertex_id_type inv  = 
        (inedges.size() > 0) ? graph_ptr->source(inedges[0]) : numv;
      vertex_id_type outv  = 
        (outedges.size() > 0) ? graph_ptr->target(outedges[0]) : numv;
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curlocked && curv < inv  && curv < outv) {
          acquire_readlock(cpuid, curv); // locks[curv].readlock();
          curlocked = true;
          curv = numv;
        } else if (inv < outv) {
          acquire_readlock(cpuid, inv); // locks[inv].readlock(); 
          ++inidx;
          inv = (inedges.size() > inidx) ? 
            graph_ptr->source(inedges[inidx]) : numv;
        } else if (outv < inv) {
          acquire_readlock(cpuid, outv); // locks[outv].readlock(); 
          ++outidx;
          outv= (outedges.size() > outidx) ? 
            graph_ptr->target(outedges[outidx]) : numv;
        } else if (inv == outv) {
          acquire_readlock(cpuid, inv); // locks[inv].readlock();
          ++inidx; ++outidx;
          inv = (inedges.size() > inidx) ? 
            graph_ptr->source(inedges[inidx]) : numv;
          outv = (outedges.size() > outidx) ? 
            graph_ptr->target(outedges[outidx]) : numv;
        }
      }
      // just in case we never got around to locking it
      if (!curlocked) {
        acquire_readlock(cpuid, curv); // locks[curv].readlock();
      }
      return context;
    }

    
    context_type& get_null_context(size_t cpuid, vertex_id_type v) {
      context_type& context = contexts[cpuid];      
      context.init(v, consistency_model::NULL_CONSISTENCY);
      return context;
    }



    // -----------------RELEASE CONTEXT-----------------------------
    void release_context(size_t cpuid, context_type& context) {
      switch(context.consistency()) {
      case consistency_model::VERTEX_CONSISTENCY:
      case consistency_model::VERTEX_READ_CONSISTENCY:
        release_vertex_context(cpuid, context);
        break;
      case consistency_model::EDGE_CONSISTENCY:
      case consistency_model::FULL_CONSISTENCY:
      case consistency_model::READ_CONSISTENCY:
        release_full_edge_context(cpuid, context);
        break;
      case consistency_model::NULL_CONSISTENCY:
        release_null_context(cpuid, context);
        break;
      case consistency_model::SINGLE_EDGE_READ_CONSISTENCY:
      case consistency_model::SINGLE_EDGE_WRITE_CONSISTENCY:
        logstream(LOG_FATAL) 
          << "Cannot free single edge scope with release context!"
          << std::endl;
        break;
      default:
        ASSERT_TRUE(false);
      }

      // \todo: shrink the cache as needed
    } // end of release context

    void release_single_edge_context(size_t cpuid, 
                                     context_type& context,
                                     edge_id_type eid) {
      vertex_id_type source = graph_ptr->source(eid);
      vertex_id_type target = graph_ptr->target(eid);
      ASSERT_NE(source, target);
      if(source > target) std::swap(source, target);
      release_lock(cpuid, source);
      release_lock(cpuid, target);
    } // end of release single edge context


    void release_full_edge_context(size_t cpuid,
                                   context_type& context) {
      const vertex_id_type v = context.vertex_id();
      const edge_list_type inedges =  graph_ptr->in_edge_ids(v);
      const edge_list_type outedges = graph_ptr->out_edge_ids(v);
      size_t inidx = inedges.size() - 1;
      size_t outidx = outedges.size() - 1;

      bool curvunlocked = false;

      vertex_id_type curv = context.vertex_id();
      vertex_id_type inv  = (inedges.size() > inidx) ?
        graph_ptr->source(inedges[inidx]) : vertex_id_type(-1);
      vertex_id_type outv  = (outedges.size() > outidx) ?
        graph_ptr->target(outedges[outidx]) : vertex_id_type(-1);
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curvunlocked && (curv > inv || inv == vertex_id_type(-1)) &&
            (curv > outv || outv == vertex_id_type(-1))) {
          release_lock(cpuid, curv); // locks[curv].unlock();
          curvunlocked = true;
        } else if ((inv+1) > (outv+1)) {
          release_lock(cpuid, inv); // locks[inv].unlock(); 
          --inidx;
          inv  = (inedges.size() > inidx) ?
            graph_ptr->source(inedges[inidx]) : vertex_id_type(-1);
        } else if ((outv+1) > (inv+1)) {
          release_lock(cpuid, outv); // locks[outv].unlock(); 
          --outidx;
          outv  = (outedges.size() > outidx) ?
            graph_ptr->target(outedges[outidx]) : vertex_id_type(-1);
        } else if (inv == outv){
          release_lock(cpuid, inv); // locks[inv].unlock();
          --inidx; --outidx;
          inv  = (inedges.size() > inidx) ?
            graph_ptr->source(inedges[inidx]) : vertex_id_type(-1);
          outv  = (outedges.size() > outidx) ?
            graph_ptr->target(outedges[outidx]) : vertex_id_type(-1);
        }
      }

      if (!curvunlocked) {
        release_lock(cpuid, curv); // locks[curv].unlock();
      }
    }

    void release_vertex_context(size_t cpuid, 
                                context_type& context) {
      vertex_id_type curv = context.vertex_id();
      release_lock(cpuid, curv); // locks[curv].unlock();
    }

    void release_null_context(size_t cpuid, 
                              context_type& context) { }

    size_t num_vertices() const { return graph_ptr->num_vertices(); }

    //    graph_type& graph() { return *graph_ptr; }
    
  }; // end of context_manager


}; // end of graphlab namespace
#include <graphlab/macros_undef.hpp>



#endif

