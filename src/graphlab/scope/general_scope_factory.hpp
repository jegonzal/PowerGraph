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



#ifndef GRAPHLAB_GENERAL_SCOPE_FACTORY_HPP
#define GRAPHLAB_GENERAL_SCOPE_FACTORY_HPP

#include <vector>

#include <graphlab/scope/iscope.hpp>
#include <graphlab/scope/iscope_factory.hpp>
#include <graphlab/scope/general_scope.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/graph/graph.hpp>


namespace graphlab {


  template<typename Graph>
  class general_scope_factory :
    public iscope_factory<Graph> {


  public:


    typedef iscope_factory<Graph> base;
    typedef typename base::iscope_type iscope_type;
    typedef typename Graph::vertex_id_type  vertex_id_type;
    typedef typename Graph::edge_id_type    edge_id_type;
    typedef typename Graph::edge_list_type  edge_list_type;
    typedef general_scope<Graph> general_scope_type;

  private:
    Graph& graph;
    std::vector<general_scope_type*> scopes;
    std::vector<rwlock> locks;
    scope_range::scope_range_enum default_scope;

  public:

    general_scope_factory(Graph& graph,
                          size_t ncpus,
                          scope_range::scope_range_enum default_scope_range 
                          = scope_range::NULL_CONSISTENCY) :
      base(graph,ncpus), graph(graph),
      default_scope(default_scope_range) {
      if (default_scope == scope_range::USE_DEFAULT)
        default_scope = scope_range::VERTEX_CONSISTENCY;
      locks.resize(graph.num_vertices());
      
      // preallocate the scopes
      scopes.resize(2 * ncpus);
      // create them to be some arbitrary scope type
      for (size_t i = 0; i < scopes.size(); ++i) {
        scopes[i] = new general_scope_type(&graph, 0, this, 
                                           scope_range::FULL_CONSISTENCY);
      }
    }

    void set_default_scope(scope_range::scope_range_enum default_scope_range) {
      default_scope = default_scope_range;
      if (default_scope == scope_range::USE_DEFAULT) 
        default_scope = scope_range::VERTEX_CONSISTENCY;
    }

    ~general_scope_factory() { 
      for (size_t i = 0;i < scopes.size(); ++i) {
        delete scopes[i];
      }
    }

    // -----------------ACQUIRE SCOPE-----------------------------
    iscope_type* get_scope(size_t cpuid,
                           vertex_id_type v,
                           scope_range::scope_range_enum scope = scope_range::USE_DEFAULT) {
      if (scope == scope_range::USE_DEFAULT) scope = default_scope;
      
      switch(scope){
      case scope_range::VERTEX_CONSISTENCY:
        return get_vertex_scope(cpuid, v);
      case scope_range::VERTEX_READ_CONSISTENCY:
        return get_vertex_read_scope(cpuid, v);
      case scope_range::READ_CONSISTENCY:
        return get_read_scope(cpuid, v);
      case scope_range::EDGE_CONSISTENCY:
        return get_edge_scope(cpuid, v);
      case scope_range::FULL_CONSISTENCY:
        return get_full_scope(cpuid, v);
      case scope_range::NULL_CONSISTENCY:
        return get_null_scope(cpuid, v);
      default:
        ASSERT_TRUE(false);
        return NULL;
      }
    }
    
    
    iscope_type* get_full_scope(size_t cpuid, vertex_id_type v) {
      // grab the scope
      general_scope_type* scope = scopes[cpuid];
      
      scope->init(&graph, v);
      scope->stype = scope_range::FULL_CONSISTENCY;

      const edge_list_type inedges =  graph.in_edge_ids(v);
      const edge_list_type outedges = graph.out_edge_ids(v);

      size_t inidx = 0;
      size_t outidx = 0;

      bool curlocked = false;
      vertex_id_type numv = (vertex_id_type)(graph.num_vertices());
      vertex_id_type curv = scope->vertex();
      vertex_id_type inv  = (inedges.size() > 0) ? graph.source(inedges[0]) : numv;
      vertex_id_type outv  = (outedges.size() > 0) ? graph.target(outedges[0]) : numv;
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curlocked && curv < inv  && curv < outv) {
          locks[curv].writelock();
          curlocked = true;
          curv = numv;
        } else if (inv < outv) {
          locks[inv].writelock(); ++inidx;
          inv = (inedges.size() > inidx) ? graph.source(inedges[inidx]) : numv;
        } else if (outv < inv) {
          locks[outv].writelock(); ++outidx;
          outv= (outedges.size() > outidx) ? graph.target(outedges[outidx]) : numv;
        } else if (inv == outv){
          locks[inv].writelock();
          ++inidx; ++outidx;
          inv = (inedges.size() > inidx) ? graph.source(inedges[inidx]) : numv;
          outv= (outedges.size() > outidx) ? graph.target(outedges[outidx]) : numv;
        }
      }
      // just in case we never got around to locking it
      if (!curlocked) {
        locks[curv].writelock();
      }
      return scope;
    }


    iscope_type* get_edge_scope(size_t cpuid, vertex_id_type v) {
      general_scope_type* scope = scopes[cpuid];
      
      scope->init(&graph, v);
      scope->stype = scope_range::EDGE_CONSISTENCY;

      const edge_list_type inedges =  graph.in_edge_ids(v);
      const edge_list_type outedges = graph.out_edge_ids(v);

      size_t inidx = 0;
      size_t outidx = 0;


      bool curlocked = false;
      vertex_id_type numv = (vertex_id_type)(graph.num_vertices());
      vertex_id_type curv = scope->vertex();
      vertex_id_type inv  = (inedges.size() > 0) ? graph.source(inedges[0]) : numv;
      vertex_id_type outv  = (outedges.size() > 0) ? graph.target(outedges[0]) : numv;
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curlocked && curv < inv  && curv < outv) {
          locks[curv].writelock();
          curlocked = true;
          curv = numv;
        } else if (inv < outv) {
          locks[inv].readlock(); ++inidx;
          inv = (inedges.size() > inidx) ? graph.source(inedges[inidx]) : numv;
        } else if (outv < inv) {
          locks[outv].readlock(); ++outidx;
          outv= (outedges.size() > outidx) ? graph.target(outedges[outidx]) : numv;
        } else if (inv == outv){
          locks[inv].readlock();
          ++inidx; ++outidx;
          inv = (inedges.size() > inidx) ? graph.source(inedges[inidx]) : numv;
          outv= (outedges.size() > outidx) ? graph.target(outedges[outidx]) : numv;
        }
      }
      // just in case we never got around to locking it
      if (!curlocked) {
        locks[curv].writelock();
      }

      return scope;
    }

    iscope_type* get_vertex_scope(size_t cpuid, vertex_id_type v) {
      general_scope_type* scope = scopes[cpuid];
      
      scope->init(&graph, v);
      scope->stype = scope_range::VERTEX_CONSISTENCY;

      vertex_id_type curv = scope->vertex();
      locks[curv].writelock();

      return scope;
    }


    void acquire_range_lock(size_t start, size_t end) {
        for (size_t i = start; i <= end; ++i) {
          locks[i].readlock();
        }
    }

    iscope_type* get_vertex_read_scope(size_t cpuid, vertex_id_type v) {
      general_scope_type* scope = scopes[cpuid];
      
      scope->init(&graph, v);
      scope->stype = scope_range::READ_CONSISTENCY;

      vertex_id_type curv = scope->vertex();
      locks[curv].readlock();
      return scope;
    }

    iscope_type* get_read_scope(size_t cpuid, vertex_id_type v) {
      general_scope_type* scope = scopes[cpuid];
      
      scope->init(&graph, v);
      scope->stype = scope_range::READ_CONSISTENCY;

      const edge_list_type inedges =  graph.in_edge_ids(v);
      const edge_list_type outedges = graph.out_edge_ids(v);

      size_t inidx = 0;
      size_t outidx = 0;


      bool curlocked = false;
      vertex_id_type numv = (vertex_id_type)(graph.num_vertices());
      vertex_id_type curv = scope->vertex();
      vertex_id_type inv  = (inedges.size() > 0) ? graph.source(inedges[0]) : numv;
      vertex_id_type outv  = (outedges.size() > 0) ? graph.target(outedges[0]) : numv;
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curlocked && curv < inv  && curv < outv) {
          locks[curv].readlock();
          curlocked = true;
          curv = numv;
        } else if (inv < outv) {
          locks[inv].readlock(); ++inidx;
          inv = (inedges.size() > inidx) ? graph.source(inedges[inidx]) : numv;
        } else if (outv < inv) {
          locks[outv].readlock(); ++outidx;
          outv= (outedges.size() > outidx) ? graph.target(outedges[outidx]) : numv;
        } else if (inv == outv){
          locks[inv].readlock();
          ++inidx; ++outidx;
          inv = (inedges.size() > inidx) ? graph.source(inedges[inidx]) : numv;
          outv= (outedges.size() > outidx) ? graph.target(outedges[outidx]) : numv;
        }
      }
      // just in case we never got around to locking it
      if (!curlocked) {
        locks[curv].readlock();
      }
      return scope;
    }

    
    iscope_type* get_null_scope(size_t cpuid, vertex_id_type v) {
      general_scope_type* scope = scopes[cpuid];
      
      scope->init(&graph, v);
      scope->stype = scope_range::NULL_CONSISTENCY;

      return scope;
    }

    // -----------------RELEASE SCOPE-----------------------------

    void release_scope(iscope_type* scopei) {
      general_scope_type* scope = (general_scope_type*)scopei;
      switch(scope->scope_type()){
      case scope_range::VERTEX_CONSISTENCY:
      case scope_range::VERTEX_READ_CONSISTENCY:
        release_vertex_scope(scope);
        break;
      case scope_range::EDGE_CONSISTENCY:
      case scope_range::FULL_CONSISTENCY:
      case scope_range::READ_CONSISTENCY:
        release_full_edge_scope(scope);
        break;
      case scope_range::NULL_CONSISTENCY:
        release_null_scope(scope);
        break;
      default:
        ASSERT_TRUE(false);
      }
    }

    void release_full_edge_scope(general_scope_type* scope) {
      vertex_id_type v = scope->vertex();
      const edge_list_type inedges =  graph.in_edge_ids(v);
      const edge_list_type outedges = graph.out_edge_ids(v);
      size_t inidx = inedges.size() - 1;
      size_t outidx = outedges.size() - 1;

      bool curvunlocked = false;

      vertex_id_type curv = scope->vertex();
      vertex_id_type inv  = (inedges.size() > inidx) ?
        graph.source(inedges[inidx]) : vertex_id_type(-1);
      vertex_id_type outv  = (outedges.size() > outidx) ?
        graph.target(outedges[outidx]) : vertex_id_type(-1);
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curvunlocked && (curv > inv || inv == vertex_id_type(-1)) &&
            (curv > outv || outv == vertex_id_type(-1))) {
          locks[curv].unlock();
          curvunlocked = true;
        } else if ((inv+1) > (outv+1)) {
          locks[inv].unlock(); --inidx;
          inv  = (inedges.size() > inidx) ?
            graph.source(inedges[inidx]) : vertex_id_type(-1);
        } else if ((outv+1) > (inv+1)) {
          locks[outv].unlock(); --outidx;
          outv  = (outedges.size() > outidx) ?
            graph.target(outedges[outidx]) : vertex_id_type(-1);
        } else if (inv == outv){
          locks[inv].unlock();
          --inidx; --outidx;
          inv  = (inedges.size() > inidx) ?
            graph.source(inedges[inidx]) : vertex_id_type(-1);
          outv  = (outedges.size() > outidx) ?
            graph.target(outedges[outidx]) : vertex_id_type(-1);
        }
      }

      if (!curvunlocked) {
        locks[curv].unlock();
      }
    }

    void release_vertex_scope(general_scope_type* scope) {
      vertex_id_type curv = scope->vertex();
      locks[curv].unlock();
    }

    void release_range_lock(size_t start, size_t end) {
      for (size_t i = start; i <= end; ++i)
        locks[i].unlock();
    }

    void release_null_scope(general_scope_type* scope) {
    }

    size_t num_vertices() const { return graph.num_vertices(); }

    Graph& get_graph() { return graph; }
    
  };

}

#endif

