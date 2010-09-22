#ifndef GRAPHLAB_RESTRICTED_VERTEX_SET_HPP
#define GRAPHLAB_RESTRICTED_VERTEX_SET_HPP


#include <string>

#include <boost/function.hpp>


#include <graphlab/graph/graph.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/schedulers/set_scheduler/ivertex_set.hpp>
#include <graphlab/schedulers/set_scheduler/set_generic_typedefs.hpp>
#include <graphlab/schedulers/set_scheduler/set_scheduler.hpp>
#include <graphlab/scope/iscope.hpp>
#include <logger/logger.hpp>

namespace graphlab {


  /**
     This represents a set of vertices stored as a bitmagic dense bit
     vector (hopefully to soon to be...  currently a
     vector<char>). This set is an incremental subset of the parent
     set using a user provided selection function
  */
  template<typename Graph>
  class restricted_vertex_set : 
    public ivertex_set<Graph> {
  public:
    typedef Graph graph_type;
    typedef typename graph_type::vertex_data_type vertex_data_type;
    typedef ivertex_set<Graph> base;
    typedef set_scheduler<Graph> set_scheduler_type;
    typedef iscope<Graph> iscope_type;

    using base::trigger_all_insert;
    using base::trigger_all_erase;
    using base::trigger_all_modify_vertex;
    using base::trigger_all_modify_edge;
    using base::trigger_all_rebuild;


    /**
       bool selector_function(vertexid, vertex data) This is a user
       provided function which returns true if this vertex is selected
       
       bool scoped_selector_function(vertexid, boolscope, bool
       &ret_value_reschedule) This is a user provided function which
       returns true if this vertex is selected If ret_value_reschedule
       is true, the vertex will be re-evaluated through the entire set
       hierachy as if it just completed an update function.
       ret_value_reschedule is by default false.
       
       \note ret_value_reschedule=true is not properly supported and
       may break in unexpected ways. Also, using it may bring a risk
       of an infinite set update loop. (consider the scoped_selector
       which keeps modifying itself)
    */
    typedef boost::function2<bool, vertex_id_t, const vertex_data_type&> 
    selector_function_type;
    
    typedef boost::function3<bool, vertex_id_t, iscope_type&, bool&> 
    scoped_selector_function_type;



    restricted_vertex_set(selector_function_type sel,
                          bool staticrestriction_ = false) {
      basicselector = true;
      selector1 = sel;
      selector2 = NULL;
      staticrestriction = staticrestriction_;
    }
    
    restricted_vertex_set(scoped_selector_function_type sel, 
                          bool staticrestriction_ =false) {
      basicselector = false;
      selector1 = NULL;
      selector2 = sel;
      staticrestriction = staticrestriction_;
    };




    ~restricted_vertex_set() {};

    /**
       This function is called by the set scheduler register an
       instance of this class with the source set, the graph, and the
       scheduler
    */
    void init(Graph* g_,
              base *src,
              set_scheduler_type *sched_,
              size_t ncpus_) {
      g = g_;
      sched = sched_;
      // register this set with the source
      if (src) src->add_event_target(this);
      // set the bitset to be the same size as the graph
      vset.resize(g->num_vertices());
      // zero the bitset
      vset.assign(g->num_vertices(), 0);
      // build the lockset
      locks.resize(g->num_vertices());
      numv.value = 0;
      ncpus = ncpus_;
    }


    /** v will be inserted into this set when this event is
        triggered */
    void insert(base* parent, vertex_id_t v) {
      locks[v].lock();
      if (!has_vertex(v) && test_vertex(v)) {
        vset[v] = 1;
        numv.inc();
        trigger_all_insert(v);
      }
      locks[v].unlock();
    }

    /** v will be erased from this set when this event is triggered */
    void erase(base* parent, vertex_id_t v) {
      locks[v].lock();
      if (has_vertex(v)) {
        vset[v] = 0;
        numv.dec();
        trigger_all_erase(v);
      }
      locks[v].unlock();
    }

    /** pass through if this set contains this vertex */
    void modify_vertex(base* parent, vertex_id_t v) {
      locks[v].lock();
      // test the vertex
      bool t = test_vertex(v);
      bool hasv = has_vertex(v);
      if (hasv && !t) {
        // if test fails, and we currently have the vertex, this is a
        // deletion
        vset[v] = 0;
        numv.dec();
        trigger_all_erase(v);
      } else if (!hasv && t) {
        // if test succeeds, and we currently do not have the vertex,
        // this is an insertion
        vset[v] = 1;
        numv.inc();
        trigger_all_insert(v);
      } else if (hasv && t) {
        // if test succeeds, but we already have the vertex, this is
        // a modification
        trigger_all_modify_vertex(v);
      }
      locks[v].unlock();
    }

    /** pass through */ 
    void modify_edge(base* parent, edge_id_t e) {
      vertex_id_t v = g->target(e);
      locks[v].lock();  
      if (basicselector == false) {
        bool t = test_vertex(v);
        bool hasv = has_vertex(v);
        if (hasv && !t) {
          // if test fails, and we currently have the vertex, this is
          // a deletion
          vset[v] = 0;
          numv.dec();
          trigger_all_erase(v);
        } else if (!hasv && t) {
          // if test succeeds, and we currently do not have the vertex,
          // this is an insertion
          vset[v] = 1;
          numv.inc();
          trigger_all_insert(v);
        }
        if (t) trigger_all_modify_edge(e);
      }
      else {
        if (has_vertex(v)) trigger_all_modify_edge(e);
      }
      locks[v].unlock();  
    }

    /**reconstructs this set from rebuildset */
    void rebuild(base* parent, const ss_set_type &rebuildset) {
      logger(LOG_INFO, "Rebuild on %s", name().c_str());
      ss_set_type_iterator i = begin(rebuildset);
      ss_set_type_iterator iend = end(rebuildset);
      ss_set_type newrebuildset;
      vset.assign(g->num_vertices(), 0);
      numv.value = 0;
      while (i != iend) {
        if (test_vertex(*i)) {
          ss_insert(newrebuildset, *i);
          vset[*i] = 1;
          numv.inc();
        }
        i++;
      }
      trigger_all_rebuild(newrebuildset);
    }

    /// Returns true if this set contains vertex 'vertexid'
    bool has_vertex(const vertex_id_t &vertexid) const {
      return vset[vertexid];
    }

    /// Returns a bitvector containing all the vertices in this set
    ss_set_type get() const {
      ss_set_type ret;
      for (size_t i = 0;i < vset.size(); ++i) {
        if (vset[i]) ss_insert(ret, i);
      }
      return ret;
    }

    /// Returns the number of elements in this set
    size_t size() const {
      return numv.value;
    }
    bool first(vertex_id_t &i, size_t cpuid) {
      i = cpuid;
      return next(i, cpuid);
    }

    /**
       Returns the next element in the set after 'i'
    */
    bool next(vertex_id_t &i, size_t cpuid) {
      size_t startpoint = i;
      i = (i + ncpus) % vset.size();
      for (;i < vset.size(); i+=ncpus) {
        if (vset[i]) return true;
      }
      for (i = cpuid ;i < startpoint; i+=ncpus) {
        if (vset[i]) return true;
      }
      return false;
    }


    /// Returns the name of the object.
    std::string name(void) const {
      return "Restricted Vertex Set";
    }

    size_t supported_events() {
      if (staticrestriction) {
        return INSERT_EVENT | ERASE_EVENT;
      } else {
        if (basicselector) {
          return INSERT_EVENT | ERASE_EVENT | MODIFY_VERTEX_EVENT;
        } else {
          return INSERT_EVENT | ERASE_EVENT | MODIFY_VERTEX_EVENT
            | MODIFY_EDGE_EVENT;
        }
      }
    }

  private:
    Graph* g;
    std::vector<char> vset;

    std::vector<spinlock> locks;
    atomic<size_t> numv;
    set_scheduler_type* sched;
    size_t ncpus;
  
    selector_function_type selector1;
    scoped_selector_function_type selector2;
    bool basicselector;
    bool staticrestriction;
  
    bool test_vertex(vertex_id_t v) {
      if (basicselector) {
        return selector1(v, g->vertex_data(v));
      } else {
        ASSERT_MSG(false, "Non-basic selector functionality has been removed.");
//         bool repropagatechanges = false;
//         // get the scope
//         iscope_type* scope = 
//           sched->scope_factory().get_scope(thread::thread_id(), v);
//         // run the selector
//         bool ret = selector2(v, *scope, repropagatechanges);
//         // commit the changes 
//         std::vector<edge_id_t> changes;
//         scope->commit();
//         // release the scope and inform the scheduler about the changes
//         // if repropagate changes is set
//         sched->scope_factory().release_scope(scope);
//         if (repropagatechanges)
//           sched->scoped_modifications(thread::thread_id(), v, changes);
//         return ret;
      }
    }


  };


}


#endif
