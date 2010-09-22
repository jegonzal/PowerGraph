#ifndef VERTEX_SET_HPP
#define VERTEX_SET_HPP

#include <string>

#include <graphlab/graph/graph.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/schedulers/set_scheduler/ivertex_set.hpp>
#include <graphlab/schedulers/set_scheduler/set_generic_typedefs.hpp>
#include <logger/logger.hpp>

namespace graphlab {

  template<typename Graph> class set_scheduler;
  /**
     The most basic instantiation of the ivertex_set. This represents a set of
     vertices stored as a dense bit vector. This set does no additional
     processing of its own. Attaching this set to another set will result
     an incremental clone of the source set (in bitvector form)
  */
  template<typename Graph>
  class vertex_set : 
    public ivertex_set<Graph> {
  public:

    typedef ivertex_set<Graph> ivertex_set_type;
    typedef set_scheduler<Graph> set_scheduler_type;


    using ivertex_set_type::trigger_all_insert;
    using ivertex_set_type::trigger_all_erase;
    using ivertex_set_type::trigger_all_modify_vertex;
    using ivertex_set_type::trigger_all_modify_edge;
    using ivertex_set_type::trigger_all_rebuild;


    vertex_set() : isrootset(false) {};
    ~vertex_set() {};

    /**
       This function is called by the set scheduler register an instance of
       this class with the source set, the graph, and the scheduler
    */
    void init(Graph* g_,
              ivertex_set_type* src,
              set_scheduler_type* sched,
              size_t ncpus_) {
      g = g_;
      // register this set with the source
      if (src) src->add_event_target(this);
      // set the bitset to be the same size as the graph
      vset.resize(g->num_vertices());
      vset.assign(g->num_vertices(), 0);
      // build the lockset
      locks.resize(g->num_vertices());
      numv.value = 0;
      ncpus = ncpus_;
    }

    void set_as_root_set() { isrootset = true; }

    /** v will be inserted into this set when this event is triggered */
    void insert(ivertex_set_type* parent, vertex_id_t v) {
      locks[v].lock();
      if (vset[v] == 0) {
        vset[v] = 1;
        numv.inc();
        trigger_all_insert(v);
      }
      locks[v].unlock();
    }

    /** v will be erased from this set when this event is triggered */
    void erase(ivertex_set_type* parent, vertex_id_t v) {
      locks[v].lock();
      if (vset[v] == 1) {
        vset[v] = 0;
        numv.dec();
        trigger_all_erase(v);
      }
      locks[v].unlock();
    }

    /** pass through if this set contains this vertex */
    void modify_vertex(ivertex_set_type* parent, vertex_id_t v) {
      if (isrootset) {
        trigger_all_modify_vertex(v);
      } else {
        locks[v].lock();
        if (vset[v]) trigger_all_modify_vertex(v);
        locks[v].unlock();
      }
    }

    /** pass through */
    void modify_edge(ivertex_set_type* parent, edge_id_t e) {
      if (isrootset) {
        trigger_all_modify_edge(e);
      }
      else {
        // send the modification only if target is here
        vertex_id_t v = g->target(e);
        locks[v].lock();  
        if (vset[v]) trigger_all_modify_edge(e);
        locks[v].unlock();   
      }
    }

    /**reconstructs this set from rebuildset */
    void rebuild(ivertex_set_type* parent, const ss_set_type &rebuildset) {
      logger(LOG_INFO, "Rebuild on %s", name().c_str());
      ss_set_type_iterator i = begin(rebuildset);
      ss_set_type_iterator iend = end(rebuildset);
      vset.resize(g->num_vertices());
      vset.assign(g->num_vertices(), 0);
      while (i != iend) {
        vset[*i] = 1;
        i++;
      }
      numv.value = ss_size(rebuildset);
      trigger_all_rebuild(rebuildset);
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
      return "Basic Vertex Set";
    }

    size_t supported_events() {
      if (!isrootset) {
        return INSERT_EVENT | ERASE_EVENT;
      } else {
        return 0;
      }
    }
  
  private:
    Graph *g;
    std::vector<unsigned char> vset;
    bool isrootset;
    size_t ncpus;
    std::vector<spinlock> locks;
    atomic<size_t> numv;
  };


}
#endif

