#ifndef ANY_INEDGE_SET_HPP
#define ANY_INEDGE_SET_HPP


#include <string>

#include <graphlab/graph/graph.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/schedulers/set_scheduler/ivertex_set.hpp>
#include <graphlab/schedulers/set_scheduler/set_generic_typedefs.hpp>
#include <graphlab/schedulers/set_scheduler/set_scheduler.hpp>
#include <graphlab/graph/iscope.hpp>
#include <graphlab/boost/function.hpp>
#include <graphlab/logger/logger.hpp>


namespace graphlab {
  /**
     bool edge_selector(srcvertex, destvertex, edgeblob) 
     This is a user provided function and should return true if this edge is selected

     bool post_update_probe_function(vertexid, vertexblob) This is a user
     provided function and should return true if the inedges of this vertex
     are to be rescheduled, and false otherwise. We guarantee that this
     function will be called immediately after an update function */
  typedef boost::function3<bool, vertex_id_t, vertex_id_t, const blob&> 
  edge_selector_function;

  typedef boost::function2<bool, vertex_id_t, const blob&> 
  post_update_probe_function;

  class set_scheduler;
  /**
     This represents a restricted set of the parent set, where any in edge satisfies
     some condition
  */
  class any_inedge_set: public ivertex_set {
  public:
    any_inedge_set(edge_selector_function sel, post_update_probe_function probe = NULL) {
      selector = sel;
      update_probe = probe;
    }
    any_inedge_set(bool (*sel)(vertex_id_t, vertex_id_t, const blob&), 
                   bool (*probe)(vertex_id_t, const blob&) = NULL) {
      selector = sel;
      update_probe = probe;
    };


    ~any_inedge_set() {};

    /**
       This function is called by the set scheduler register an instance of
       this class with the source set, the graph, and the scheduler
    */
    void init(blob_graph *g_,
              ivertex_set *src,
              set_scheduler *sched_,
              size_t ncpus_) {
      g = g_;
      sched = sched_;
      // register this set with the source
      if (src) src->add_event_target(this);
      // set the bitset to be the same size as the graph
      vset.resize(g->num_vertices());
      vset.assign(g->num_vertices(), 0);

      vset_inedgecount.resize(g->num_vertices());
      vset_inedgecount.assign(g->num_vertices(), 0);


      edgeset.resize(g->num_edges());
      // zero the bitset
      edgeset.clear();


      // build the lockset
      locks.resize(g->num_vertices());
      numv.value = 0;
      ncpus = ncpus_;
    }


    /** v will be inserted into this set when this event is triggered */
    void insert(ivertex_set* parent, vertex_id_t v) {
      locks[v].lock();
      if (!has_vertex(v) && test_vertex(v)) {
        vset[v] = 1;
        numv.inc();
        trigger_all_insert(v);
      }
      locks[v].unlock();
    }

    /** v will be erased from this set when this event is triggered */
    void erase(ivertex_set* parent, vertex_id_t v) {
      locks[v].lock();
      if (vset[v] != 0) {
        vset[v] = 0;
        vset_inedgecount[v] = 0;
        numv.dec();
        trigger_all_erase(v);
      }
      locks[v].unlock();
    }

    /** pass through if this set contains this vertex */
    void modify_vertex(ivertex_set* parent, vertex_id_t v) {
      bool reprobeall = false;
      locks[v].lock();
      // if there is no update_probe, we have to retest all the edges
      if (update_probe == NULL) {
        reprobeall = true;
      }
      else {
        // If update probe says yes, we have to retest all the edges
        reprobeall = update_probe(v, g->vertex_data(v));
      }
    
      if (reprobeall) {
        const std::vector<edge_id_t> &inedges = g->in_edge_ids(v);
        for(size_t i = 0; i < inedges.size(); ++i) {
          unsync_modify_edge(parent, inedges[i]);
        }
      }
      else {
        // if update probe says no, we delete the vertex from the set
        // unset all the inedges of v
        const std::vector<edge_id_t> &inedges = g->in_edge_ids(v);
        for(size_t i = 0; i < inedges.size(); ++i) {
          edgeset.clear_bit(inedges[i]);
        }
        vset_inedgecount[v] = 0;
        if (vset[v]) {
          numv.dec();
          vset[v] = 0;
          trigger_all_erase(v);
        }
      }
      // trigger children
      if (vset[v]) trigger_all_modify_vertex(v);
      locks[v].unlock();

    }

    /** pass through */
    void modify_edge(ivertex_set* parent, edge_id_t e) {
      vertex_id_t v = g->target(e);
      // lock the destination
      locks[v].lock();    
      unsync_modify_edge(parent, e);
      locks[v].unlock();
    }


    /** pass through */
    void unsync_modify_edge(ivertex_set* parent, edge_id_t e) {
      vertex_id_t u = g->source(e);
      vertex_id_t v = g->target(e);

      // if the parent has the destination vertex, we need to some processing..
      // since an inedge has changed

      //check the new edge value
      bool prevval = edgeset.get(e);
      bool newval = selector(u, v, g->edge_data(e));
      if (prevval && !newval) {
        // we need to retest this vertex
        --vset_inedgecount[v];
        edgeset.clear_bit(e);
        bool newinset = vset_inedgecount[v] > 0;
        if(!newinset) {
          // if it should not be in the set
          vset[v] = 0;
          numv.dec();
          trigger_all_erase(v);
        }
        else {
          trigger_all_modify_edge(e);
        }
      }
      else if (!prevval && newval) {
        // we need to insert the dest vertex if its not already in the set
        bool wasnotinset = vset_inedgecount[v] == 0;
        ++vset_inedgecount[v];
        edgeset.set_bit(e);
        if(wasnotinset) {
          vset[v] = 1;
          numv.inc();
          trigger_all_insert(v);
        }
        trigger_all_modify_edge(e);
        // else: if its already in the set, we don't need to do anything
      }
    }
    /**reconstructs this set from rebuildset */
    void rebuild(ivertex_set* parent, const ss_set_type &rebuildset) {
      logger(LOG_INFO, "Rebuild on %s", name().c_str());
      ss_set_type_iterator i = begin(rebuildset);
      ss_set_type_iterator iend = end(rebuildset);
      ss_set_type newrebuildset;

      vset.resize(g->num_vertices());
      vset.assign(g->num_vertices(), 0);

      vset_inedgecount.resize(g->num_vertices());
      vset_inedgecount.assign(g->num_vertices(), 0);

      edgeset.resize(g->num_edges());
      edgeset.clear();

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
      return "Any Inedge Set";
    }

    size_t supported_events() {
      return INSERT_EVENT | ERASE_EVENT | MODIFY_VERTEX_EVENT | MODIFY_EDGE_EVENT;
    }

  private:
    blob_graph *g;
    dense_bitset edgeset; // this edgeset caches the value of the selector
    // on inedges which are already in the set
    std::vector<unsigned char> vset;
    std::vector<uint32_t> vset_inedgecount;
    std::vector<mutex> locks;
    atomic<size_t> numv;

    set_scheduler* sched;

    edge_selector_function selector;
    post_update_probe_function update_probe;
  
    size_t ncpus;
  
    inline bool test_edge(edge_id_t e) {
      // this function keeps the edgeset up to date
      bool v = selector(g->source(e), g->target(e), g->edge_data(e));
      edgeset.set(e, v);
      return v;
    }

    bool test_vertex(vertex_id_t v) {
      // we don't have the data for this edge, so we have to check it
      // vertex tests take a while. so we acquire the scope
      // to make sure no one else gets here
      const std::vector<edge_id_t>& inedges = g->in_edge_ids(v);
      bool inset = false;
      vset_inedgecount[v] = 0;
      for (size_t i = 0;i < inedges.size(); ++i) {
        bool testresult = test_edge(inedges[i]);
        vset_inedgecount[v] += testresult;
        edgeset.set(inedges[i], testresult);
        inset |= testresult;
      }
      return inset;
    }


  };


}


#endif
