#ifndef MULTINOMIAL_VERTEX_SET_HPP
#define MULTINOMIAL_VERTEX_SET_HPP

#include <string>

#include <graphlab/graph/graph.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/fast_multinomial.hpp>
#include <graphlab/schedulers/set_scheduler/set_generic_typedefs.hpp>
#include <graphlab/schedulers/set_scheduler/ivertex_set.hpp>
namespace graphlab {


  

  /**
     A vertex set built around the bitmagic bitvector
  */
  template<typename Graph>
  class multinomial_vertex_set : 
    public ivertex_set<Graph> {
  public:

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
       double priority_function(vertexid, vertexblob) 
       This is a user provided function which returns the priority of the vertex
    */
    typedef boost::function2<double, vertex_id_t, const vertex_data_type&> 
    priority_function_type;


    multinomial_vertex_set(priority_function_type p) {
      pfunction = p;
      vset = NULL;
    }

    multinomial_vertex_set(double (*p)(vertex_id_t, const vertex_data_type&)) {
      pfunction = p;
      vset = NULL;
    };
  
    ~multinomial_vertex_set() {
      if (vset) delete vset;
    }
  
  
    void init(Graph *g_,
              base *src,
              set_scheduler_type *sched,
              size_t ncpus_) {
      g = g_;
      // register this set with the source
      if (src) src->add_event_target(this);
      // set the bitset to be the same size as the graph
      vset = new fast_multinomial(g->num_vertices(), sched->num_cpus());
      // build the lockset
      locks.resize(g->num_vertices());
    }

    /** v will be inserted into this set when this event is triggered */
    void insert(base* parent, vertex_id_t v) {
      locks[v].lock();
      if (!has_vertex(v)) {
        vset->set(v, get_priority(v));
        trigger_all_insert(v);
      }
      locks[v].unlock();
    }

    /** v will be erased from this set when this event is triggered */
    void erase(base* parent, vertex_id_t v) {
      locks[v].lock();
      if (has_vertex(v)) {
        vset->zero(v);
        trigger_all_erase(v);
      }
      locks[v].unlock();
    }

    /** pass through if this set contains this vertex */
    void modify_vertex(base* parent, vertex_id_t v) {
      locks[v].lock();
      // test the vertex
      double pr = get_priority(v);
      bool hasv = has_vertex(v);
      if (hasv && pr == 0) {
        // if priority is 0, and we currently have the vertex, this is a deletion
        vset->zero(v);
        trigger_all_erase(v);
      } else if (!hasv && pr > 0) {
        // if priority is non zero, and we currently do not have the vertex,
        // this is an insertion
        vset->set(v, pr);
        trigger_all_insert(v);
      } else if (hasv && pr > 0) {
        // if priority is non zero, but we already have the vertex, this is
        // a modification
        vset->set(v, pr);
        trigger_all_modify_vertex(v);
      }
      locks[v].unlock();
    }

    /** pass through */
    void modify_edge(base* parent, edge_id_t e) {
      // send the modification only if source/target is here
      vertex_id_t v = g->target(e);
      locks[v].lock();  
      if (has_vertex(v)) trigger_all_modify_edge(e);
      locks[v].unlock();  
    }

    /**reconstructs this set from rebuildset */
    void rebuild(base* parent, const ss_set_type &rebuildset) {
      logger(LOG_INFO, "Rebuild on %s", name().c_str());
      ss_set_type_iterator i = begin(rebuildset);
      ss_set_type_iterator iend = end(rebuildset);
      ss_set_type newrebuildset;
      vset->clear();
      while (i != iend) {
        double pr = get_priority(*i);
        if (pr > 0) {
          vset->set(*i, get_priority(*i));
          ss_insert(newrebuildset, *i);
        }
        i++;
      }
      trigger_all_rebuild(newrebuildset);
    }


    size_t size() const {
      return vset->positive_support();
    }

    /// Returns true if this set has vertex vertexid
    bool has_vertex(const vertex_id_t &vertexid) const {
      return vset->has_support(vertexid);
    }
  
    /// Returns this set as a bitvector
    ss_set_type get() const {
      // TODO:
      ss_set_type ret;
      for(size_t i = 0;i < g->num_vertices(); ++i) {
        if (vset->has_support(i)) ss_insert(ret, i);
      }
      return ret;
    }

    std::string name(void) const{
      return "Multinomial Priority";
    }
  
    bool first(vertex_id_t &i, size_t cpuid) {
      return this->next(i, cpuid);
    }
    /// Returns the next vertex in the set (in increasing value), looping around
    /// when we reach the last vertex. Returns false if there are no elements
    /// in the set
    bool next(vertex_id_t &i, size_t cpuid) {
      size_t s = 0;
      bool ret = vset->sample(s, cpuid);
      i = vertex_id_t(s);
      return ret;
    }

    size_t supported_events() {
      return INSERT_EVENT | ERASE_EVENT | MODIFY_VERTEX_EVENT;
    }


  private:
    Graph *g;
    fast_multinomial *vset;
    priority_function_type pfunction;
    std::vector<spinlock> locks;

    double get_priority(vertex_id_t v) {
      return pfunction(v, g->vertex_data(v));
    }
  };

};

#endif
