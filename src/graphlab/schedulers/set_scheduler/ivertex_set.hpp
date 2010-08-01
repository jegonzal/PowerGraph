#ifndef IVERTEX_SET_HPP
#define IVERTEX_SET_HPP

#include <graphlab/graph/graph.hpp>
#include <graphlab/schedulers/set_scheduler/set_generic_typedefs.hpp>

namespace graphlab {

  template<typename Graph> class set_scheduler;

  const size_t INSERT_EVENT = 1;
  const size_t ERASE_EVENT = 2;
  const size_t MODIFY_VERTEX_EVENT = 4;
  const size_t MODIFY_EDGE_EVENT = 8;

  /**
     A generic base class representing an vertex set which can be manipulated by 
     the set scheduler. All operations must be atomic, (including the triggering
     of changes), and should be made as fine grained as possible. 
     i.e. individual vertex locks are preferred over complete locks over the set.
  */
  template<typename Graph>
  class ivertex_set {
  public:
    
    typedef set_scheduler<Graph> set_scheduler_type;

    ivertex_set() {};
    virtual ~ivertex_set() {};

    /**
       This function is called by the set scheduler register an instance of 
       this class with the source set, the graph, and the scheduler
    */
    virtual void init(Graph* g_,
                      ivertex_set* src,
                      set_scheduler_type* sched,
                      size_t ncpus) = 0;


    /** This event is triggered when 'v' is inserted into the 
        parent set and changes should be made to this set as necessary. This
        function should be atomic*/
    virtual void insert(ivertex_set* parent, vertex_id_t v) = 0;
  
    /** This event is triggered when 'v' is removed from the 
        parent set and changes should be made to this set as necessary. This
        function should be atomic.*/
    virtual void erase(ivertex_set* parent, vertex_id_t v) = 0;
  
    /** This event is triggered when 'v' is modified. and
        changes should be made to this set as necessary. This function
        should be atomic*/
    virtual void modify_vertex(ivertex_set* parent, vertex_id_t v) = 0;

    /** This event is triggered when 'e' is modified. and
        changes should be made to this set as necessary. This function
        should be atomic*/
    virtual void modify_edge(ivertex_set* parent, edge_id_t e) = 0;

    /** This an event is triggered when the parent set is 
        reconstructed. rebuildset is the entire contents of the parent set.
        This function does not have to be threadsafe, nor atomic. */
    virtual void rebuild(ivertex_set* parent, const ss_set_type &rebuildset) = 0;
 
    /// Returns true if this set contains vertex 'vertexid'
    virtual bool has_vertex(const vertex_id_t &vertexid) const = 0;
  
    /// Returns a bitvector containing all the vertices in this set
    virtual ss_set_type get() const = 0;
  
    /// Returns the number of elements in this set
    virtual size_t size() const = 0;

    virtual bool first(vertex_id_t &i, size_t cpuid) = 0;

    /**
       This function provides an iterator over the set.
    
       Where 'i' is the last vertex id returned by next, this should return a 
       "next" vertex ID for some definition of "next".
       That is to say, repeated calls to next should eventually reach all elements
       and next() should loop around as long as there are still elements in the 
       set. This function should be atomic.
    */
    virtual bool next(vertex_id_t &i, size_t cpuid) = 0;
  
    /// Returns the name of the object.
    virtual std::string name(void) const = 0;
  
    /// Adds an event handler
    virtual void add_event_target(ivertex_set *ev) {
      if (ev != NULL) {
        all_children.push_back(ev);
      }
    }

    virtual size_t supported_events() = 0;

    virtual size_t all_dependent_events() {
      unsigned char dependentevs = supported_events();
      for (size_t i = 0; i < all_children.size(); ++i) {
        dependentevs |= all_children[i]->all_dependent_events();
      }
      return dependentevs;
    }

    virtual void resolve_event_handlers() {
      for (size_t i = 0; i < all_children.size(); ++i) {
        ivertex_set *ev = all_children[i];
        size_t supportevs = ev->all_dependent_events();
        if (supportevs & INSERT_EVENT) insert_events.push_back(ev);
        if (supportevs & ERASE_EVENT) erase_events.push_back(ev);
        if (supportevs & MODIFY_VERTEX_EVENT) modify_vertex_events.push_back(ev);
        if (supportevs & MODIFY_EDGE_EVENT) modify_edge_events.push_back(ev);
        ev->resolve_event_handlers();
      }
    }

  
  protected:
    // the set of all connected sets
    std::vector<ivertex_set*> all_children;
    // subsets of all_children, depending on what events are supported
    std::vector<ivertex_set*> insert_events;
    std::vector<ivertex_set*> erase_events;
    std::vector<ivertex_set*> modify_vertex_events;
    std::vector<ivertex_set*> modify_edge_events;

    inline void trigger_all_insert(vertex_id_t v) {
      for(size_t i = 0;i < insert_events.size(); ++i) 
        insert_events[i]->insert(this, v);
    }

    inline void trigger_all_erase(vertex_id_t v) {
      for(size_t i = 0;i < erase_events.size(); ++i) 
        insert_events[i]->erase(this, v);
    }

    inline void trigger_all_modify_vertex(vertex_id_t v) {
      for(size_t i = 0;i < modify_vertex_events.size(); ++i) 
        insert_events[i]->modify_vertex(this, v);
    }

    inline void trigger_all_modify_edge(edge_id_t e) {
      for(size_t i = 0;i < modify_edge_events.size(); ++i) 
        insert_events[i]->modify_edge(this, e);
    }

    inline void trigger_all_rebuild(const ss_set_type &rebuildset) {
      for(size_t i = 0;i < all_children.size(); ++i) 
        all_children[i]->rebuild(this, rebuildset);
    }

  };


}
#endif
