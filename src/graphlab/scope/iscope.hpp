/** \file 
 *
 * This file describes the iscope interface as well as the the
 * scope_range_enum.
 *
 */
#ifndef GRAPHLAB_SCOPE_HPP
#define GRAPHLAB_SCOPE_HPP

#include <set>
#include <vector>
#include <cassert>

#include <graphlab/graph/graph.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  /** \brief defines the types of scope consistency guarantees provided  
   
      There are several choices for consistency mechanisms in the
      graphlab framework.  Each choice determines to what extent
      adjacent vertices can be operated on in parallel.  

      <ul> 

        <li> Null Consistency: provides no guarantees allowing update
        functions to operate on the same vertex concurrently </li>

        <li> Vertex Read Consistency: On ensures that that you can
        read from the local vertex correctly</li>

        <li> Vertex Consistency: Ensures that a scope is aquired by
        only one processor at a time </li>

        <li> Edge Consistency: Ensures that adjacent vertices are not
        updated simultaneoulsy. If the update function only modifies
        the data on the scope vertex and its adjacent edges then this
        consistency model is sufficient to guarantee sequential
        consistency </li>

        <li> Fully Consistency: This consistency models guarantees
        sequential consistency but may limit the available
        parallelism.  Effectively, this consistency model ensures that
        overlapping scopes cannot be executed simultaneously.</li>

      </ul>

      The scope_range_enum is passed to the engine through the iengine
      interface or set using the engine factory.
 
   */
  struct scope_range {
    /// \brief scope types
    enum scope_range_enum {
      NULL_CONSISTENCY = 0,    ///< no locks
      VERTEX_READ_CONSISTENCY, ///< read only from self
      READ_CONSISTENCY,        ///< read from self and adjacent structures
      VERTEX_CONSISTENCY,      ///< write to self. no lock on adjacent
      EDGE_CONSISTENCY,        ///< write to self, read from adjacent structures
      FULL_CONSISTENCY,        ///< write to self and adjacent structures
      USE_DEFAULT
    };
    
    enum lock_type_enum {
      NO_LOCK = 0,
      READ_LOCK = 1,
      WRITE_LOCK = 2
    };
  };

  inline scope_range::lock_type_enum central_vertex_lock_type(scope_range::scope_range_enum srange) {
    switch (srange) {
      case scope_range::NULL_CONSISTENCY:
        return scope_range::NO_LOCK;
      case scope_range::VERTEX_READ_CONSISTENCY:
      case scope_range::READ_CONSISTENCY:
        return scope_range::READ_LOCK;
      case scope_range::VERTEX_CONSISTENCY:
      case scope_range::EDGE_CONSISTENCY:
      case scope_range::FULL_CONSISTENCY:
        return scope_range::WRITE_LOCK;
      default:
        assert(false);
        // unreachable
        return scope_range::NO_LOCK;
    }
  }

  inline scope_range::lock_type_enum adjacent_vertex_lock_type(scope_range::scope_range_enum srange) {
    switch (srange) {
      case scope_range::NULL_CONSISTENCY:
      case scope_range::VERTEX_READ_CONSISTENCY:
      case scope_range::VERTEX_CONSISTENCY:
        return scope_range::NO_LOCK;
      case scope_range::READ_CONSISTENCY:
      case scope_range::EDGE_CONSISTENCY:
        return scope_range::READ_LOCK;
      case scope_range::FULL_CONSISTENCY:
        return scope_range::WRITE_LOCK;
      default:
        assert(false);
        // unreachable
        return scope_range::NO_LOCK;
    }
  }


  inline bool scope_is_subset_of(scope_range::scope_range_enum A,
                                 scope_range::scope_range_enum B) {
    /*
    if (A==scope_range::READ_CONSISTENCY && B == scope_range::VERTEX_CONSISTENCY) return false;
    else return (A < B);*/
    
    return (!(A==scope_range::READ_CONSISTENCY && B == scope_range::VERTEX_CONSISTENCY)) && (A < B);
  }

  /**
   * \brief represents the data associated with a vertex its adjacent
   * edges and neighbors.
   *
   * The update function is passed an instance of iscope and is uses
   * that instance to obtain information about the graph structure and
   * to read and modify the graph data.
   */
  template<typename Graph>
  class iscope {
  public:
    //! The type of graph that the iscope operates on
    typedef Graph graph_type;

    //! The vertex data type associated with the graph
    typedef typename Graph::vertex_data_type vertex_data_type;

    //! The edge data type associated with the graph
    typedef typename Graph::edge_data_type   edge_data_type;

    //! The edge data type associated with the graph
    typedef typename Graph::edge_list_type   edge_list_type;



  public:    

    /** 
     * \brief construct an iscope from a graph This is called by the
     * engine when creating an iscope to be passed into an update
     * function.
     */
    iscope(Graph* graph_ptr = NULL, vertex_id_t vertex = -1) : 
      _graph_ptr(graph_ptr), _vertex(vertex) {
    }
    
    /** iscope destructor */
    virtual ~iscope() { }
    
    /**
     * \brief commits all changes.
     * This is called by the engine after the update function returns.
     */
    virtual void commit() = 0;

    /**
     * Get the number of vertices in the graph
     */
    size_t num_vertices() const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->num_vertices();  
    }

    
    /**
     * \brief Returns the vertex id of the base vertex in this scope.
     *
     * This method is used by the update function to get the base
     * vertex of the scope.  The base vertex is the vertex that the
     * update function is being applied to.
     */
    vertex_id_t vertex() const { return _vertex; }

    /** 
     * \brief edge lookup from source target pair to edge id.
     *
     * This is used to get structural information from the graph.  If
     * the edge is not present this method will fail. 
     */
    edge_id_t edge(vertex_id_t source,
                   vertex_id_t target) const {
      assert(_graph_ptr != NULL);
      // No cheating
      // assert(source == _vertex || target == _vertex);
      return _graph_ptr->edge_id(source, target);
    }

    /**
     * \brief test whether an edge is present
     *
     * This method tests whether the edge exists.  If the edge exists
     * this method returns true.  
     */
    bool edge_exists(vertex_id_t source,
                     vertex_id_t target) const {
      assert(_graph_ptr != NULL);
      // No cheating
      // assert(source == _vertex || target == _vertex);
      return _graph_ptr->find(source, target).first;
    }
    

    /**
     * \brief Get the reverse edge.
     *
     * Get the reverse edge id.  If no such edge exists this method
     * will fail. 
     */
    edge_id_t reverse_edge(edge_id_t eid) const {      
      assert(_graph_ptr != NULL);      
      //       // No cheating
      //       assert(_graph_ptr->source(eid) == _vertex ||
      //              _graph_ptr->target(eid) == _vertex);
      return _graph_ptr->rev_edge_id(eid);
    }

    /** 
     * \brief get all in edges to the base vertex of this scope. 
     * 
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    edge_list_type in_edge_ids() const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->in_edge_ids(_vertex);
    }

    /** 
     * \brief get all in edge ids to the vertex argument
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    edge_list_type in_edge_ids(vertex_id_t v) const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->in_edge_ids(v);
    }


    /** 
     * \brief get all out edge ids to the base vertex of this scope
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    edge_list_type out_edge_ids() const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->out_edge_ids(_vertex);
    }

    /** 
     * \brief get all out ede ids to the vertex argument.
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    edge_list_type out_edge_ids(vertex_id_t v) const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->out_edge_ids(v);
    }

    //! Get the source vertex of the edge id argument
    vertex_id_t source(edge_id_t edge_id) const {
      //      assert(_graph_ptr != NULL);
      return _graph_ptr->source(edge_id);
    }

    //! get the target vertex of the edge id argument
    vertex_id_t target(edge_id_t edge_id) const {
      //       assert(_graph_ptr != NULL);
      return _graph_ptr->target(edge_id);
    }
    
    /**
     * \brief Get a mutable reference to the data associated with the
     * base vertex
     *
     * Certain optimizations may be made in future version of graphlab
     * that could benefit from immutable vertex_data requests.
     * Therefore if the vertex data does not need to be mutable use
     * the const reference version of vertex_data.
     */
    virtual vertex_data_type& vertex_data() = 0;

    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     * \deprecated use const_vertex_data
     * This should be called if the data does not need to be modified.
     */
    virtual const vertex_data_type& vertex_data() const = 0;
    
    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     *
     * This should be called if the data does not need to be modified.
     */    
    virtual const vertex_data_type& const_vertex_data() const = 0;
    
    
    /**
     * \brief Get a mutable reference to the data associated with the
     * edge.
     *
     * This should only be invoked on edges that are adjacent to the
     * base vertex.  If the edge_data is not going to be modified the
     * const version of this function should be used to permit further
     * optimization.
     */
    virtual edge_data_type& edge_data(edge_id_t eid) = 0;

    /**
     * \brief Get an immutable reference to the data associated with
     * the edge.
     * \deprecated use const_edge_data
     * This should only be invoked on edges that are adjacent to the
     * base vertex. 
     */
    virtual const edge_data_type& edge_data(edge_id_t eid) const = 0;
    
    
    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     *
     * This should be called if the data does not need to be modified.
     */    
    virtual const edge_data_type& const_edge_data(edge_id_t eid) const = 0;
    
    /**
     * \brief get a mutable reference to the data associated with a
     * neighboring vertex.
     *
     * This function should only be invoked on neighboring
     * vertices. Unfortunately, due to the Log(d) lookup required to
     * enforce the adjacency constraint we do not check at this time.
     * If the neighboring vertex data is not going to be modified the
     * const version of this function should be called to permit
     * further optimization by the graphlab engine.
     */
    virtual vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) = 0;

    /**
     * \brief get an immutable reference to the data associated with a
     * neighboring vertex.
     * \deprecated Use const_neighbor_vertex_data
     * This function should only be invoked on neighboring
     * vertices. Unfortunately, due to the Log(d) lookup required to
     * enforce the adjacency constraint we do not check at this time.
     */
    virtual const vertex_data_type& 
    neighbor_vertex_data(vertex_id_t vertex) const = 0;
        

        /**
     * \brief get an immutable reference to the data associated with a
     * neighboring vertex.
     *
     * This function should only be invoked on neighboring
     * vertices. Unfortunately, due to the Log(d) lookup required to
     * enforce the adjacency constraint we do not check at this time.
     */
    virtual const vertex_data_type& 
    const_neighbor_vertex_data(vertex_id_t vertex) const = 0;


    /**
    Experimental scope upgrade scheme. Returns true if scope upgrade is 
    successful. If this ever returns false, you are hosed. Should work
    with general_scope. Note that after scope_upgrade is called, any graph
    data within the scope may change due to a race between releasing and 
    reacquiring the upgraded scope.
    */
    virtual bool 
    experimental_scope_upgrade(scope_range::scope_range_enum newrange) { 
      return false;
    }

  protected:
    /** A pointer to the underlying graph datastructure */
    Graph* _graph_ptr;

    /** The vertex that this graph represents*/
    vertex_id_t _vertex;

  }; // end of iscope
  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif
