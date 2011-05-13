/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_SYNCHRONOUS_SCOPE_HPP
#define GRAPHLAB_SYNCHRONOUS_SCOPE_HPP

#include <boost/bind.hpp>


#include <graphlab/scope/iscope.hpp>


namespace graphlab {
  
  /**
   * This defines a scope type which is meant for "synchronous" type of 
   * algorithms. This type of scope should only be used by synchronous_engine
   */
  template<typename Graph>
  class synchronous_scope : 
    public iscope<Graph> {
  public:
    typedef iscope<Graph> base;
    typedef typename Graph::vertex_data_type vertex_data_type;
    typedef typename Graph::edge_data_type edge_data_type;

    using base::_vertex;
    using base::_graph_ptr;

  public:
    synchronous_scope() : base(NULL,0) { }

    // TODO: Why 3 graph arguments?  
    synchronous_scope(Graph* srcgraph, Graph* destgraph,
                      Graph* vertexdatagraph, 
                      vertex_id_t vertex) : 
      base(srcgraph, vertex)  {
      _srcgraph = srcgraph;
      _destgraph = destgraph;
      _vertexdatagraph = vertexdatagraph;
    }

    
    ~synchronous_scope() { }
    
    void commit() {};
    
    void init(Graph* srcgraph, 
              Graph* destgraph,
              Graph* vertexdatagraph, vertex_id_t vertex) {
      base::_graph_ptr = srcgraph;
      base::_vertex = vertex;
      _srcgraph = srcgraph;
      _destgraph = destgraph;
      _vertexdatagraph = vertexdatagraph;
    }
    
    /// Returns the data on the base vertex
    vertex_data_type& vertex_data()  {
      return (_vertexdatagraph->vertex_data( _vertex ));
    }

    const vertex_data_type& vertex_data() const  {
      return const_vertex_data();
    }

    const vertex_data_type& const_vertex_data() const  {
      return (_vertexdatagraph->vertex_data( _vertex ));
    }

    /// Direct calls to access edge data
    const edge_data_type& edge_data(edge_id_t eid) const { 
      return const_edge_data(eid);
    }

    /// Direct calls to access edge data
    const edge_data_type& const_edge_data(edge_id_t eid) const { 
      // TODO: make sure edge is associated with this vertex
      if (_srcgraph->target(eid) == _vertex ) {
        // in edge
        return (_srcgraph->edge_data(eid));
      } else {
        // outedge
        return (_destgraph->edge_data(eid));
      }
    }

    edge_data_type& edge_data(edge_id_t eid) {
      // TODO: make sure edge is associated with this vertex
      if (_srcgraph->target(eid) == _vertex) {
        // in edge
        return (_srcgraph->edge_data(eid));
      } else {
        // outedge
        return (_destgraph->edge_data(eid));
      }
    }
    
    const vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) const {
      return const_neighbor_vertex_data(vertex);
    }

    const vertex_data_type& const_neighbor_vertex_data(vertex_id_t vertex) const {
      return _srcgraph->vertex_data(vertex);
    }

    vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) {
      // this totally does not make sense for synchronous execution
      assert(false);
      return _srcgraph->vertex_data(vertex);
    }


      
  private:
    Graph* _srcgraph;
    Graph* _destgraph;
    Graph* _vertexdatagraph;
  }; // end of synchronous_scope
  
} // end of graphlab namespace


#endif

