#ifndef GRAPHLAB_VERTEX_CHANNEL_HPP
#define GRAPHLAB_VERTEX_CHANNEL_HPP
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/parallel/mutex.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class vertex_channel {
    typedef VertexData vertex_data_type;
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    typedef typename graph_type::mirror_type mirror_type;

   public:
    vertex_id_type vid;
    mirror_type mirrors;
    vertex_data_type vdata;
    size_t num_in_edges, num_out_edges;
    procid_t owner;
    bool data_is_set;
    mutex mtx;

    vertex_channel() : 
        vid(-1), vdata(vertex_data_type()), num_in_edges(0), num_out_edges(0), 
        owner(-1), data_is_set(false) { }

    vertex_channel& operator+=(const vertex_channel& other) {
      ASSERT_EQ(vid, other.vid);
      num_in_edges += other.num_in_edges;
      num_out_edges += other.num_out_edges;
      // mirrors += other.mirrors
      typename mirror_type::const_iterator it = other.mirrors.begin();
      while (it != other.mirrors.end()) {
        mirrors.set_bit(*it);
        ++it;
      }
      return *this;
    }

    void load(iarchive& arc) { 
      arc >> vid >> num_in_edges >> num_out_edges >> owner >> mirrors >> data_is_set;
      if (data_is_set) {
       arc >> vdata;
      }
    }

    void save(oarchive& arc) const {
      arc << vid << num_in_edges << num_out_edges << owner << mirrors << data_is_set;
      if (data_is_set) {
        arc << vdata;
      }
    }
  }; // end of vertex_channel
} // end of namespace
#endif
