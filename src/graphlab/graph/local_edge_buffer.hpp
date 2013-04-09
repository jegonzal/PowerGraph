#ifndef GRAPHLAB_LOCAL_EDGE_BUFFER
#define GRAPHLAB_LOCAL_EDGE_BUFFER

#include <vector>
#include <graphlab/graph/graph_basic_types.hpp>

namespace graphlab {    

    template<typename VertexData, typename EdgeData>
    // Edge class for temporary storage. Will be finalized into the CSR+CSC form.
    class local_edge_buffer {
    public:
      std::vector<EdgeData> data;
      std::vector<lvid_type> source_arr;
      std::vector<lvid_type> target_arr;
    public:
      local_edge_buffer() {}
      void reserve_edge_space(size_t n) {
        data.reserve(n);
        source_arr.reserve(n);
        target_arr.reserve(n);
      }
      // \brief Add an edge to the temporary storage.
      void add_edge(lvid_type source, lvid_type target, EdgeData _data) {
        data.push_back(_data);
        source_arr.push_back(source);
        target_arr.push_back(target);
      }
      // \brief Add edges in block to the temporary storage.
      void add_block_edges(const std::vector<lvid_type>& src_arr, 
                           const std::vector<lvid_type>& dst_arr, 
                           const std::vector<EdgeData>& edata_arr) {
        data.insert(data.end(), edata_arr.begin(), edata_arr.end());
        source_arr.insert(source_arr.end(), src_arr.begin(), src_arr.end());
        target_arr.insert(target_arr.end(), dst_arr.begin(), dst_arr.end());
      }
      // \brief Remove all contents in the storage. 
      void clear() {
        std::vector<EdgeData>().swap(data);
        std::vector<lvid_type>().swap(source_arr);
        std::vector<lvid_type>().swap(target_arr);
      }
      // \brief Return the size of the storage.
      size_t size() const {
        return source_arr.size();
      }
      // \brief Return the estimated memory footprint used.
      size_t estimate_sizeof() const {
        return data.capacity()*sizeof(EdgeData) + 
          source_arr.capacity()*sizeof(lvid_type)*2 + 
          sizeof(data) + sizeof(source_arr)*2 + sizeof(local_edge_buffer);
      }
    }; // end of class local_edge_buffer.
} // end of namespace
#endif
