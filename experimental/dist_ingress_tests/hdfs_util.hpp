#include <iostream>
#include <fstream>
#include <string>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


#include <graphlab/util/hdfs.hpp>


namespace graphlab{
  typedef double vertex_data;
  typedef double edge_data;
  typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

//   /**
//   A boost source device which can attach to hdfs file
//   */
//   struct hdfs_file_source {
//     hdfs_file_source(hdfsFS* fs, hdfsFile* f, size_t maxlen = size_t(-1)) :
//        f(f), maxlen(maxlen)  {}
    
//     hdfsFS* fs;
//     hdfsFile* f;
//     size_t maxlen;
//     typedef char      char_type;
//     struct category : public boost::iostreams::source_tag { };

//    /** to satisfy the optimally buffered tag. Since this is an
//         in-memory buffer. the optimal buffer size (for any wrapping 
//         stream) is 0. */
//     inline std::streamsize optimal_buffer_size() const { return 0; }

//     inline std::streamsize read(char* s, std::streamsize n) {
//       if ((size_t)(n) > maxlen) n = (std::streamsize)maxlen;
//       maxlen -= (size_t)n;
//       if (n == 0) return -1;
//       else return hdfsRead(*fs, *f, s, n);
//     }
//   };

 
//   bool hdfs_load_adj_structure (
//       hdfsFS* fs, hdfsFile* f, graph_type& graph, tOffset fsize, bool gzip=false) {
//       logstream(LOG_INFO) << "Loading adjacency file" << std::endl;
//         namespace bios = boost::iostreams;
//         bios::stream<hdfs_file_source> hdfs_stream(fs, f, fsize);
//         bios::filtering_stream<bios::input> fin;  
//         // Using gzip filter
//         if (gzip)
//           fin.push(bios::gzip_decompressor());

//         fin.push(hdfs_stream);
//         assert(fin.good());

//       logstream(LOG_INFO) << "file open successful" << std::endl;

//       size_t self_edges = 0;     
//       size_t ctr = 0;
//       // Loop over the contents
//       while(fin.good()) {
//         // Load a vertex
//         size_t source = 0, nneighbors = 0;
//         try { fin >> source >> nneighbors; } catch ( ... ) { 
//           logstream(LOG_WARNING) 
//             << "Error reading source." << std::endl; return false;
//         }
//         if(!fin.good()) break;
//         // Resize the graph if needed
//         if(source >= graph.num_vertices()) graph.resize(source + 1);
//         // add neighbors
//         for(size_t i = 0; i < nneighbors; ++i) {
//           size_t target = 0;
//           try { fin >> target; } catch ( ... ) {
//             logstream(LOG_WARNING) 
//               << "Error reading neighbor" << std::endl; 
//             return false;
//           }
//           if(!fin.good()) {
//             logstream(LOG_WARNING) 
//               << "Error reading neighbor" << std::endl; 
//             return false;
//           }
//           // Resize the graph if needed
//           if(target >= graph.num_vertices()) graph.resize(target + 1);
    
//           if(source != target) graph.add_edge(source, target);
//           else if(self_edges++ == 0) 
//             logstream(LOG_WARNING) 
//               << "Self edge encountered but not supported!" << std::endl
//               << "\t Further warnings will be surpressed." << std::endl;
//         } // end of loop over neighbors
//         if (++ctr % 1000000 == 0) 
//           logstream(LOG_INFO) 
//             << "Added edata for " << ctr << " vertices: " 
//             << source << std::endl; 
//       } // end of loop over file
//         // fin.close();
//       return true;
//   }
};
