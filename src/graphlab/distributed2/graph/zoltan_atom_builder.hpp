#ifndef GRAPHLAB_ZOLTAN_ATOM_BUILDER
#define GRAPHLAB_ZOLTAN_ATOM_BUILDER


#include <map>


#include <boost/filesystem.hpp>


#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/distributed2/graph/graph_fragment.hpp>






#include <graphlab/macros_def.hpp>
namespace graphlab {



  template<typename VertexData, typename EdgeData>
  class zoltan_atom_builder {
  public:
    typedef VertexData vertex_data_type;
    typedef EdgeData edge_data_type;
    typedef fragment::description description_type;
    typedef fragment::vertex_record<VertexData, EdgeData> vertex_record_type;

    struct fragment_details {
      description_type desc;
      std::vector<vertex_id_t> vrecord_offset;
    };

  private: // data members
    size_t machine_id;
    std::string path;
    std::map<size_t, fragment_details> fragments;
    std::map<size_t, size_t> fragment2machine;
    std::vector<size_t> vrecord_offset;
    std::vector<size_t>
    


  public:
    
    zoltan_atom_builder(size_t machine_id, 
                        const std::string& path) : 
      machine_id(machine_id), path(path) {

      // Get all the local files
      std::vector<std::string> local_files; 
      namespace bfs = boost::filesystem;
      bfs::path path(pathname);
      assert(bfs::exists(path));
      for(bfs::directory_iterator iter( path ), end_iter; 
          iter != end_iter; ++iter) {
        if( !bfs::is_directory(iter->status()) ) {
          local_files.push_back(iter->path().filename());
        }
      }
      std::sort(local_files.begin(), local_files.end());

      // Fill out the fragment details vector
      foreach(const std::string& filename, local_files) {
        std::string abs_filename(path + "/" + filename);
        std::ifstream fin(abs_filename.c_str(),
                              std::ios::binary | std::ios::in);
        assert(fin.good());
        iarchive iarc(fstream);        
        description_type desc;
        iarc >> desc;
      
        fragment_details& details(fragments[desc.id]);
        details.desc = desc;
        fragment2machine[desc.id] = machine_id;

      }

      // Load all the fragment structural information into memory
      

      

    }
    




    



  }; // End of class zoltan_atom_builder


      

//         Std::vector<vertex_id_t> tmpvec;
//         // read in all the objects in the fragment file and update the detials
//         for(size_t i = 0; i < desc.num_local_verts(); ++i) {
//           assert(fin.good());
//           size_t offset(fin.tellg());
//           vertex_record_type vrec;
//           iarc >> vrec;
//           assert(vrec.vid == i + desc.begin_vertex);
//           // update the neighbors map
//           std::sort(vrec.in_neighbors.begin(), vrec.in_neighbors.end());
//           std::sort(vrec.out_neighbors.begin() vrec.out_neighbors.end());
//           std::merge(vrec.in_neighbors.begin(), vrec.in_neighbors.end(),
//                      vrec.out_neighbors.begin(), vrec.out_neighbors.end(),
//                      tmpvec.begin());
//         }

























}; // end graphlab namspace





#include <graphlab/macros_undef.hpp>
#endif

