


// #include <sys/types.h>
// #include <dirent.h>


#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <graphlab.hpp>
#include <boost/filesystem.hpp>






#include <boost/filesystem.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/distributed2/graph/graph_fragment.hpp>
#include <graphlab/distributed2/graph/zoltan_atom_builder.hpp>



#include <graphlab/macros_def.hpp>

using namespace graphlab;

void list_structure_files(const std::string& pathname, 
                          std::vector<std::string>& files) {
  namespace fs = boost::filesystem;
  fs::path path(pathname);
  assert(fs::exists(path));
  for(fs::directory_iterator iter( path ), end_iter; 
      iter != end_iter; ++iter) {
    if( ! fs::is_directory(iter->status()) ) {
      std::string filename(iter->path().filename());
      std::string ending(filename.substr(filename.rfind('.')));
      if(ending == graph_fragment::structure_suffix) {
        files.push_back(iter->path().filename());
      }
    }
  }
  std::sort(files.begin(), files.end());
} // end of list_structure_file



void load_structures(const std::string& path,
                     const std::vector<std::string>& structure_fnames,
                     std::vector<graph_fragment::structure_description>& structures) {
  structures.resize(structure_fnames.size());
  for(size_t i = 0; i < structure_fnames.size(); ++i) {
    std::string abs_filename(path + "/" + structure_fnames[i]);
    std::cout << abs_filename << std::endl;
    std::ifstream fin(abs_filename.c_str(), std::ios::binary | std::ios::in);
    assert(fin.good());
    graphlab::iarchive iarc(fin);
    graph_fragment::file_description desc;
    graph_fragment::structure_description structure_desc;
    iarc >> desc >> structures[i];
    fin.close();
  }
} // end of load structures



void graphlab::build_atom_files(int arc, char** argv, 
                      const std::string& path) {
  // Load the filenames
  std::vector<std::string> structure_fnames;
  list_structure_files(path, structure_fnames);
  // Load the structures files
  std::vector<graph_fragment::structure_description> structures;
  load_structures(path, structure_fnames, structures);
  
  

} // build atom files

      

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



//   class zoltan_atom_builder {
//   public:
//     typedef graph_fragment::file_description      file_description;
//     typedef graph_fragment::structure_description structure_description;

//   private:

//     size_t machine_id;
//     std::vector<structure_description> structures;

//     std::string path;
//     std::map<size_t, fragment_details> fragments;
//     std::map<size_t, size_t> fragment2machine;
//     std::vector<size_t> vrecord_offset;
//     std::vector<size_t>
    


//   public:
    
//     zoltan_atom_builder(size_t machine_id, 
//                         const std::string& path) : 
//       machine_id(machine_id), path(path) {

//       // Get all the local files
//       std::vector<std::string> local_files; 
//       namespace bfs = boost::filesystem;
//       bfs::path path(pathname);
//       assert(bfs::exists(path));
//       for(bfs::directory_iterator iter( path ), end_iter; 
//           iter != end_iter; ++iter) {
//         if( !bfs::is_directory(iter->status()) ) {
//           local_files.push_back(iter->path().filename());
//         }
//       }
//       std::sort(local_files.begin(), local_files.end());

//       // Fill out the fragment details vector
//       foreach(const std::string& filename, local_files) {
//         std::string abs_filename(path + "/" + filename);
//         std::ifstream fin(abs_filename.c_str(),
//                               std::ios::binary | std::ios::in);
//         assert(fin.good());
//         iarchive iarc(fstream);        
//         description_type desc;
//         iarc >> desc;
      
//         fragment_details& details(fragments[desc.id]);
//         details.desc = desc;
//         fragment2machine[desc.id] = machine_id;

//       }

//       // Load all the fragment structural information into memory
      

      

//     }       



//   }; // End of class zoltan_atom_builder
 




















