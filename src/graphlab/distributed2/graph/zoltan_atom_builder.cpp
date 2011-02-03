


// #include <sys/types.h>
// #include <dirent.h>


#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <graphlab.hpp>
#include <boost/filesystem.hpp>

#include <mpi.h>
#include <zoltan_cpp.h>




#include <boost/filesystem.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/distributed2/graph/graph_fragment.hpp>
#include <graphlab/distributed2/graph/zoltan_atom_builder.hpp>
#include <graphlab/util/charstream.hpp>


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


void distribute_part2proc_map(std::map<vertex_id_t, vertex_id_t>& part2proc) {
  // Get the mpi rank and size
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Serialize the local map
  graphlab::charstream cstrm(128);
  graphlab::oarchive arc(cstrm);
  arc << part2proc;

  // All gather the size of the maps that must be received
  int buffer_size = cstrm->size();
  std::vector<int> sizes(mpi_size);
  // Compute the sizes
  MPI_Allgather(&buffer_size,  // Send buffer
                1,             // send count
                MPI_INT,       // send type
                &sizes[0],     // recvbuffer
                1,             // recvcount
                MPI_UNSIGNED,  // recvtype
                MPI_COMM_WORLD);  


  // Construct offsets
  std::vector<int> offsets = sizes;
  int sum = 0, tmp = 0;
  for(size_t i = 0; i < offsets.size(); ++i) {
    tmp = offsets[i]; offsets[i] = sum; sum += tmp; 
  }

  if(mpi_rank == 0) {
    for(size_t i = 0; i < sizes.size(); ++i) 
      std::cout << sizes[i] << '\t';
    std::cout << std::endl;
    for(size_t i = 0; i < offsets.size(); ++i) 
      std::cout << offsets[i] << '\t';
    std::cout << std::endl;
    std::cout << "Sum: " << sum << std::endl;
  }
    
  // recv all the maps
//   std::vector<char> recv_buffer(sum);
//   MPI_Allgatherv(cstrm->c_str(),  // sned buffer
//                  buffer_size,                // how much to send
//                  MPI_CHAR,                   // send type
//                  &recv_buffer[0],            // recv buffer
//                  &sizes[0],                  // amount to recv for each process
//                  &offsets[0],                // where to place data
//                  MPI_CHAR,
//                  MPI_COMM_WORLD);
                 



  



 
  
  
  
  
}


void graphlab::build_atom_files(int argc, char** argv, 
                                const std::string& path) {
  // Get the mpi rank and size
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


  // Load the filenames
  std::vector<std::string> structure_fnames;
  list_structure_files(path, structure_fnames);
  // Load the structures files
  std::vector<graph_fragment::structure_description> structures;
  load_structures(path, structure_fnames, structures);


  // Initialize the map
  std::map<vertex_id_t, vertex_id_t> part2proc;
  for(size_t i = 0; i < structures.size(); ++i) {
    part2proc[structures[i].desc.id] = mpi_rank;
  } // end of for
  distribute_part2proc_map(part2proc);
  





  // // if(mpi_rank == 0) {
  // //   std::cout << "Size: " << mpi_size << std::endl;
  // // }

  // // Initialize Zoltan
  // float zoltan_version;
  // int error = Zoltan_Initialize(argc, argv, &zoltan_version);
  // if(error != ZOLTAN_OK) {
  //   std::cout << "Unable to launch Zoltan!" << std::endl;
  //   exit(EXIT_FAILURE);
  // }
  // if(mpi_rank == 0)
  //   std::cout << "Zoltan version: " << zoltan_version << std::endl;

  
  // // Create the Zoltan object
  // Zoltan* zolt_ptr = new Zoltan(MPI::COMM_WORLD);
  // assert(zolt_ptr != NULL);
  // // Zoltan& zolt(*zolt_ptr);


  


  // // Destroy the zoltan object
  // delete zolt_ptr;

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
 




// class local_graph {

//   int rank, ncpus;
//   int nverts, nedges;
//   int block_size, block_remainder;
//   int begin_vertex, end_vertex;
//   std::vector< std::vector<int> > edges;
//   std::vector< std::string > vdata;

// public:


//   // Load the graph from file
//   local_graph(const int rank, 
//               const int ncpus, 
//               const char* filename) :
//     rank(rank), ncpus(ncpus) {

//     // Open the file
//     assert(filename != NULL);
//     std::ifstream fin(filename);
//     assert(fin.good());
    
//     // read in the vertices
//     fin >> nverts >> nedges;
//     assert(fin.good());
//     if(rank == 0) {
//       std::cout << "Vertices: " <<  nverts << std::endl
//                 << "Edges:    " <<  nedges
//                 << std::endl;
//     }
    
//     // Compute the vertex range
//     block_size = nverts / ncpus;
//     assert(block_size >= 1);
//     block_remainder = nverts % ncpus;

//     begin_vertex = block_size * rank + std::min(rank, block_remainder);
//     end_vertex = block_size * (rank + 1) + std::min(rank+1, block_remainder);

//     assert(begin_vertex < end_vertex);
//     std::cout << "Rank(" << rank << ") : " 
//               << begin_vertex << " to " << end_vertex 
//               << " with " << num_local_objects() << " local objects."
//               << std::endl;

//     // Resize the edges vector
//     edges.resize(num_local_objects());

//     // if(rank == 0) {
//     //   std::cout << "Blocksize: " << block_size << std::endl;
//     //   for(int i = 0; i < nverts; ++i) 
//     //     std::cout << owning_proc(i) << " owns " << i << std::endl;      
//     // }


//     // Read in vertices in range
//     int source = 0, target = 0;
//     while(fin.good()) {
//       fin >> source >> target;
//       if(fin.good()) {
//         assert( source >= 0 && source < nverts);
//         assert( target >= 0 && target < nverts);
//         assert( !is_local(source) || owning_proc(source) == rank );
//         assert( !is_local(target) || owning_proc(target) == rank );
//         if(is_local(source))  {
//           std::stringstream strm;
//           strm << "Source: " << source;
//           vdata[local_id(source)] = strm.str();
//           edges[local_id(source)].push_back(target);
//         }
//         if(is_local(target)) {
//           edges[local_id(target)].push_back(source);       
//         }
//       }
//     }
//     // File read completed and graph should be loaded
    
//   } // end of constructor


  
//   bool is_local(int vid) const {
//     return vid >= begin_vertex && vid < end_vertex;
//   }

//   int owning_proc(int vid) const {   
//     if(vid < (block_size+1) * block_remainder) 
//       return vid / (block_size+1);
//     else 
//       return (vid - ((block_size+1) * block_remainder)) / block_size 
//         + block_remainder;
//   }

//   int num_local_objects() const {
//     return end_vertex - begin_vertex;
//   }

//   int local_id(int gid) const {
//     int id = gid - begin_vertex;
//     assert(id >= 0);
//     assert(id < num_local_objects());
//     return id;
//   }

//   int global_id(int lid) const {
//     int id = lid + begin_vertex;
//     assert(id >= 0);
//     assert(id < nverts);
//     return id;
//   }


//   int num_edges(int lid) const {
//     assert(lid < int(edges.size()));
//     return edges[lid].size();
//   }


//   static int zoltan_num_obj_fun(void* data, int* ierr) {
//     local_graph& graph(*reinterpret_cast<local_graph*>(data));
//     *ierr = ZOLTAN_OK;
//     return graph.num_local_objects();    
//   }


//   static void zoltan_obj_list_fun(void* data, 
//                                  int num_gid_entries,
//                                  int num_lid_entries,
//                                  ZOLTAN_ID_PTR global_ids,
//                                  ZOLTAN_ID_PTR local_ids,
//                                  int wgt_dim,
//                                  float* obj_wgts,
//                                  int* ierr) {
//     local_graph& graph(*reinterpret_cast<local_graph*>(data));
//     *ierr = ZOLTAN_OK;

//     assert(num_gid_entries == 1);
//     assert(num_lid_entries == 1);

//     // set the local and global ids;
//     for(int i = 0; i < graph.num_local_objects(); ++i) {
//       local_ids[i] = i;
//       global_ids[i] = graph.global_id(i);
//     }
//   }


//   static void zoltan_num_edges_multi_fun(void* data,
//                                          int num_gid_entries,
//                                          int num_lid_entries,
//                                          int num_obj,
//                                          ZOLTAN_ID_PTR global_ids,
//                                          ZOLTAN_ID_PTR local_ids, 
//                                          int* num_edges, 
//                                          int* ierr) {
//     local_graph& graph(*reinterpret_cast<local_graph*>(data));
//     *ierr = ZOLTAN_OK;
    
//     assert(num_gid_entries == 1);
//     assert(num_lid_entries == 1);

//     assert(num_obj == graph.num_local_objects());

//     // set the num_edges field for all local edges
//     for(int i = 0; i < graph.num_local_objects(); ++i) {
//       num_edges[i] = graph.num_edges(local_ids[i]);
//     }
//   } 


//   static void zoltan_edge_list_multi_fun(void* data,
//                                          int num_gid_entries,
//                                          int num_lid_entries,
//                                          int num_obj,
//                                          ZOLTAN_ID_PTR global_ids,
//                                          ZOLTAN_ID_PTR local_ids, 
//                                          int* num_edges, 
//                                          ZOLTAN_ID_PTR nbor_global_id,
//                                          int* nbor_procs,
//                                          int wgt_dim,
//                                          float* ewgts,
//                                          int* ierr) {
//     local_graph& graph(*reinterpret_cast<local_graph*>(data));
//     *ierr = ZOLTAN_OK;
    
//     assert(num_gid_entries == 1);
//     assert(num_lid_entries == 1);
//     assert(num_obj == graph.num_local_objects());

//     // set the num_edges field for all local edges
//     for(int i = 0, sum_index = 0; i < graph.num_local_objects(); ++i) {
//       assert(num_edges[i] == graph.num_edges(local_ids[i]));
//       for(size_t j = 0; j < graph.edges[i].size(); ++j, ++sum_index) {
//         nbor_global_id[sum_index] = graph.edges[i][j];
//         nbor_procs[sum_index] = graph.owning_proc(nbor_global_id[sum_index]);
//       }
//     }
//   } 


//   static void zoltan_obj_size_multi_fun(void* data,
//                                         int num_gid_entries,
//                                         int num_lid_entries,
//                                         int num_ids,
//                                         ZOLTAN_ID_PTR global_ids,
//                                         ZOLTAN_ID_PTR local_ids, 
//                                         int* sizes,
//                                         int* ierr) {
//     local_graph& graph(*reinterpret_cast<local_graph*>(data));
//     *ierr = ZOLTAN_OK;
    
//     assert(num_gid_entries == 1);
//     assert(num_lid_entries == 1);
//     assert(num_ids == graph.num_local_objects());
//     for(int i = 0; i < graph.num_local_objects(); ++i) {
//       sizes[i] = graph.vdata[i].size();
//     }
//   }

//   static void zoltan_pack_obj_fun(void* data,
//                                   int num_gid_entries,
//                                   int num_lid_entries,
//                                   ZOLTAN_ID_PTR global_id,
//                                   ZOLTAN_ID_PTR local_id, 
//                                   int dest,
//                                   int size,
//                                   char* buffer,
//                                   int* ierr) {
//     local_graph& graph(*reinterpret_cast<local_graph*>(data));   
    
//     assert(num_gid_entries == 1);
//     assert(num_lid_entries == 1);
//     assert(global_id != NULL);
//     assert(local_id != NULL);
    
//     int lid = *local_id;
//     assert(lid > 0);
//     assert(size_t(lid) < graph.vdata.size());

//     // Do the copy
//     memcpy(buffer, graph.vdata[lid].c_str(), size);
    
//     *ierr = ZOLTAN_OK;

//   }






// };


// int main(int argc, char** argv) {
//   // Initialize MPI comm layer
//   MPI::Init(argc, argv);
//   int mpi_rank = MPI::COMM_WORLD.Get_rank();
//   int mpi_size = MPI::COMM_WORLD.Get_size();
//   if(mpi_rank == 0) {
//     std::cout << "Size: " << mpi_size << std::endl;
//   }

//   // Initialize Zoltan
//   float zoltan_version;
//   int error = Zoltan_Initialize(argc, argv, &zoltan_version);
//   if(error != ZOLTAN_OK) {
//     std::cout << "Unable to launch Zoltan!" << std::endl;
//     return EXIT_FAILURE;
//   }
//   if(mpi_rank == 0)
//     std::cout << "Zoltan version: " << zoltan_version << std::endl;

//   //  std::cout << "Loading file: " << argv[1] << std::endl;
//   local_graph graph(mpi_rank, mpi_size, argv[1]);


//   // Create the Zoltan object
//   Zoltan* zolt_ptr = new Zoltan(MPI::COMM_WORLD);
//   assert(zolt_ptr != NULL);
//   Zoltan& zolt(*zolt_ptr);


//   // Set the partitioning method to graph
//   error = zolt.Set_Param("LB_METHOD", "GRAPH");
//   assert(error == ZOLTAN_OK);

//   error = zolt.Set_Param( "NUM_GID_ENTRIES", "1");  /* global ID is 1 integer */
//   assert(error == ZOLTAN_OK);

//   error = zolt.Set_Param( "NUM_LID_ENTRIES", "1");  /* local ID is 1 integer */
//   assert(error == ZOLTAN_OK);

//   error = zolt.Set_Param( "OBJ_WEIGHT_DIM", "0");   /* we omit object weights */
//   assert(error == ZOLTAN_OK);

//   // Set the package to parmetis
//   // For more details see:
//   // http://www.cs.sandia.gov/zoltan/ug_html/ug_alg_parmetis.html
//   error = zolt.Set_Param("GRAPH_PACKAGE", "Parmetis");
//   assert(error == ZOLTAN_OK);
//   error = zolt.Set_Param("LB_APPROACH", "PARTITION");
//   assert(error == ZOLTAN_OK);

//   // http://www.cs.sandia.gov/zoltan/ug_html/ug_alg.html#RETURN_LISTS
//   if(mpi_rank == 0 && argc >= 3) {
//     std::cout << "Number of parts: " << argv[2] << std::endl;  
//     error = zolt.Set_Param("NUM_GLOBAL_PARTS", argv[2]);
//     assert(error == ZOLTAN_OK);
//   }

//   error = zolt.Set_Param("RETURN_LISTS", "PARTS");
//   assert(error == ZOLTAN_OK);


//   // Auto migration?
//   error = zolt.Set_Param("AUTO_MIGRATE", "TRUE");
//   assert(error == ZOLTAN_OK);



//   // Set functions
//   error = zolt.Set_Num_Obj_Fn(local_graph::zoltan_num_obj_fun, &graph);
//   assert(error == ZOLTAN_OK);

//   error = zolt.Set_Obj_List_Fn(local_graph::zoltan_obj_list_fun, &graph);
//   assert(error == ZOLTAN_OK);

//   error = zolt.Set_Num_Edges_Multi_Fn(local_graph::
//                                       zoltan_num_edges_multi_fun, &graph);
//   assert(error == ZOLTAN_OK);

//   error = zolt.Set_Edge_List_Multi_Fn(local_graph::
//                                       zoltan_edge_list_multi_fun, &graph);
//   assert(error == ZOLTAN_OK);



//   // Do the partitioning 
//   int changes = -1;
//   int num_gid_entries = -1;
//   int num_lid_entries = -1;
//   int num_import = -1;
//   ZOLTAN_ID_PTR import_global_ids = NULL;
//   ZOLTAN_ID_PTR import_local_ids = NULL;
//   int* import_procs = NULL;
//   int* import_to_part = NULL;
//   int num_export = -1;
//   ZOLTAN_ID_PTR export_global_ids = NULL;
//   ZOLTAN_ID_PTR export_local_ids = NULL;
//   int* export_procs = NULL;
//   int* export_to_part = NULL;

//   error = 
//     zolt.LB_Partition(changes,
//                       num_gid_entries,
//                       num_lid_entries,
//                       num_import,
//                       import_global_ids,
//                       import_local_ids,
//                       import_procs,
//                       import_to_part,
//                       num_export,
//                       export_global_ids,
//                       export_local_ids,
//                       export_procs,
//                       export_to_part);
//   assert(error == ZOLTAN_OK);


//   // Do something with the partitioning 
//   if(mpi_rank == 1) {
//     assert(export_global_ids != NULL);
//     assert(export_local_ids != NULL);
//     assert(export_procs != NULL);
//     assert(export_to_part != NULL);
//     for(int i = 0; i < num_export; ++i) {
//       std::cout << export_global_ids[i] << '\t'
//                 << export_local_ids[i] << '\t'
//                 << export_procs[i] << '\t'
//                 << export_to_part[i] << std::endl;
//     }
//   }
  


//   // Free the arrays after doing something
//   error = zolt.LB_Free_Part(&import_global_ids,
//                             &import_local_ids,
//                             &import_procs,
//                             &import_to_part);
//   assert(error == ZOLTAN_OK);

//   error = zolt.LB_Free_Part(&export_global_ids,
//                             &export_local_ids,
//                             &export_procs,
//                             &export_to_part);
//   assert(error == ZOLTAN_OK);


//   // Destroy the zoltan object
//   delete zolt_ptr;

//   // Kill the mpi layer
//   MPI::Finalize();
//   return EXIT_SUCCESS;
// }



















