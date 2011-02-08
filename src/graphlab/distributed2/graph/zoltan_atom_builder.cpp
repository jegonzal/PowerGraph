


#include <sys/types.h>
#include <ifaddrs.h>
#include <netinet/in.h>



#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <graphlab.hpp>
#include <boost/filesystem.hpp>

#include <mpi.h>
#include <zoltan_cpp.h>


#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/filesystem.hpp>


#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/distributed2/graph/graph_fragment.hpp>
#include <graphlab/distributed2/graph/zoltan_atom_builder.hpp>
#include <graphlab/util/mpi_tools.hpp>


#include <graphlab/macros_def.hpp>

using namespace graphlab;

typedef std::vector<graph_fragment::file_description> 
description_vec_type;

typedef std::vector<graph_fragment::structure_description> 
structure_vec_type;

typedef std::map<vertex_id_t, vertex_id_t> part2cpu_type;


struct vertex_info {
  vertex_id_t vid;
  uint16_t atom_id;
};


struct graph_data {
  description_vec_type desc_vec;
  structure_vec_type structure_vec;
  part2cpu_type part2cpu;
};

void list_structure_files(const std::string& pathname, 
                          std::vector<std::string>& files) {
  namespace fs = boost::filesystem;
  fs::path path(pathname);
  assert(fs::exists(path));
  for(fs::directory_iterator iter( path ), end_iter; 
      iter != end_iter; ++iter) {
    if( ! fs::is_directory(iter->status()) ) {
      std::string filename(iter->path().filename());
      size_t period_loc = filename.rfind('.');
      if(period_loc != std::string::npos) {
        std::string ending(filename.substr(period_loc));
        if(ending == graph_fragment::structure_suffix) {
          files.push_back(iter->path().filename());
        }
      }
    }
  }
  std::sort(files.begin(), files.end());
} // end of list_structure_file


void load_file_descs(const std::string& path,
                     const std::vector<std::string>& fnames,
                     description_vec_type& file_desc_vec) {
  file_desc_vec.resize(fnames.size());
  for(size_t i = 0; i < fnames.size(); ++i) {
    std::string abs_filename(path + "/" + fnames[i]);
    std::ifstream fin(abs_filename.c_str(), std::ios::binary | std::ios::in);
    assert(fin.good());
    graphlab::iarchive iarc(fin);
    iarc >> file_desc_vec[i];
    fin.close();
  }
} // end of load structures




void load_structures(const std::string& path,
                     const std::vector<std::string>& structure_fnames,
                     structure_vec_type& structures) {
  structures.resize(structure_fnames.size());
  for(size_t i = 0; i < structure_fnames.size(); ++i) {
    std::string abs_filename(path + "/" + structure_fnames[i]);
    std::cout << abs_filename << std::endl;
    std::ifstream fin(abs_filename.c_str(), std::ios::binary | std::ios::in);
    assert(fin.good());
    graphlab::iarchive iarc(fin);
    graph_fragment::file_description desc;
    iarc >> desc >> structures[i];
    fin.close();
  }
} // end of load structures











void compute_local_fnames(std::vector<std::string>& fnames) {
  int mpi_rank(-1), mpi_size(-1);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  assert(mpi_rank >= 0);
  assert(mpi_size >= 0);
  size_t num_cpus(mpi_size);

  std::vector< std::vector<std::string> > cpu2fnames;
  mpi_tools::all_gather(fnames, cpu2fnames);
  assert(num_cpus == cpu2fnames.size());

 
  // Compute set of all fnames
  std::set<std::string> unassigned_fnames;
  for(size_t cpuid = 0; cpuid < num_cpus; ++cpuid) {
    for(size_t j = 0; j < cpu2fnames[cpuid].size(); ++j) {
      unassigned_fnames.insert(cpu2fnames[cpuid][j]);
    }    
  }

  // Assign fnames to machines that own them
  fnames.clear();
  for(size_t cpuid = 0; !unassigned_fnames.empty(); 
      cpuid = (cpuid + 1) % num_cpus) {
    assert(cpuid < cpu2fnames.size());
    while(!cpu2fnames[cpuid].empty()) {
      std::string lastfname = cpu2fnames[cpuid].back();
      cpu2fnames[cpuid].pop_back();
      if(unassigned_fnames.count(lastfname) > 0) {
        unassigned_fnames.erase(lastfname);
        if(cpuid == size_t(mpi_rank)) {
          fnames.push_back(lastfname);
        }        
        break;
      }
      
    }
  }
  assert(unassigned_fnames.empty());
}




uint32_t get_local_ip() {
  uint32_t ip;
  // code adapted from
  struct ifaddrs * ifAddrStruct = NULL;
  getifaddrs(&ifAddrStruct);
  struct ifaddrs * firstifaddr = ifAddrStruct;
  ASSERT_NE(ifAddrStruct, NULL);
  while (ifAddrStruct != NULL) {
    if (ifAddrStruct->ifa_addr != NULL && 
        ifAddrStruct->ifa_addr->sa_family == AF_INET) {
      char* tmpAddrPtr = NULL;
        // check it is IP4 and not lo0.
      tmpAddrPtr = (char*)&((struct sockaddr_in *)ifAddrStruct->ifa_addr)->sin_addr;
      ASSERT_NE(tmpAddrPtr, NULL);
      if (tmpAddrPtr[0] != 127) {
        memcpy(&ip, tmpAddrPtr, 4);
        break;
      }
      //break;
    }
    ifAddrStruct=ifAddrStruct->ifa_next;
  }
  freeifaddrs(firstifaddr);
  return ip;
}


void gather_partition_results(const std::vector<int>& root_rank,
                              const int num_export,
                              const ZOLTAN_ID_PTR export_global_ids,
                              const int* export_to_part) {
  std::vector< vertex_info > local_vertex_info(num_export);
  for(size_t i = 0; i < local_vertex_info.size(); ++i) {
    vertex_info& vinfo(local_vertex_info[i]);
    vinfo.vid = export_global_ids[i];
    vinfo.atom_id = export_to_part[i];
  }
  
}




void compute_part2cpu(const description_vec_type& local_files,
                      part2cpu_type& part2cpu) {

  std::vector<vertex_id_t> local_parts(local_files.size());
  for(size_t i = 0; i < local_files.size(); ++i)
    local_parts[i] = local_files[i].id;

  std::vector< std::vector<vertex_id_t> > cpu2parts;
  mpi_tools::all_gather(local_parts, cpu2parts);
  for(size_t i = 0; i < cpu2parts.size(); ++i) {
    for(size_t j = 0; j < cpu2parts[i].size(); ++j) {
      part2cpu[cpu2parts[i][j]] = i;
    }
  }
}






int zoltan_num_obj_fun(void* data, int* ierr) {
  assert(data != NULL);
  graph_data& zgdata(*reinterpret_cast<graph_data*>(data));
  int num_objects(0);
  foreach(const graph_fragment::file_description& desc, zgdata.desc_vec) {
    num_objects += desc.num_local_verts;
  }
  *ierr = ZOLTAN_OK;
  return num_objects;
} // end of zoltan_num_obj_fun




void zoltan_obj_list_fun(void* data, 
                         int num_gid_entries,
                         int num_lid_entries,
                         ZOLTAN_ID_PTR global_ids,
                         ZOLTAN_ID_PTR local_ids,
                         int wgt_dim,
                         float* obj_wgts,
                         int* ierr) {
  assert(data != NULL);
  graph_data& zgdata(*reinterpret_cast<graph_data*>(data));
  *ierr = ZOLTAN_OK;
  assert(num_gid_entries == 1);
  assert(num_lid_entries == 0);
  assert(local_ids == NULL);
  assert(obj_wgts == NULL);
  size_t index(0);
  foreach(const graph_fragment::structure_description& structure, 
          zgdata.structure_vec) {
    for(vertex_id_t i = 0; i < structure.desc.num_local_verts; ++i, ++index) {
      global_ids[index] = i + structure.desc.begin_vertex;
    }
  }
}




static void zoltan_num_edges_multi_fun(void* data,
                                       int num_gid_entries,
                                       int num_lid_entries,
                                       int num_objs,
                                       ZOLTAN_ID_PTR global_ids,
                                       ZOLTAN_ID_PTR local_ids, 
                                       int* num_edges, 
                                       int* ierr) {
  assert(data != NULL);
  graph_data& zgdata(*reinterpret_cast<graph_data*>(data));
  *ierr = ZOLTAN_OK;
  assert(num_gid_entries == 1);
  assert(num_lid_entries == 0);
  assert(local_ids == NULL);
  assert(num_objs >= 0);
  size_t index(0);
  foreach(const graph_fragment::structure_description& structure, zgdata.structure_vec) {
    for(vertex_id_t i = 0; i < structure.desc.num_local_verts; ++i, ++index) {
      assert(index < size_t(num_objs));
      assert(global_ids[index] == i + structure.desc.begin_vertex);
      num_edges[index] = structure.neighbor_ids[i].size();
    }
  }
} 






static void zoltan_edge_list_multi_fun(void* data,
                                       int num_gid_entries,
                                       int num_lid_entries,
                                       int num_objs,
                                       ZOLTAN_ID_PTR global_ids,
                                       ZOLTAN_ID_PTR local_ids, 
                                       int* num_edges, 
                                       ZOLTAN_ID_PTR nbor_global_id,
                                       int* nbor_cpus,
                                       int wgt_dim,
                                       float* ewgts,
                                       int* ierr) {
  assert(data != NULL);
  graph_data& zgdata(*reinterpret_cast<graph_data*>(data));
  *ierr = ZOLTAN_OK;
  assert(num_gid_entries == 1);
  assert(num_lid_entries == 0);
  assert(local_ids == NULL);
  assert(num_objs >= 0);

  size_t vindex(0), eindex(0);
  foreach(const graph_fragment::structure_description& structure, 
          zgdata.structure_vec) {
    for(vertex_id_t i = 0; i < structure.desc.num_local_verts; ++i, ++vindex) {
      assert(vindex < size_t(num_objs));
      assert(global_ids[vindex] == i + structure.desc.begin_vertex);
      assert(num_edges[vindex] >= 0);
      assert(size_t(num_edges[vindex]) == structure.neighbor_ids[i].size());
      // save all the edges
      for(vertex_id_t j = 0; 
          j < structure.neighbor_ids[i].size(); ++j, ++eindex) {
        nbor_global_id[eindex] = structure.neighbor_ids[i][j];
        // compute the owning cpuessor
        nbor_cpus[eindex] = 
          zgdata.part2cpu[structure.desc.owning_fragment( nbor_global_id[eindex] ) ];

      }
    }
  }
}




void graphlab::construct_partitioning(int argc, char** argv,
                                      int numparts,
                                      const std::string& path) {
  // Get the mpi rank and size
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


  // Load the filenames
  std::vector<std::string> fnames;
  list_structure_files(path, fnames);
  compute_local_fnames(fnames);

  // global structure with all necessary information
  graph_data zgdata;
  

  // Load all the descriptions for the filenames we own
  load_file_descs(path, fnames, zgdata.desc_vec);
  
  // compute part2cpus map
  compute_part2cpu(zgdata.desc_vec, zgdata.part2cpu);

  // load the actual structure files
  load_structures(path, fnames, zgdata.structure_vec);

  // Initialize Zoltan
  float zoltan_version;
  int error = Zoltan_Initialize(argc, argv, &zoltan_version);
  if(error != ZOLTAN_OK) {
    std::cout << "Unable to launch Zoltan!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(mpi_rank == 0)
    std::cout << "Zoltan version: " << zoltan_version << std::endl;

  
  // Create the Zoltan object
  Zoltan* zolt_ptr = new Zoltan(MPI::COMM_WORLD);
  assert(zolt_ptr != NULL);
  Zoltan& zolt(*zolt_ptr);


  // Set the partitioning method to graph
  error = zolt.Set_Param("LB_METHOD", "GRAPH");
  assert(error == ZOLTAN_OK);

  error = zolt.Set_Param("DEBUG_LEVEL", "0");
  assert(error == ZOLTAN_OK);


  error = zolt.Set_Param( "NUM_GID_ENTRIES", "1");  /* global ID is 1 integer */
  assert(error == ZOLTAN_OK);

  error = zolt.Set_Param( "NUM_LID_ENTRIES", "0");  /* local ID is 1 integer */
  assert(error == ZOLTAN_OK);

  error = zolt.Set_Param( "OBJ_WEIGHT_DIM", "0");   /* we omit object weights */
  assert(error == ZOLTAN_OK);

  // Set the package to parmetis
  // For more details see:
  // http://www.cs.sandia.gov/zoltan/ug_html/ug_alg_parmetis.html
  error = zolt.Set_Param("GRAPH_PACKAGE", "Parmetis");
  assert(error == ZOLTAN_OK);
  error = zolt.Set_Param("LB_APPROACH", "PARTITION");
  assert(error == ZOLTAN_OK);

  // http://www.cs.sandia.gov/zoltan/ug_html/ug_alg.html#RETURN_LISTS 
  {
    std::stringstream strm;
    strm << numparts;
    std::string str(strm.str());
    if(mpi_rank == 0) std::cout << "Parts: " << str << std::endl;
    error = zolt.Set_Param("NUM_GLOBAL_PARTS", str.c_str());
    assert(error == ZOLTAN_OK);
  } 

  
  // Documentation for return lists
  // http://www.cs.sandia.gov/zoltan/ug_html/ug_alg.html#RETURN_LISTS
  error = zolt.Set_Param("RETURN_LISTS", "PARTS");
  assert(error == ZOLTAN_OK);

//   // Auto migration?
//   error = zolt.Set_Param("AUTO_MIGRATE", "TRUE");
//   assert(error == ZOLTAN_OK);



  // Set functions
  error = zolt.Set_Num_Obj_Fn(zoltan_num_obj_fun, &zgdata);
  assert(error == ZOLTAN_OK);
  
  error = zolt.Set_Obj_List_Fn(zoltan_obj_list_fun, &zgdata);
  assert(error == ZOLTAN_OK);
  
  error = zolt.Set_Num_Edges_Multi_Fn(zoltan_num_edges_multi_fun, &zgdata);
  assert(error == ZOLTAN_OK);
  
  error = zolt.Set_Edge_List_Multi_Fn(zoltan_edge_list_multi_fun, &zgdata);
  assert(error == ZOLTAN_OK);


  ///////////////////////////////////////////////////////////////////////
  // Do the partitioning 
  //////
  // Documentation:
  // http://www.cs.sandia.gov/zoltan/ug_html/ug_interface_lb.html#Zoltan_LB_Partition
  
  int changes = -1;
  int num_gid_entries = -1;
  int num_lid_entries = -1;
  int num_import = -1;
  ZOLTAN_ID_PTR import_global_ids = NULL;
  ZOLTAN_ID_PTR import_local_ids = NULL;
  int* import_cpus = NULL;
  int* import_to_part = NULL;
  int num_export = -1;
  ZOLTAN_ID_PTR export_global_ids = NULL;
  ZOLTAN_ID_PTR export_local_ids = NULL;
  int* export_cpus = NULL;
  int* export_to_part = NULL;

  if(mpi_rank == 0) 
    std::cout << "Running partitioner." << std::endl;

  error = 
    zolt.LB_Partition(changes,
                      num_gid_entries,
                      num_lid_entries,
                      num_import,
                      import_global_ids,
                      import_local_ids,
                      import_cpus,
                      import_to_part,
                      num_export,
                      export_global_ids,
                      export_local_ids,
                      export_cpus,
                      export_to_part);
  assert(error == ZOLTAN_OK);

  if(mpi_rank == 0) 
    std::cout << "Finished Successfully." << std::endl;
  

  // Check return values
  assert(changes == true);
  assert(num_gid_entries == 1);
  assert(num_lid_entries == 0);
  
  assert(num_import == -1);
  assert(import_global_ids == NULL);
  assert(import_local_ids == NULL);
  assert(import_cpus == NULL);
  assert(import_to_part == NULL);
  

  assert(export_global_ids != NULL);
  assert(export_local_ids == NULL);
  assert(export_cpus != NULL);
  assert(export_to_part != NULL);
  


  // Export the partition file on each machine
  {
    if(mpi_rank == 0)
      std::cout << "Saving partitioning." << std::endl;
    std::string partition_filename;
    std::stringstream strm;
    strm << path << "/partition_" << std::setw(3) << std::setfill('0')
         << mpi_rank << ".txt";
    partition_filename = strm.str();
    std::ofstream fout(partition_filename.c_str());
    assert(fout.good());
    for(int i = 0; i < num_export; ++i) {
      fout << export_global_ids[i] << '\t'
           << export_to_part[i] << '\n';
    }
    fout.close();
    if(mpi_rank == 0)
      std::cout << "Finished saving partitioning." << std::endl;
  }

  // Evaluate the partitioning
  {     
    if(mpi_rank == 0)
      std::cout << "Evaluating Partitioning" << std::endl;
    const int print_stats(1);
    ZOLTAN_GRAPH_EVAL graph_info;
    // http://www.cs.sandia.gov/zoltan/ug_html/ug_interface_lb.html#Zoltan_LB_Eval
    error = zolt.LB_Eval_Graph(print_stats,   &graph_info);
    assert(error == ZOLTAN_OK);
//     if(mpi_rank == 0) {
//       std::cout << "Imbalance: " << graph_info.obj_imbalance
//                 << std::endl;
//       std::cout << "Cuts:\t";
//       for(size_t i = 1; i < EVAL_SIZE; ++i) {
//         std::cout << graph_info.cuts[i] << '\t';
//       }
//       std::cout << std::endl;
//       std::cout << "nnborparts:\t";
//       for(size_t i = 1; i < EVAL_SIZE; ++i) {
//         std::cout << graph_info.nnborparts[i] << '\t';
//       }
//       std::cout << std::endl;
//     }
    
    if(mpi_rank == 0)
      std::cout << "Finished evaluating partitioning." << std::endl;
  }





  ///////////////////////////////////////////////////////////////////////
  // cleanup
  //////

  // Free the arrays after doing something
  error = zolt.LB_Free_Part(&import_global_ids,
                            &import_local_ids,
                            &import_cpus,
                            &import_to_part);
  assert(error == ZOLTAN_OK);

  error = zolt.LB_Free_Part(&export_global_ids,
                            &export_local_ids,
                            &export_cpus,
                            &export_to_part);
  assert(error == ZOLTAN_OK);



  


  // Destroy the zoltan object
  delete zolt_ptr;

} // build atom files



















































































