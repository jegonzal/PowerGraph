




#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <graphlab.hpp>

#include <mpi.h>
#include <zoltan_cpp.h>


#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/filesystem.hpp>


#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

#include <graphlab/distributed2/graph/partitioning/adjacency_list.hpp>
#include <graphlab/distributed2/graph/partitioning/partitioning_tools.hpp>



#include <graphlab/macros_def.hpp>

using namespace graphlab;



struct vertex_info {
  vertex_id_t vid;
  uint16_t atom_id;
};
SERIALIZABLE_POD(vertex_info);




struct graph_data {
  adjacency_list alist;
  std::vector<vertex_id_t> vertex2cpu;
  vertex_id_t nverts;
  
  graph_data(const std::string& path,
             const std::vector< std::string >& fnames) {
    std::cout << "Loading " << fnames.size() << " files." << std::endl;
    // load all the adjacency lists
    for(size_t i = 0; i < fnames.size(); ++i) {
      std::string abs_fname(path + "/" + fnames[i]);
      std::cout << "Loading: " << abs_fname << std::endl;
      // Merge all the adjacency info
      alist.load(abs_fname);
    }
    std::cout << "Finished Load!" << std::endl;
    const size_t ROOT_NODE(0);
    // Compute the initial fragmentation
    if(mpi_tools::rank() == ROOT_NODE) {
      std::vector< std::vector<vertex_id_t> > cpu2vertex;
      mpi_tools::gather(alist.local_vertices, cpu2vertex);

      std::cout << "Gather on root info: " << std::endl;
      for(size_t i = 0; i < cpu2vertex.size(); ++i) {
        std::cout << "Proc " << i << ": " 
                  << cpu2vertex[i].size() << std::endl;
      }

      // Determine the maximum  vertex_id_value
      vertex_id_t maxid(0);
      size_t vertices(0);
      for(size_t i = 0; i < cpu2vertex.size(); ++i) {
        for(size_t j = 0; j < cpu2vertex[i].size(); ++j) { 
          vertices++;
          maxid = std::max(maxid, cpu2vertex[i][j]);
        }
      }
      std::cout << "Vertices:\t" << vertices << std::endl;
      std::cout << "Maxid:\t" << maxid << std::endl;
      assert(vertices > 0);
      assert(vertices == maxid+1);
      vertex2cpu.resize(vertices, vertex_id_t(-1));
      // Fill out the map
      for(size_t i = 0; i < cpu2vertex.size(); ++i) {
        for(size_t j = 0; j < cpu2vertex[i].size(); ++j) {
          vertex_id_t vid = cpu2vertex[i][j];
          assert(vid < vertex2cpu.size());
          assert(vertex2cpu[vid] == vertex_id_t(-1));
          vertex2cpu[vid] = i;
        }
      }
      // Check that all vertices were assigned to a cpu
      for(size_t i = 0; i < vertex2cpu.size(); ++i) 
        assert(vertex2cpu[i] != vertex_id_t(-1));
    } else {
      mpi_tools::gather(ROOT_NODE, alist.local_vertices);
    }
    // scatter the map to all machines
    mpi_tools::bcast(ROOT_NODE, vertex2cpu);
    nverts = vertex2cpu.size();
    
    // A fairly costly ceck of local structures
    std::cout << "Checking local structures." << std::endl;
    alist.check_local_structures(nverts);
    std::cout << "Finished." << std::endl;
  } // end of constructor            
}; // end of graph_data




void compute_local_fnames(std::vector<std::string>& fnames) {
  size_t mpi_rank = graphlab::mpi_tools::rank();
  size_t mpi_size = graphlab::mpi_tools::size();
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






void save_partition_results(const size_t num_vertices,
                            const int num_export,
                            const ZOLTAN_ID_PTR export_global_ids,
                            const int* export_to_part,
                            const std::string& path) {

  std::vector< vertex_info > local_vertex_info(num_export);
  for(size_t i = 0; i < local_vertex_info.size(); ++i) {
    vertex_info& vinfo(local_vertex_info[i]);
    vinfo.vid = export_global_ids[i];
    vinfo.atom_id = export_to_part[i];
  }

  
  // build up the local vertex info table
  std::vector< std::vector<vertex_info> > all_vertex_info;
  if(mpi_tools::rank() == 0) { //root code
    mpi_tools::gather(local_vertex_info, all_vertex_info);
  } else {
    mpi_tools::gather(0, local_vertex_info);
  }


  // compute the masters
  std::set<size_t> master_ranks;
  mpi_tools::get_master_ranks(master_ranks);
  assert(master_ranks.count(0) > 0);

  // Construct global atom table
  std::vector<uint16_t> global_atom_table;  

  // on node zero gather up a single giant vector
  if(mpi_tools::rank() == 0) {
    assert(all_vertex_info.size() == mpi_tools::size());
    // Construct global atom table
    global_atom_table.resize(num_vertices);
    for(size_t i = 0; i < all_vertex_info.size(); ++i) {
      for(size_t j = 0; j < all_vertex_info[i].size(); ++j) {
        global_atom_table[all_vertex_info[i][j].vid] = 
          all_vertex_info[i][j].atom_id;
      }
    }
    // Send the giant vector to the master ranks
    foreach(size_t dest, master_ranks) {
      if(dest != 0) {    
        mpi_tools::send(global_atom_table, dest);
      }
    }   
  } else if(master_ranks.count(mpi_tools::rank()) > 0) {
    mpi_tools::recv(global_atom_table, 0);
  }

  // save the master files
  if(!global_atom_table.empty()) {
    assert(master_ranks.count(mpi_tools::rank()) > 0);
    // save a raw text file of the global atom table
    std::string absfname = path + "/partitioning.txt";
    std::ofstream fout(absfname.c_str());
    for(size_t i = 0; i < global_atom_table.size(); ++i) {
      fout << global_atom_table[i] << '\n';
    }
    fout.close();
  }
}









int zoltan_num_obj_fun(void* data, int* ierr) {
  assert(data != NULL);
  const graph_data& zgdata(*reinterpret_cast<graph_data*>(data));
  *ierr = ZOLTAN_OK;  
  return zgdata.alist.local_vertices.size();
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
  const graph_data& zgdata(*reinterpret_cast<graph_data*>(data));
  *ierr = ZOLTAN_OK;
  assert(num_gid_entries == 1);
  assert(num_lid_entries == 0);
  assert(local_ids == NULL);
  assert(obj_wgts == NULL);

  for(size_t i = 0; i < zgdata.alist.local_vertices.size(); ++i) {
    global_ids[i] = zgdata.alist.local_vertices.at(i);
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
  const graph_data& zgdata(*reinterpret_cast<graph_data*>(data));
  *ierr = ZOLTAN_OK;
  assert(num_gid_entries == 1);
  assert(num_lid_entries == 0);
  assert(local_ids == NULL);
  assert(num_objs >= 0);
  for(size_t i = 0; i < zgdata.alist.in_neighbor_ids.size(); ++i) {
    assert(i < size_t(num_objs));
    assert(i < zgdata.alist.local_vertices.size());
    assert(global_ids[i] == zgdata.alist.local_vertices[i]);
    num_edges[i] = zgdata.alist.in_neighbor_ids.at(i).size();
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
  const graph_data& zgdata(*reinterpret_cast<graph_data*>(data));
  *ierr = ZOLTAN_OK;
  assert(num_gid_entries == 1);
  assert(num_lid_entries == 0);
  assert(local_ids == NULL);
  assert(num_objs >= 0);
  size_t sum_num_edges = 0;
  size_t eindex = 0;
  for(size_t i = 0; i < zgdata.alist.in_neighbor_ids.size(); ++i) {
    assert(i < size_t(num_objs));
    assert(i < zgdata.alist.local_vertices.size());
    assert(global_ids[i] == zgdata.alist.local_vertices[i]);
    assert(size_t(num_edges[i]) == zgdata.alist.in_neighbor_ids[i].size());
    sum_num_edges += num_edges[i];
    assert(num_edges[i] == zgdata.alist.in_neighbor_ids[i].size());
    for(size_t j = 0; j < zgdata.alist.in_neighbor_ids[i].size(); ++j, ++eindex) {
      vertex_id_t vid = zgdata.alist.in_neighbor_ids[i][j];
      vertex_id_t cpu = zgdata.vertex2cpu[vid];
      assert(cpu < mpi_tools::size());
      nbor_global_id[eindex] = vid;
      nbor_cpus[eindex] = cpu;
    }
  }
  assert(sum_num_edges == eindex); 
}




void graphlab::partitioning_tools::
construct_partitioning(int argc, char** argv,
                       int numparts,
                       const std::string& path) {
  // Get the mpi rank and size
  size_t mpi_rank = graphlab::mpi_tools::rank();
  //size_t mpi_size = graphlab::mpi_tools::size();



  // Load the filenames
  std::vector<std::string> fnames;
  adjacency_list::list_vlist_files(path, fnames);
  compute_local_fnames(fnames);

  // construct the local graph data
  graph_data zgdata(path, fnames);
  
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

  error = zolt.Set_Param("DEBUG_LEVEL", "1");
  assert(error == ZOLTAN_OK);

  //  error = zolt.Set_Param("IMBALANCE_TOL", "2");
  //  assert(error == ZOLTAN_OK);


  error = zolt.Set_Param("GRAPH_SYMMETRIZE", "TRANSPOSE");
  assert(error == ZOLTAN_OK);

  error = zolt.Set_Param("CHECK_GRAPH", "2");
  assert(error == ZOLTAN_OK);

  error = zolt.Set_Param("GRAPH_BUILD_TYPE", "NORMAL");
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




  // Do coloring on machine zero online
  if(false){
    //     error = zolt.Set_Param("COLORING_PROBLEM", "distance-2");
    //     assert(error == ZOLTAN_OK);
    
    // error = zolt.Set_Param("COMM_PATTERN", "A");
    // assert(error == ZOLTAN_OK);

    size_t super_step_size = 
      std::max(size_t(100),
               zgdata.nverts / mpi_tools::size());
    std::stringstream strm;
    strm << super_step_size;
    std::string super_step_size_str(strm.str());
    

    // error = zolt.Set_Param("SUPERSTEP_SIZE", super_step_size_str.c_str());
    // assert(error == ZOLTAN_OK);


    std::set<size_t> master_ranks;
    mpi_tools::get_master_ranks(master_ranks);
    if( master_ranks.count(mpi_tools::rank()) > 0) {
      int nverts = zgdata.nverts;
      int num_gid_entries = 1;
      std::vector<ZOLTAN_ID_TYPE> globalids(nverts);
      std::vector<int> colors(nverts, -1);
      for(size_t i = 0; i < globalids.size(); ++i) 
        globalids[i] = i;
      error = zolt.Color(num_gid_entries, 
                         nverts,
                         &(globalids[0]),
                         &(colors[0]));
      std::string absfname = path + "/coloring.txt";
      std::ofstream fout(absfname.c_str());
      for(size_t i = 0; i < colors.size(); ++i) {
        assert(colors[i] >= 1);
        --colors[i];
        fout << colors[i] << '\n';
      }
      fout.close();
    } else {
//       int ierr(-1);
//       int nverts = zoltan_num_obj_fun(&zgdata, &ierr);
//       std::cout << "local objects: " << nverts << std::endl;
//       std::vector<ZOLTAN_ID_TYPE> globalids(nverts);
//       std::vector<int> colors(nverts, -1);
//       zoltan_obj_list_fun(&zgdata, 1, 0, &(globalids[0]), NULL, 0, NULL, &ierr);
//       int num_gid_entries = 1;
//       error = zolt.Color(num_gid_entries, 
//                          nverts,
//                          &(globalids[0]), 
//                          &(colors[0]));
      int num_gid_entries = 1;
      int nverts = 0;
      error = zolt.Color(num_gid_entries, 
                         nverts,
                         NULL, NULL);


    }
  }
                     


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
  

  save_partition_results(zgdata.nverts,
                         num_export,
                         export_global_ids,
                         export_to_part,
                         path);



  // // Export the partition file on each machine
  // {
  //   if(mpi_rank == 0)
  //     std::cout << "Saving partitioning." << std::endl;
  //   std::string partition_filename;
  //   std::stringstream strm;
  //   strm << path << "/partition_" << std::setw(3) << std::setfill('0')
  //        << mpi_rank << ".txt";
  //   partition_filename = strm.str();
  //   std::ofstream fout(partition_filename.c_str());
  //   assert(fout.good());
  //   for(int i = 0; i < num_export; ++i) {
  //     fout << export_global_ids[i] << '\t'
  //          << export_to_part[i] << '\n';
  //   }
  //   fout.close();
  //   if(mpi_rank == 0)
  //     std::cout << "Finished saving partitioning." << std::endl;
  // }

  // Evaluate the partitioning
//   {     
//     if(mpi_rank == 0)
//       std::cout << "Evaluating Partitioning" << std::endl;
//     const int print_stats(1);
//     ZOLTAN_GRAPH_EVAL graph_info;
//     // http://www.cs.sandia.gov/zoltan/ug_html/ug_interface_lb.html#Zoltan_LB_Eval
//     error = zolt.LB_Eval_Graph(print_stats,   &graph_info);
//     assert(error == ZOLTAN_OK);
// //     if(mpi_rank == 0) {
// //       std::cout << "Imbalance: " << graph_info.obj_imbalance
// //                 << std::endl;
// //       std::cout << "Cuts:\t";
// //       for(size_t i = 1; i < EVAL_SIZE; ++i) {
// //         std::cout << graph_info.cuts[i] << '\t';
// //       }
// //       std::cout << std::endl;
// //       std::cout << "nnborparts:\t";
// //       for(size_t i = 1; i < EVAL_SIZE; ++i) {
// //         std::cout << graph_info.nnborparts[i] << '\t';
// //       }
// //       std::cout << std::endl;
// //     }
    
//     if(mpi_rank == 0)
//       std::cout << "Finished evaluating partitioning." << std::endl;
//   }





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



















































































