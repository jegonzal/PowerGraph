#ifndef FROM_DISTRIBUTED_GRAPH_INCLUDE

#warning "distributed_graph_incremental_loader.hpp should not be included directly."
#warning "You should include only distributed_graph.hpp"
#warning "I will fix this for you now, but don't do it again!"

#include <graphlab/distributed2/graph/distributed_graph.hpp>


#else
#include <graphlab/util/generics/shuffle.hpp>

  /** From the atoms listed in the atom index file, construct the local store
  * using all the atoms in the current partition. Uses the dump files.
  * This is possibly the most memory efficient loading method, requiring very
  * little excess memory.
  */
template <typename VertexData, typename EdgeData>
void distributed_graph<VertexData,EdgeData>::construct_local_fragment_playback(const atom_index_file &atomindex,
                                        std::vector<std::vector<size_t> > partitiontoatom,
                                        size_t curpartition,
                                        bool do_not_load_data) {
  timer loadtimer;
  loadtimer.start();
  // first make a map mapping atoms to machines
  // we will need this later
  std::vector<procid_t> atom2machine;
  for (size_t i = 0 ;i< partitiontoatom.size(); ++i) {
    for (size_t j = 0 ; j < partitiontoatom[i].size(); ++j) {
      if (atom2machine.size() <= partitiontoatom[i][j]) 
        atom2machine.resize(partitiontoatom[i][j] + 1);
      atom2machine[partitiontoatom[i][j]] = i;
    }
  }



  // for convenience take a reference to the list of atoms in this partition
  std::vector<size_t>& atoms_in_curpart = partitiontoatom[curpartition];
  // create the atom file readers
  logstream(LOG_INFO) << "Atoms on this machine: " << atoms_in_curpart.size() << std::endl;
  
  // initiate playback
  #pragma omp parallel for
  for (int i = 0;i < (int)(atoms_in_curpart.size()); ++i) {
    std::string fname = atomindex.atoms[atoms_in_curpart[i]].file;
    logstream(LOG_DEBUG) << "Loading file: " << fname << '\t'
                         << localstore.num_vertices() << std::endl;
    playback_dump(fname + ".dump", atoms_in_curpart[i], atom2machine, do_not_load_data);
  }
  
  // playback complete.
  // now to shuffle IDs around to rearrange all of the 
  // owned vertices to the start
  shuffle_local_vertices_to_start();
  construct_ghost_auxiliaries();
 /* std::cout << "Owned:";
  for (size_t i = 0;i < ownedvertices.size(); ++i) std::cout << ownedvertices[i] << "\t";
  std::cout << "\n";

  std::cout << "Ghost:";
  for (size_t i = 0;i < ghostvertices.size(); ++i) std::cout << ghostvertices[i] << "\t";
  std::cout << "\n";
*/
  rmi.barrier();
  logger(LOG_INFO, "Synchronizing ghost data...");
  // shuffle for all the ghost data
  push_all_owned_vertices_to_replicas();
  rmi.dc().full_barrier();
  logger(LOG_INFO, "vertices synchronized.");

  push_all_owned_edges_to_replicas();
  rmi.dc().full_barrier();
  logger(LOG_INFO, "edges synchronized.");

  logger(LOG_INFO, "Synchronization complete.");
  rmi.dc().barrier();
  logger(LOG_INFO, "Performing data verification.");
  for (size_t i = 0;i < localstore.num_vertices(); ++i) {
    ASSERT_EQ(localstore.vertex_version(i), 1);
  }
  for (size_t i = 0;i < localstore.num_edges(); ++i) {
    ASSERT_EQ(localstore.edge_version(i), 1);
  }
  
  logger(LOG_INFO, "Finalize");
  localstore.finalize();
  logger(LOG_INFO, "Load complete.");
  rmi.comm_barrier();
  std::cout << "Load complete in " << loadtimer.current_time() << std::endl;  
}


template <typename VertexData, typename EdgeData>
typename distributed_graph<VertexData,EdgeData>::vertex_id_type 
  distributed_graph<VertexData,EdgeData>::create_vertex_if_missing(vertex_id_type globalvid,
                                                                      procid_t machine,
                                                                      uint16_t sourceatom,
                                                                      bool overwritedata,
                                                                      const VertexData &vdata) {
    typename global2localvid_type::const_iterator iter = global2localvid.find(globalvid);
    // nope. Insert it.
    if (iter == global2localvid.end()) {
      // add the vertex
      vertex_id_type localvid = localstore.add_vertex(vdata);
      // update the vid mappings
      global2localvid[globalvid] = localvid;
      local2globalvid.push_back(globalvid);
      ASSERT_EQ(local2globalvid.size(), localvid+1);
      // update the owner tables
      localvid2owner.push_back(machine);
      localvid2atom.push_back(sourceatom);
      if (machine == rmi.procid()) {
        globalvid2owner.set(globalvid, rmi.procid());
      }
      return localvid;
    }
    else {
      if (overwritedata) {
        localstore.vertex_data(iter->second) = vdata;
      }
      return iter->second;
    }
}

template <typename VertexData, typename EdgeData>
void distributed_graph<VertexData,EdgeData>::playback_dump(std::string filename,
                                      size_t atomid,
                                      std::vector<procid_t> atom2machine,
                                      bool do_not_load_data) {

  std::ifstream in_file(filename.c_str(), std::ios::binary);

  boost::iostreams::filtering_stream<boost::iostreams::input> fin; 
  fin.push(boost::iostreams::zlib_decompressor());
  fin.push(in_file);
  // flush the commands
  iarchive iarc(fin);

  while(fin.good()) {
    char command;
    fin >> command;
    if (fin.fail()) break;
    if (command == 'a' || command == 'b') {
      // add vertex skip
      vertex_id_type vid; uint16_t owner;
      iarc >> vid >> owner;
      alldatalock.lock();
      create_vertex_if_missing(vid, atom2machine[owner], owner);
      alldatalock.unlock();
    } else if (command == 'c') {
      vertex_id_type vid; uint16_t owner; std::string data;
      iarc >> vid >> owner >> data;
      // deserialize it
      VertexData vd;
      if (data.size() > 0) deserialize_from_string(data, vd);
      alldatalock.lock();
      vertex_id_type localvid = create_vertex_if_missing(vid, atom2machine[owner], owner, data.size() > 0, vd);
      localstore.set_vertex_version(localvid, 1);
      alldatalock.unlock();
    } else if (command == 'd') {
      vertex_id_type src; vertex_id_type target; std::string data;
      uint16_t srcowner, targetowner;
      iarc >> src >> srcowner >> target >> targetowner >> data;
      
      EdgeData ed;
      if (data.size() > 0) deserialize_from_string(data, ed);
      alldatalock.lock();
      vertex_id_type localsrcvid =  create_vertex_if_missing(src, atom2machine[srcowner], srcowner, false);
      vertex_id_type localtargetvid = create_vertex_if_missing(target, atom2machine[targetowner], targetowner, false);
      std::pair<bool, edge_id_type> hasedge = localstore.find(localsrcvid, localtargetvid);
      edge_id_type eid;
      if (hasedge.first) eid = hasedge.second;
      else eid = localstore.add_edge(localsrcvid, localtargetvid);
      if (data.size() > 0)  {
        localstore.edge_data(eid) = ed;
        localstore.set_edge_version(eid, 1);
      }
      alldatalock.unlock();
    } else if (command == 'k') {
      vertex_id_type vid; vertex_color_type color;
      iarc >> vid >> color;
      vertex_id_type localvid = globalvid_to_localvid(vid);
      alldatalock.lock();
      localstore.color(localvid) = color;
      alldatalock.unlock();
    } else if (command == 'l') {
      // ignored
      vertex_id_type vid; uint16_t owner;
      iarc >> vid >> owner;
      // ignored
    }
  }
  in_file.close();
}
      

template <typename VertexData, typename EdgeData>
void distributed_graph<VertexData,EdgeData>::shuffle_local_vertices_to_start() {
  // ok. There are a whole bunch of stuff to update.
  // the entire localstore has to be shuffled
  // and in this class, local2globalvid and localvid2owner 
  // has to be shuffled
  
  std::vector<size_t> renumber(localvid2owner.size(), 0);

  for (size_t i = 0;i < localvid2owner.size(); ++i) renumber[i] = i;
  
    // advance cur through everything which is already in the right place
  {
    int j = (int)(localvid2owner.size()) - 1;
    int i = 0;
    while (j > i) {
      while(j > i && localvid2owner[i] == rmi.procid()) ++i;
      while(j > i && localvid2owner[j] != rmi.procid()) --j;
      if ( j <= i) break;
      renumber[i] = j;
      renumber[j] = i;
      ++i;
      --j;
    }
  }
  // for local2globalvid and localvid2owner, they are small 
  // and we just do the renumbering out of place
  outofplace_shuffle(local2globalvid, renumber);
  outofplace_shuffle(localvid2owner, renumber);
  outofplace_shuffle(localvid2atom, renumber);
  // check that invariance is good
  bool isowned = true;
  for (size_t i = 0;i < localvid2owner.size(); ++i) {
    if (isowned && localvid2owner[i] != rmi.procid()) {
      isowned = false;
    }
    if (isowned == false) {
      ASSERT_MSG(localvid2owner[i] != rmi.procid(), "local VID invariant not preserved");
    }
  }
  // shuffle global2localvid
  typename global2localvid_type::iterator iter = global2localvid.begin();
  while(iter != global2localvid.end()) {
    iter->second = renumber[iter->second];
    ++iter;
  }
  localstore.shuffle_vertex_ids(renumber);
  // renumber array should now be sorted
  for (size_t i = 0;i < renumber.size(); ++i) ASSERT_EQ(renumber[i], i);
  
}
#endif
