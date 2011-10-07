#ifndef FROM_DISTRIBUTED_GRAPH_INCLUDE

#warning "distributed_graph_incremental_loader.hpp should not be included directly."
#warning "You should include only distributed_graph.hpp"
#warning "I will fix this for you now, but don't do it again!"

#include <graphlab/distributed2/graph/distributed_graph.hpp>


#else
#include <graphlab/util/generics/shuffle.hpp>

inline void count_vertices_and_edges(std::string filename,
                                      std::vector<procid_t> atom2machine,
                                      procid_t mymachine,
                                      boost::unordered_set<vertex_id_t> &vertices,
                                      size_t &localedges) {
  localedges = 0;

  boost::unordered_set<vertex_id_t> ret;
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
      vertex_id_t vid; uint16_t owner;
      iarc >> vid >> owner;
      vertices.insert(vid);
    } else if (command == 'c') {
      vertex_id_t vid; uint16_t owner; std::string data;
      iarc >> vid >> owner >> data;
      vertices.insert(vid);
    } else if (command == 'd') {
      vertex_id_t src; vertex_id_t target; std::string data;
      uint16_t srcowner, targetowner;
      iarc >> src >> srcowner >> target >> targetowner >> data;
      // if there is data, it must be local
      if (data.size() > 0) localedges++;
      else {
        // if it is outedge to some other machine. I need to add it to
        // also, it will only appear once
        localedges += (atom2machine[targetowner] != mymachine);
      }
      vertices.insert(src);       
      vertices.insert(target);
    } else if (command == 'k') {
      vertex_id_t vid; vertex_color_type color;
      iarc >> vid >> color;
    } else if (command == 'l') {
      // ignored
      vertex_id_t vid; uint16_t owner;
      iarc >> vid >> owner;
      // ignored
    }
  }
  fin.pop();
  fin.pop();
  in_file.close();
}

template <typename VertexData, typename EdgeData>
typename distributed_graph<VertexData,EdgeData>::vertex_id_type
  distributed_graph<VertexData,EdgeData>::incremental_loader_add_vertex(vertex_id_type globalvid,
                                                                      procid_t machine,
                                                                      uint16_t sourceatom,
                                                                      bool overwritedata,
                                                                      const VertexData &vdata) {
    typename global2localvid_type::const_iterator iter = global2localvid.find(globalvid);
    ASSERT_TRUE(iter != global2localvid.end());
    if (overwritedata) {
      localstore.vertex_data(iter->second) = vdata;
    }
    localvid2atom[iter->second] = sourceatom;
    if (machine == rmi.procid()) {
      globalvid2owner.set(globalvid, rmi.procid());
    }
      return iter->second;
}

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
  
  logstream(LOG_INFO) << "First pass: Counting size of local store " << std::endl;
  
  std::vector<boost::unordered_set<vertex_id_type> > vertexset(omp_get_max_threads());
  size_t numedges = 0;
  
  #pragma omp parallel for reduction(+ : numedges)
  for (int i = 0;i < (int)(atoms_in_curpart.size()); ++i) {
    std::cout << ".";
    std::cout.flush();
    std::string fname = atomindex.atoms[atoms_in_curpart[i]].file;
    count_vertices_and_edges(fname + ".dump",
                              atom2machine,
                              rmi.procid(),
                              vertexset[omp_get_thread_num()],
                              numedges);
  }
  std::cout << std::endl;
  // create the vertex mapping
  for (size_t i = 0;i < vertexset.size() ; ++i) {
    std::copy(vertexset[i].begin(), vertexset[i].end(), 
              std::inserter(local2globalvid, local2globalvid.begin()));
  }
  for (size_t i = 0; i < local2globalvid.size(); ++i) global2localvid[local2globalvid[i]] = i;
  localvid2atom.resize(local2globalvid.size());
  
  logstream(LOG_INFO) <<  "Creating:" << local2globalvid.size() << " vertices," << numedges << " edges" << std::endl;
  // now lets construct the graph structure
  localstore.create_store(local2globalvid.size(), numedges);

  // create all the vertices

  // initiate playback
  logstream(LOG_INFO) << "Second pass: Loading data " << std::endl;
  std::vector<mutex> edgelockset;
  edgelockset.resize(1 << 14);
  atomic<edge_id_t> edgecount(0);
  #pragma omp parallel for
  for (int i = 0;i < (int)(atoms_in_curpart.size()); ++i) {
    std::cout << ".";
    std::cout.flush();
    std::string fname = atomindex.atoms[atoms_in_curpart[i]].file;
    playback_dump(fname + ".dump", 
                  atoms_in_curpart[i], 
                  atom2machine, 
                  rmi.procid(),
                  do_not_load_data, edgelockset, edgecount);
  }
  std::cout << std::endl;
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
void distributed_graph<VertexData,EdgeData>::playback_dump(std::string filename,
                                      size_t atomid,
                                      std::vector<procid_t> atom2machine,
                                      procid_t mymachine,
                                      bool do_not_load_data,
                                      std::vector<mutex>& edgelockset,
                                      atomic<edge_id_type>& edgecounter) {

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
      incremental_loader_add_vertex(vid, atom2machine[owner], owner);
    } else if (command == 'c') {
      vertex_id_type vid; uint16_t owner; std::string data;
      iarc >> vid >> owner >> data;
      // deserialize it
      VertexData vd;
      if (data.size() > 0) deserialize_from_string(data, vd);
      vertex_id_type localvid = incremental_loader_add_vertex(vid, atom2machine[owner], owner, data.size() > 0, vd);
      localstore.set_vertex_version(localvid, 1);
    } else if (command == 'd') {
      vertex_id_type src; vertex_id_type target; std::string data;
      uint16_t srcowner, targetowner;
      iarc >> src >> srcowner >> target >> targetowner >> data;
      
      EdgeData ed;
      if (data.size() > 0) deserialize_from_string(data, ed);
      vertex_id_type localsrcvid =  incremental_loader_add_vertex(src, atom2machine[srcowner], srcowner, false);
      vertex_id_type localtargetvid = incremental_loader_add_vertex(target, atom2machine[targetowner], targetowner, false);
      size_t locka = localsrcvid % edgelockset.size();
      size_t lockb = localtargetvid % edgelockset.size();
      if (data.size() > 0 || atom2machine[targetowner] != mymachine) {
        edge_id_t eid = edgecounter.inc_ret_last();
        edgelockset[std::min(locka, lockb)].lock();
        if (locka != lockb)  edgelockset[std::max(locka, lockb)].lock();
        localstore.add_edge(eid, localsrcvid, localtargetvid);
        if (locka != lockb)  edgelockset[std::max(locka, lockb)].unlock();
        edgelockset[std::min(locka, lockb)].unlock();
        if (data.size() > 0)  {
          localstore.edge_data(eid) = ed;
          localstore.set_edge_version(eid, 1);
        }
      }
    } else if (command == 'k') {
      vertex_id_type vid; vertex_color_type color;
      iarc >> vid >> color;
      vertex_id_type localvid = globalvid_to_localvid(vid);
      localstore.color(localvid) = color;
    } else if (command == 'l') {
      // ignored
      vertex_id_type vid; uint16_t owner;
      iarc >> vid >> owner;
      // ignored
    }
  }
  fin.pop();
  fin.pop();
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
