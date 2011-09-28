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
  for (int i = 0;i < (int)(atoms_in_curpart.size()); ++i) {
    std::string fname = atomindex.atoms[atoms_in_curpart[i]].file;
    playback_dump(fname + ".dump", atoms_in_curpart[i], atom2machine, do_not_load_data);
  }
  
  // playback complete.
  // now to shuffle IDs around to rearrange all of the 
  // owned vertices to the start
  shuffle_local_vertices_to_start();
  construct_ghost_auxiliaries();
}

template <typename VertexData, typename EdgeData>
void distributed_graph<VertexData,EdgeData>::playback_dump(std::string filename,
                                      size_t atomid,
                                      std::vector<procid_t> atom2machine,
                                      bool do_not_load_data) {

  std::ifstream in_file(filename.c_str(), std::ios::binary);

  boost::iostreams::filtering_stream<boost::iostreams::input> fin; 
  fin.push(boost::iostreams::gzip_decompressor());
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
      // have we seen this vertex before?
      typename global2localvid_type::const_iterator iter = global2localvid.find(vid);
      // nope. Insert it.
      if (iter == global2localvid.end()) {
        vertex_id_t localvid = localstore.add_vertex(VertexData());
        global2localvid[vid] = localvid;
        local2globalvid.push_back(vid);
        localvid2owner.push_back(atom2machine[owner]);
        globalvid2owner.set(vid, rmi.procid());
        if (owner == atomid) ownedvertices.push_back(vid);
        else ghostvertices.push_back(vid);
      }
    } else if (command == 'c') {
      vertex_id_type vid; uint16_t owner; std::string data;
      iarc >> vid >> owner >> data;
      // deserialize it
      VertexData vd;
      deserialize_from_string(data, vd);

      // have we seen this vertex before?
      typename global2localvid_type::const_iterator iter = global2localvid.find(vid);
      // nope. Insert it.
      if (iter == global2localvid.end()) {
        vertex_id_t localvid = localstore.add_vertex(vd);
        global2localvid[vid] = localvid;
        local2globalvid.push_back(vid);
        localvid2owner.push_back(atom2machine[owner]);
        globalvid2owner.set(vid, rmi.procid());
        if (owner == atomid) ownedvertices.push_back(vid);
        else ghostvertices.push_back(vid);
      }
      else {
        localstore.vertex_data(iter->second) = vd;
      }
    } else if (command == 'd') {
      vertex_id_type src; vertex_id_type target; std::string data;
      uint16_t srcowner, targetowner;
      iarc >> src >> srcowner >> target >> targetowner >> data;
      
      EdgeData ed;
      deserialize_from_string(data, ed);

      vertex_id_t localsrcvid = globalvid_to_localvid(src);
      vertex_id_t localtargetvid = globalvid_to_localvid(target);
      edge_id_t eid = localstore.add_edge(localsrcvid, localtargetvid);
      localstore.edge_data(eid) = ed;
    } else if (command == 'k') {
      vertex_id_type vid; vertex_color_type color;
      iarc >> vid >> color;
      vertex_id_t localvid = globalvid_to_localvid(vid);
      localstore.color(localvid) = color;
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
  size_t cur = 0; 
  for (size_t i = 0;i < localvid2owner.size(); ++i) {
    if (localvid2owner[i] != rmi.procid()) {
      renumber[cur] = i;
      ++cur;
      // advance cur through everything which is already in the right place
      while(cur < localvid2owner.size() && 
            localvid2owner[cur] == rmi.procid()) {
        renumber[cur] = cur;
        ++cur;
      }
    }
  }
  // for local2globalvid and localvid2owner, they are small 
  // and we just do the renumbering out of place
  outofplace_shuffle(local2globalvid, renumber);
  outofplace_shuffle(localvid2owner, renumber);
  // check that invariance is good
  bool isowned = true;
  for (size_t i = 0;i < localvid2owner.size(); ++i) {
    if (isowned && localvid2owner[i] != rmi.procid()) {
      isowned = false;
    }
    else if (isowned == false) {
      ASSERT_MSG(localvid2owner[i] == rmi.procid(), "local VID invariant not preserved");
    }
  }
  
  localstore.shuffle_vertex_ids(renumber);
  // renumber array should now be sorted
  for (size_t i = 0;i < renumber.size(); ++i) ASSERT_EQ(renumber[i], i);
}
#endif