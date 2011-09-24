#include <sstream>
#include <map>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/memory_atom.hpp>
#include <graphlab/logger/logger.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace graphlab {

memory_atom::memory_atom(std::string filename, uint16_t atomid):atomid(atomid),filename(filename) {
  std::ifstream in_file(filename.c_str(), std::ios::binary);

  if (in_file.good() && in_file.is_open()) {
    boost::iostreams::filtering_stream<boost::iostreams::input> fin; 
    fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);

    iarchive iarc(fin);
    uint64_t nv,ne;
    iarc >> nv >> ne
         >> vertices >> vidmap >> vid2owner_segment;
    numv.value = nv;
    nume.value = ne;
    mutated = false;
  }
  else {
    mutated = true;
  }
}


void memory_atom::add_vertex(vertex_id_type vid, uint16_t owner) {
  if (!add_vertex_skip(vid, owner)) {
    mut.lock();
    boost::unordered_map<vertex_id_t, size_t>::const_iterator iter = vidmap.find(vid);
    vertices[iter->second].owner = owner;
    mutated = true;
    mut.unlock();
  }
}

/**
 * \brief Inserts vertex 'vid' into the file without data.
 * If the vertex already exists, nothing will be done.
 * Returns true if vertex was added.
 */
bool memory_atom::add_vertex_skip(vertex_id_type vid, uint16_t owner) {
  mutated = true;
  bool ret = false;
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator iter = vidmap.find(vid);
  if (iter == vidmap.end()) {
    if (owner == atom_id()) numv.inc();
    vertices.push_back(vertex_entry(vid, owner));
    vidmap[vid] = vertices.size() - 1;
    ret = true;
  }
  mut.unlock();
  return ret;
}


void memory_atom::add_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string &vdata) {
  mutated = true;
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator iter = vidmap.find(vid);
  if (iter == vidmap.end()) {
    if (owner == atom_id()) numv.inc();
    vertices.push_back(vertex_entry(vid, owner, -1, vdata));
    vidmap[vid] = vertices.size() - 1;
  }
  else {
    vertices[iter->second].owner = owner;
    vertices[iter->second].vdata = vdata;
  }
  mut.unlock();
}


void memory_atom::add_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string &edata) {
  mutated = true;
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator srciter = vidmap.find(src);

  if (srciter == vidmap.end()) {
    mut.unlock();
    if (edata.size() > 0) add_vertex_skip(src, uint16_t(-1));     // source may or may not be a ghost
    else add_vertex_skip(src, atom_id()); // source is definitely local
    mut.lock();
    srciter = vidmap.find(src);
  }
  size_t vsrc = srciter->second;
  
  boost::unordered_map<vertex_id_t, size_t>::const_iterator destiter = vidmap.find(target);
  if (destiter == vidmap.end()) {
    mut.unlock();
    if (edata.size() > 0) add_vertex_skip(target, atom_id());     // target is definitely local
    else add_vertex_skip(target, uint16_t(-1)); // target is definitely a ghost
    mut.lock();
    destiter = vidmap.find(target);
  }
  size_t vtarget = destiter->second;
  
  if (vertices[vsrc].outedges.count(target) == 0) {
    if (edata.length() > 0) nume.inc();
    vertices[vsrc].outedges.insert(target);
  }
  vertices[vtarget].inedges.insert(std::make_pair(src, edata));
  mut.unlock();

}

void memory_atom::add_edge_with_data(vertex_id_type src, uint16_t srcowner,
                                     vertex_id_type target, uint16_t targetowner, 
                                     const std::string &edata) {
  mutated = true;
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator srciter = vidmap.find(src);

  if (srciter == vidmap.end()) {
    mut.unlock();
    add_vertex_skip(src, srcowner);
    mut.lock();
    srciter = vidmap.find(src);
  }
  size_t vsrc = srciter->second;
  
  boost::unordered_map<vertex_id_t, size_t>::const_iterator destiter = vidmap.find(target);
  if (destiter == vidmap.end()) {
    mut.unlock();
    add_vertex_skip(target, targetowner);
    mut.lock();
    destiter = vidmap.find(target);
  }
  size_t vtarget = destiter->second;

  if (vertices[vsrc].outedges.count(target) == 0) {
    if (edata.length() > 0) nume.inc();
    vertices[vsrc].outedges.insert(target);
  }
  vertices[vtarget].inedges.insert(std::make_pair(src, edata));
  mut.unlock();
}



void memory_atom::set_vertex(vertex_id_type vid, uint16_t owner) {
  mutated = true;
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator iter = vidmap.find(vid);
  ASSERT_TRUE(iter != vidmap.end());
  size_t v = iter->second;
  mut.unlock();
  vertices[v].owner = owner;
}


void memory_atom::set_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string &vdata) {
  mutated = true;
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator iter = vidmap.find(vid);
  ASSERT_TRUE(iter != vidmap.end());
  size_t v = iter->second;
  mut.unlock();
  vertices[v].owner = owner;
  vertices[v].vdata = vdata;
}




void memory_atom::set_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string &edata) {
  mutated = true;
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator srciter = vidmap.find(src);
  ASSERT_TRUE(srciter != vidmap.end());
  size_t vsrc = srciter->second;
  
  boost::unordered_map<vertex_id_t, size_t>::const_iterator destiter = vidmap.find(target);
  ASSERT_TRUE(destiter != vidmap.end());
  size_t vtarget = destiter->second;

  ASSERT_EQ(vertices[vsrc].outedges.count(target), 1);
  vertices[vtarget].inedges[src] = edata;
  mut.unlock();

}





bool memory_atom::get_vertex(vertex_id_type vid, uint16_t &owner) {
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator iter = vidmap.find(vid);
  if (iter == vidmap.end()) {
    mut.unlock();
    return false;
  }
  size_t v = iter->second;
  mut.unlock();
  owner = vertices[v].owner;
  return true;
}



bool memory_atom::get_vertex_data(vertex_id_type vid, uint16_t &owner, std::string &vdata) {
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator iter = vidmap.find(vid);
  
  if(iter == vidmap.end()) {
    mut.unlock();
    return false;
  }
  size_t v = iter->second;
  mut.unlock();
  
  owner = vertices[v].owner;
  vdata = vertices[v].vdata;
  return true;
}

/**
 * \brief Reads a edge from the file returning results in 'owner' and 'vdata'.
 * Returns true if edge exists and false otherwise.
 * If there is no edge data stored, edata will not be modified.
 */
bool memory_atom::get_edge_data(vertex_id_type src, vertex_id_type target, std::string &edata) {
  mut.lock();
  boost::unordered_map<vertex_id_t, size_t>::const_iterator targetiter = vidmap.find(target);
  if(targetiter == vidmap.end()) {
    mut.unlock();
    return false;
  }
  size_t vtarget = targetiter->second;
  mut.unlock();
  std::map<vertex_id_type, std::string>::const_iterator iter = vertices[vtarget].inedges.find(src);
  if (iter == vertices[vtarget].inedges.end()) return false;
  edata = iter->second;
  return true;
}


std::vector<memory_atom::vertex_id_type> memory_atom::enumerate_vertices() {
  std::vector<vertex_id_type> ret;
  ret.resize(vertices.size());
  for (size_t i = 0; i < ret.size(); ++i) {
    ret[i] = vertices[i].vid;
  }
  return ret;
}


std::map<uint16_t, uint32_t> memory_atom::enumerate_adjacent_atoms() {
  std::map<uint16_t, uint32_t> ret;
  for (size_t i = 0; i < vertices.size(); ++i) {
    if (vertices[i].owner != atomid) {
      ++ret[vertices[i].owner];
    }
  }
  return ret;
}

std::vector<memory_atom::vertex_id_type> memory_atom::get_out_vertices(vertex_id_type vid) {
  std::vector<vertex_id_type> ret;
  boost::unordered_map<vertex_id_t, size_t>::const_iterator vidmap_iter = vidmap.find(vid);
  if (vidmap_iter == vidmap.end()) return ret;
  size_t v = vidmap_iter->second;
  
  std::copy(vertices[v].outedges.begin(), 
            vertices[v].outedges.end(), 
            std::inserter(ret, ret.begin()));
  return ret;
}


std::vector<memory_atom::vertex_id_type> memory_atom::get_in_vertices(vertex_id_type vid) {
  std::vector<vertex_id_type> ret;
  boost::unordered_map<vertex_id_t, size_t>::const_iterator vidmap_iter = vidmap.find(vid);
  if (vidmap_iter == vidmap.end()) return ret;
  size_t v = vidmap_iter->second;
  
  std::map<vertex_id_type, std::string>::const_iterator iter = vertices[v].inedges.begin();
  ret.reserve(vertices[v].inedges.size());
  
  while(iter != vertices[v].inedges.end()) {
    ret.push_back(iter->first);
    ++iter;
  }
  return ret;
}


memory_atom::vertex_color_type memory_atom::get_color(vertex_id_type vid) {
  mut.lock();
  vertex_color_type ret = vertex_color_type(-1);
  boost::unordered_map<vertex_id_t, size_t>::const_iterator viter = vidmap.find(vid);
  if (viter != vidmap.end()) {
    ret = vertices[viter->second].color;
  }
  mut.unlock();
  return ret;
}


void memory_atom::set_color(vertex_id_type vid, vertex_color_type color) {
  mutated = true;
  mut.lock();
  std::vector<vertex_id_type> ret;
  boost::unordered_map<vertex_id_t, size_t>::const_iterator viter = vidmap.find(vid);
  if (viter != vidmap.end()) {
    vertices[viter->second].color = color;
  }
  mut.unlock();

}

vertex_color_type memory_atom::max_color() {
  vertex_color_type m = 0;
  for (size_t i = 0;i < vertices.size(); ++i) {
    if (vertices[i].color > m) m = vertices[i].color;
  }
  return m;
}

uint16_t memory_atom::get_owner(vertex_id_type vid) {
  maplock.lock();
  boost::unordered_map<vertex_id_t, uint16_t>::const_iterator iter = vid2owner_segment.find(vid);
  uint16_t ret = 0;
  if (iter != vid2owner_segment.end()) ret = iter->second;
  else ret = uint16_t(-1);
  maplock.unlock();
  return ret;
}

void memory_atom::set_owner(vertex_id_type vid, uint16_t owner) {
  mutated = true;
  maplock.lock();
  vid2owner_segment[vid] = owner;
  maplock.unlock();
}

void memory_atom::clear() {
  vid2owner_segment.clear();
  vertices.clear();
  vidmap.clear();
}

void memory_atom::synchronize() {
  if (mutated) {
    std::ofstream out_file(filename.c_str(), std::ios::binary);
    boost::iostreams::filtering_stream<boost::iostreams::output> fout; 
    fout.push(boost::iostreams::gzip_compressor(boost::iostreams::zlib::best_compression));
    fout.push(out_file);

    oarchive oarc(fout);
    uint64_t nv,ne;
    nv = numv.value;
    ne = nume.value;
  
    oarc << nv << ne
         << vertices << vidmap << vid2owner_segment;
  }
}
    
}
