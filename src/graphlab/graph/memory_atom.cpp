#include <sstream>
#include <map>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/memory_atom.hpp>
#include <graphlab/logger/logger.hpp>

namespace graphlab {

memory_atom::memory_atom(std::string filename, uint16_t atomid):atomid(atomid),filename(filename) {
  std::ifstream fin(filename.c_str(), std::ios::binary);
  if (fin.good() && fin.is_open()) {
    iarchive iarc(fin);
    iarc >> vertices >> vidmap >> vid2owner_segment;
  }
}

void memory_atom::add_vertex(vertex_id_type vid, uint16_t owner) {
  mut.lock();
  boost::unordered_map<uint64_t, size_t>::const_iterator iter = vidmap.find(vid);
  if (iter == vidmap.end()) {
    vertices.push_back(vertex_entry(vid, owner));
    vidmap[vid] = vertices.size() - 1;
  }
  else {
    vertices[iter->second].owner = owner;
  }
  mut.unlock();
}

/**
 * \brief Inserts vertex 'vid' into the file without data.
 * If the vertex already exists, nothing will be done.
 * Returns true if vertex was added.
 */
bool memory_atom::add_vertex_skip(vertex_id_type vid, uint16_t owner) {
  bool ret = false;
  mut.lock();
  boost::unordered_map<uint64_t, size_t>::const_iterator iter = vidmap.find(vid);
  if (iter == vidmap.end()) {
    vertices.push_back(vertex_entry(vid, owner));
    vidmap[vid] = vertices.size() - 1;
    ret = true;
  }
  mut.unlock();
  return ret;
}


void memory_atom::add_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string &vdata) {
  mut.lock();
  boost::unordered_map<uint64_t, size_t>::const_iterator iter = vidmap.find(vid);
  if (iter == vidmap.end()) {
    vertices.push_back(vertex_entry(vid, owner, -1, vdata));
    vidmap[vid] = vertices.size() - 1;
  }
  else {
    vertices[iter->second].owner = owner;
    vertices[iter->second].vdata = vdata;
  }
  mut.unlock();
}
  
/**
 * \brief Inserts vertex 'vid' into the file without data.
 * If the vertex already exists, nothing will be done.
 * Returns true if vertex was added.
 */
void memory_atom::add_edge(vertex_id_type src, vertex_id_type target) {
  mut.lock();
  boost::unordered_map<uint64_t, size_t>::const_iterator srciter = vidmap.find(src);
  ASSERT_TRUE(srciter != vidmap.end());
  size_t vsrc = srciter->second;
  boost::unordered_map<uint64_t, size_t>::const_iterator destiter = vidmap.find(target);
  ASSERT_TRUE(destiter != vidmap.end());
  size_t vtarget = destiter->second;
  mut.unlock();
  
  edgemut[vsrc % 511].lock();
  vertices[vsrc].outedges.insert(std::make_pair(target, std::string()));
  edgemut[vsrc % 511].unlock();
  edgemut[vtarget % 511].lock();
  vertices[vtarget].inedges.insert(src);
  edgemut[vtarget % 511].unlock();
}


/**
 * \brief Inserts vertex 'vid' into the file without data.
 * If the vertex already exists, nothing will be done.
 * Returns true if vertex was added.
 */
bool memory_atom::add_edge_skip(vertex_id_type src, vertex_id_type target) {
  mut.lock();
  boost::unordered_map<uint64_t, size_t>::const_iterator srciter = vidmap.find(src);
  ASSERT_TRUE(srciter != vidmap.end());
  size_t vsrc = srciter->second;
  boost::unordered_map<uint64_t, size_t>::const_iterator destiter = vidmap.find(target);
  ASSERT_TRUE(destiter != vidmap.end());
  size_t vtarget = destiter->second;
  mut.unlock();
  
  edgemut[vsrc % 511].lock();
  if (vertices[vsrc].outedges.count(target) > 0) {
    edgemut[vsrc % 511].unlock();
    return false;
  }
  vertices[vsrc].outedges.insert(std::make_pair(target, std::string()));
  edgemut[vsrc % 511].unlock();
  edgemut[vtarget % 511].lock();
  vertices[vtarget].inedges.insert(src);
  edgemut[vtarget % 511].unlock();
  return true;
}


void memory_atom::add_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string &edata) {
  mut.lock();
  boost::unordered_map<uint64_t, size_t>::const_iterator srciter = vidmap.find(src);
  ASSERT_TRUE(srciter != vidmap.end());
  size_t vsrc = srciter->second;
  boost::unordered_map<uint64_t, size_t>::const_iterator destiter = vidmap.find(target);
  ASSERT_TRUE(destiter != vidmap.end());
  size_t vtarget = destiter->second;
  mut.unlock();
  
  edgemut[vsrc % 511].lock();
  vertices[vsrc].outedges[target] = edata;
  edgemut[vsrc % 511].unlock();
  edgemut[vtarget % 511].lock();
  vertices[vtarget].inedges.insert(src);
  edgemut[vtarget % 511].unlock();
}



void memory_atom::set_vertex(vertex_id_type vid, uint16_t owner) {
  mut.lock();
  boost::unordered_map<uint64_t, size_t>::iterator iter = vidmap.find(vid);
  ASSERT_TRUE(iter != vidmap.end());
  size_t v = iter->second;
  mut.unlock();
  vertices[v].owner = owner;
}


void memory_atom::set_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string &vdata) {

  mut.lock();
  boost::unordered_map<uint64_t, size_t>::const_iterator iter = vidmap.find(vid);
  ASSERT_TRUE(iter != vidmap.end());
  size_t v = iter->second;
  mut.unlock();
  vertices[v].owner = owner;
  vertices[v].vdata = vdata;
}



void memory_atom::set_edge(vertex_id_type src, vertex_id_type target) {
  add_edge(src, target);
}


void memory_atom::set_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string &edata) {
  add_edge_with_data(src, target, edata);
}





bool memory_atom::get_vertex(vertex_id_type vid, uint16_t &owner) {
  mut.lock();
  boost::unordered_map<uint64_t, size_t>::const_iterator iter = vidmap.find(vid);
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
  boost::unordered_map<uint64_t, size_t>::const_iterator iter = vidmap.find(vid);
  
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
  boost::unordered_map<uint64_t, size_t>::const_iterator srciter = vidmap.find(src);
  if(srciter == vidmap.end()) {
    mut.unlock();
    return false;
  }
  size_t vsrc = srciter->second;
  mut.unlock();
  std::map<vertex_id_type, std::string>::const_iterator iter = vertices[vsrc].outedges.find(target);
  if (iter == vertices[vsrc].outedges.end()) return false;
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
  for (size_t i = 0; i < ret.size(); ++i) {
    if (vertices[i].owner != atomid) {
      ++ret[vertices[i].owner];
    }
  }
  return ret;
}

std::vector<memory_atom::vertex_id_type> memory_atom::get_out_vertices(vertex_id_type vid) {
  std::vector<vertex_id_type> ret;
  boost::unordered_map<uint64_t, size_t>::const_iterator vidmap_iter = vidmap.find(vid);
  if (vidmap_iter == vidmap.end()) return ret;
  size_t v = vidmap_iter->second;
  
  std::map<vertex_id_type, std::string>::const_iterator iter = vertices[v].outedges.begin();
  ret.reserve(vertices[v].outedges.size());
  
  while(iter != vertices[v].outedges.end()) {
    ret.push_back(iter->first);
    ++iter;
  }
  return ret;
}


std::vector<memory_atom::vertex_id_type> memory_atom::get_in_vertices(vertex_id_type vid) {
  std::vector<vertex_id_type> ret;
  boost::unordered_map<uint64_t, size_t>::const_iterator vidmap_iter = vidmap.find(vid);
  if (vidmap_iter == vidmap.end()) return ret;
  size_t v = vidmap_iter->second;
  
  std::set<vertex_id_type>::const_iterator iter = vertices[v].inedges.begin();
  ret.reserve(vertices[v].inedges.size());
  
  while(iter != vertices[v].inedges.end()) {
    ret.push_back(*iter);
    ++iter;
  }
  return ret;
}


memory_atom::vertex_color_type memory_atom::get_color(vertex_id_type vid) {
  boost::unordered_map<uint64_t, size_t>::const_iterator viter = vidmap.find(vid);
  if (viter == vidmap.end()) return vertex_color_type(-1);
  
  return vertices[viter->second].color;
}


void memory_atom::set_color(vertex_id_type vid, vertex_color_type color) {
  std::vector<vertex_id_type> ret;
  boost::unordered_map<uint64_t, size_t>::const_iterator viter = vidmap.find(vid);
  if (viter == vidmap.end()) return;
  
  vertices[viter->second].color = color;

}

vertex_color_type memory_atom::max_color() {
  vertex_color_type m = 0;
  for (size_t i = 0;i < vertices.size(); ++i) {
    if (vertices[i].color > m) m = vertices[i].color;
  }
  return m;
}

uint16_t memory_atom::get_owner(vertex_id_type vid) {
  boost::unordered_map<uint64_t, uint16_t>::const_iterator iter = vid2owner_segment.find(vid);
  if (iter != vid2owner_segment.end()) return iter->second;
  else return uint16_t(-1);
}

void memory_atom::set_owner(vertex_id_type vid, uint16_t owner) {
  vid2owner_segment[vid] = owner;
}

void memory_atom::clear() {
  vid2owner_segment.clear();
  vertices.clear();
}

void memory_atom::synchronize() {
  std::ofstream fout(filename.c_str(), std::ios::binary);
  oarchive oarc(fout);
  oarc << vertices << vidmap << vid2owner_segment;
}
    
}
