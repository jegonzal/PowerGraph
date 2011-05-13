/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <sstream>
#include <map>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/graph/disk_atom.hpp>
#include <kchashdb.h>


namespace graphlab {

disk_atom::disk_atom(std::string filename, 
                            uint16_t atomid):atomid(atomid),filename(filename) {
  db.tune_options(kyotocabinet::HashDB::TSMALL);
  db.tune_buckets(100 * 1000);
#if __LP64__
  db.tune_map(256 * 1024 * 1024); // 256MB
#endif
  cache_invalid = true;
  ASSERT_TRUE(db.open(filename));
  // get the pointers to the linked list of vertices
  if (db.get("head_vid", 7, (char*)&head_vid, sizeof(head_vid)) == -1) head_vid = (uint64_t)(-1);
  if (db.get("tail_vid", 7, (char*)&tail_vid, sizeof(tail_vid)) == -1) tail_vid = (uint64_t)(-1);
  if (db.get("numv", 4, (char*)&numv.value, sizeof(numv.value)) == -1) numv = 0;
  if (db.get("nume", 4, (char*)&nume.value, sizeof(nume.value)) == -1) nume = 0;
  if (db.get("numlocalv", 9, (char*)&numlocalv.value, sizeof(numlocalv.value)) == -1) numlocalv = 0;
  if (db.get("numlocale", 9, (char*)&numlocale.value, sizeof(numlocale.value)) == -1) numlocale = 0;
}

void disk_atom::precache() {
  kyotocabinet::HashDB::Cursor* cur = db.cursor();
  cur->jump();
  std::string key, val;
  while (cur->get(&key, &val, true)) {
    cache.set(key, val);
  }
  delete cur;
  cache_invalid = false;
}

disk_atom::~disk_atom() {
  synchronize();
  db.close();
}

void disk_atom::synchronize() {
  // update the linked list
  db.set("head_vid", 7, (char*)&head_vid, sizeof(head_vid));
  db.set("tail_vid", 7, (char*)&tail_vid, sizeof(tail_vid));
  db.set("numv", 4, (char*)&numv.value, sizeof(numv.value));
  db.set("nume", 4, (char*)&nume.value, sizeof(nume.value));
  db.set("numlocalv", 9, (char*)&numlocalv.value, sizeof(numlocalv.value));
  db.set("numlocale", 9, (char*)&numlocale.value, sizeof(numlocale.value));
  db.synchronize();
}

void disk_atom::add_vertex(vertex_id_t vid, uint16_t owner) {
  if (!add_vertex_skip(vid, owner)) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner;
    strm.flush();
    db.set("v"+id_to_str(vid), strm.str());
    cache_invalid = true;
  }
}


bool disk_atom::add_vertex_skip(vertex_id_t vid, uint16_t owner) {
  std::stringstream strm;
  oarchive oarc(strm);    
  oarc << owner;
  strm.flush();
  if (db.add("v"+id_to_str(vid), strm.str())) {
    // first entry in linked list
    mut.lock();
    if (head_vid == (uint64_t)(-1)) {
      head_vid = vid;
      tail_vid = vid;
      mut.unlock();
    }
    else {
      // update linked list
      std::string tail_next_key = "ll" + id_to_str(tail_vid);
      tail_vid = vid;
      mut.unlock();
      db.set(tail_next_key.c_str(), tail_next_key.length(), 
             (char*)&vid, sizeof(vid));
      cache_invalid = true;
    }
    numv.inc();
    if (owner == atomid) numlocalv.inc();
    return true;
  }
  return false;
}


void disk_atom::add_edge(vertex_id_t src, vertex_id_t target) {
  if (!add_edge_skip(src, target)) {
    db.set("e"+id_to_str(src)+"_"+id_to_str(target), std::string(""));
    cache_invalid = true;
  }
}


bool disk_atom::add_edge_skip(vertex_id_t src, vertex_id_t target) {
  if (db.add("e"+id_to_str(src)+"_"+id_to_str(target), std::string(""))) {
    // increment the number of edges
    nume.inc();
    // append to the adjacency entries
    std::string oadj_key = "o"+id_to_str(src);
    uint64_t target64 = (uint64_t)target;
    db.append(oadj_key.c_str(), oadj_key.length(), (char*)&target64, sizeof(target64));
    
    std::string iadj_key = "i"+id_to_str(target);
    uint64_t src64 = (uint64_t)src;
    db.append(iadj_key.c_str(), iadj_key.length(), (char*)&src64, sizeof(src64));
    return true;
  }
  return false;
}

std::vector<vertex_id_t> disk_atom::enumerate_vertices() {
  std::vector<vertex_id_t> ret;
  if (head_vid == (uint64_t)(-1)) return ret;
  else {
    uint64_t curvid = head_vid;
    while(1) {
      ret.push_back((vertex_id_t)(curvid));
      std::string next_key = "ll" + id_to_str(curvid);
      if (cache_invalid || 
          cache.get(next_key.c_str(), next_key.length(), (char*)&curvid, sizeof(curvid)) == -1) {
        if (db.get(next_key.c_str(), next_key.length(), (char*)&curvid, sizeof(curvid)) == -1) {
          break;
        }
      }
    }
  }
  return ret;
}


bool disk_atom::get_vertex(vertex_id_t vid, uint16_t &owner) {
  std::string val;
  std::string key = "v"+id_to_str(vid);
  if (cache_invalid || cache.get(key, &val) == false) {
    if (db.get("v"+id_to_str(vid), &val) == false) return false;
  }
  
  boost::iostreams::stream<boost::iostreams::array_source> 
                              istrm(val.c_str(), val.length());   
  iarchive iarc(istrm);
  iarc >> owner;
  return true;
}


std::map<uint16_t, uint32_t> disk_atom::enumerate_adjacent_atoms() {
  std::map<uint16_t, uint32_t> ret;
  if (head_vid == (uint64_t)(-1)) return ret;
  else {
    uint64_t curvid = head_vid;
    while(1) {
      uint16_t owner;
      get_vertex((vertex_id_t)curvid, owner);
      if (owner != atomid) ret[owner]++;
      std::string next_key = "ll" + id_to_str(curvid);
      if (db.get(next_key.c_str(), next_key.length(), (char*)&curvid, sizeof(curvid)) == -1) {
        break;
      }
    }
  }
  return ret;
}

uint32_t disk_atom::max_color() {
  uint32_t mcolor = 0;
  if (head_vid == (uint64_t)(-1)) return mcolor;
  else {
    uint64_t curvid = head_vid;
    while(1) {

      uint32_t c = get_color(curvid);
      if (c != (uint32_t)(-1)) {
        mcolor = std::max(mcolor, c);
      }
      std::string next_key = "ll" + id_to_str(curvid);
      if (db.get(next_key.c_str(), next_key.length(), (char*)&curvid, sizeof(curvid)) == -1) {
        break;
      }
    }
  }
  return mcolor;
}


std::vector<vertex_id_t> disk_atom::get_in_vertices(vertex_id_t vid) {
  std::vector<vertex_id_t> ret;
  std::string val;
  std::string key = "i"+id_to_str(vid);
  if ((cache_invalid == false && cache.get(key, &val)) || db.get(key, &val)) {
    const uint64_t* v = reinterpret_cast<const uint64_t*>(val.c_str());
    ASSERT_TRUE(val.length() % 8 == 0);
    size_t numel = val.length() / 8;
    ret.resize(numel);
    for (size_t i = 0;i < numel; ++i) ret[i] = v[i];
  }
  return ret;    
}
   
   

std::vector<vertex_id_t> disk_atom::get_out_vertices(vertex_id_t vid) {
  std::vector<vertex_id_t> ret;
  std::string val;
  std::string key = "o"+id_to_str(vid);
  if ((cache_invalid == false && cache.get(key, &val)) || db.get(key, &val)) {
    const uint64_t* v = reinterpret_cast<const uint64_t*>(val.c_str());
    size_t numel = val.length() / 8;
    ASSERT_TRUE(val.length() % 8 == 0);
    ret.resize(numel);
    for (size_t i = 0;i < numel; ++i) ret[i] = v[i];
  }
  return ret;    
}



uint32_t disk_atom::get_color(vertex_id_t vid) {
  std::string key = "c" + id_to_str(vid);
  uint32_t ret;
  if (cache_invalid == false && 
      cache.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) != -1) return ret;

  if (db.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) == -1) ret = (uint32_t)(-1);
  return ret;
}


void disk_atom::set_color(vertex_id_t vid, uint32_t color) {
  std::string key = "c" + id_to_str(vid);
  db.set(key.c_str(), key.length(), (char*)&color, sizeof(color));
  cache_invalid = true;
}


uint16_t disk_atom::get_owner(vertex_id_t vid) {
  std::string key = "h" + id_to_str(vid);
  uint16_t ret;
  if (cache_invalid == false && 
      cache.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) != -1) return ret;
  
  if (db.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) == -1) ret = (uint16_t)(-1);
  return ret;
}


void disk_atom::set_owner(vertex_id_t vid, uint16_t owner) {
  std::string key = "h" + id_to_str(vid);
  db.set(key.c_str(), key.length(), (char*)&owner, sizeof(owner));
  cache_invalid = true;
}


void disk_atom::clear() {
  head_vid = (uint64_t)(-1);   
  tail_vid = (uint64_t)(-1);
  numv.value = 0;
  nume.value = 0;
  numlocalv.value = 0;
  numlocale.value = 0;
  db.clear();
}


} // namespace graphlab

