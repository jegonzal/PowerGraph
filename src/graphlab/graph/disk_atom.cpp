/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */

/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */



#include <sstream>
#include <map>
#include <ios>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/graph/memory_atom.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/graph/disk_atom.hpp>
#include <kchashdb.h>


namespace graphlab {

  disk_atom::disk_atom(std::string filename, 
                       uint16_t atomid):
                          atomid(atomid),filename(filename) {
    open_db();
  }

  void disk_atom::open_db() {
    db.tune_options(storage_type::TLINEAR);
    db.tune_buckets(1000);
    db.tune_page(32768);
#if __LP64__
    db.tune_map(256 * 1024 * 1024); // 256MB
#endif
    ASSERT_TRUE(db.open(filename));
    // get the pointers to the linked list of vertices
    if (db.get("numv", 4, (char*)&numv.value, sizeof(numv.value)) == -1) numv = 0;
    if (db.get("nume", 4, (char*)&nume.value, sizeof(nume.value)) == -1) nume = 0;
  }
  

  disk_atom::~disk_atom() {
    synchronize();
    db.close();
  }

  void disk_atom::synchronize() {
    // update the linked list
    db.set("numv", 4, (char*)&numv.value, sizeof(numv.value));
    db.set("nume", 4, (char*)&nume.value, sizeof(nume.value));
    db.synchronize();
  }

  void disk_atom::add_vertex(disk_atom::vertex_id_type vid, uint16_t owner) {
    if (!add_vertex_skip(vid, owner)) {
      std::stringstream strm;
      oarchive oarc(strm);    
      oarc << owner;
      strm.flush();
      db.set("v"+id_to_str(vid), strm.str());
    }
  }


  bool disk_atom::add_vertex_skip(disk_atom::vertex_id_type vid, uint16_t owner) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner;
    strm.flush();
    if (db.add("v"+id_to_str(vid), strm.str())) {
      uint64_t v64 = (uint64_t)vid;
      db.append("_vidlist", 8, (char*)&v64, sizeof(v64));
      if (owner == atom_id()) numv.inc();
      return true;
    }
    return false;
  }

  void disk_atom::add_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string &vdata) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner << vdata;
    strm.flush();
    if (db.add("v"+id_to_str(vid), strm.str())) {
      uint64_t v64 = (uint64_t)vid;
      db.append("_vidlist", 8, (char*)&v64, sizeof(v64));
      if (owner == atom_id()) numv.inc();
    }
    else {
      db.set("v"+id_to_str(vid), strm.str());
    }
  }

  void disk_atom::add_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string& edata) {
    if (db.add("e"+id_to_str(src)+"_"+id_to_str(target), edata)) {
      // increment the number of edges
      
      // target is definitely local. source may or may not be a ghost
      if (edata.size() > 0) {
        nume.inc();
        add_vertex_skip(target, atom_id());
        add_vertex_skip(src, (uint16_t)(-1));
      }
      else {
        // src is local target is definitely a ghost
        add_vertex_skip(src, atom_id());
        add_vertex_skip(target, (uint16_t)(-1));
      }
      
      // append to the adjacency entries
      std::string oadj_key = "o"+id_to_str(src);
      uint64_t target64 = (uint64_t)target;
      db.append(oadj_key.c_str(), oadj_key.length(), (char*)&target64, sizeof(target64));
    
      std::string iadj_key = "i"+id_to_str(target);
      uint64_t src64 = (uint64_t)src;
      db.append(iadj_key.c_str(), iadj_key.length(), (char*)&src64, sizeof(src64));
    }
    else {
      db.set("e"+id_to_str(src)+"_"+id_to_str(target), edata);
    }
  }


  void disk_atom::add_edge_with_data(vertex_id_type src, uint16_t srcowner,
                                     vertex_id_type target, uint16_t targetowner, const std::string& edata) {
    if (db.add("e"+id_to_str(src)+"_"+id_to_str(target), edata)) {
      // increment the number of edges
      
      // target is definitely local. source may or may not be a ghost
      if (edata.size() > 0) {
        nume.inc();
      }
      add_vertex_skip(src, srcowner);
      add_vertex_skip(target, targetowner);

      // append to the adjacency entries
      std::string oadj_key = "o"+id_to_str(src);
      uint64_t target64 = (uint64_t)target;
      db.append(oadj_key.c_str(), oadj_key.length(), (char*)&target64, sizeof(target64));
    
      std::string iadj_key = "i"+id_to_str(target);
      uint64_t src64 = (uint64_t)src;
      db.append(iadj_key.c_str(), iadj_key.length(), (char*)&src64, sizeof(src64));
    }
    else {
      db.set("e"+id_to_str(src)+"_"+id_to_str(target), edata);
    }
  }

  std::vector<disk_atom::vertex_id_type> disk_atom::enumerate_vertices() {
    std::vector<disk_atom::vertex_id_type> ret;
    // read the entire vertex list
    std::string vidlist;

    if (db.get(std::string("_vidlist"), &vidlist) == false) {
      return ret;
    }
    
    ASSERT_EQ(vidlist.size() % sizeof(uint64_t), 0);
    ret.resize(vidlist.size() / sizeof(uint64_t));
    const uint64_t* arr = reinterpret_cast<const uint64_t*>(vidlist.c_str());
    std::copy(arr, arr + ret.size(), ret.begin());
    return ret;
  }


  void disk_atom::set_vertex(vertex_id_type vid, uint16_t owner) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner;
    strm.flush();
    db.set("v"+id_to_str(vid), strm.str());
  }


  void disk_atom::set_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string &vdata) {
      std::stringstream strm;
      oarchive oarc(strm);    
      oarc << owner << vdata;
      strm.flush();
      db.set("v"+id_to_str(vid), strm.str());
    }
  

  void disk_atom::set_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string &edata) {
    std::string s;
    db.get("e"+id_to_str(src)+"_"+id_to_str(target), &s);
    db.set("e"+id_to_str(src)+"_"+id_to_str(target), edata);
  }

  bool disk_atom::get_vertex(disk_atom::vertex_id_type vid, uint16_t &owner) {
    std::string val;
    std::string key = "v"+id_to_str(vid);
    
    if (db.get(key, &val) == false) {
      return false;
    }

    boost::iostreams::stream<boost::iostreams::array_source> 
      istrm(val.c_str(), val.length());   
    iarchive iarc(istrm);
    iarc >> owner;
    return true;
  }

  bool disk_atom::get_vertex_data(vertex_id_type vid, uint16_t &owner, std::string &vdata) {
    std::string val;
    std::string key = "v"+id_to_str(vid);
    if (db.get(key, &val) == false) {
      return false; 
    }
    
    boost::iostreams::stream<boost::iostreams::array_source> 
            istrm(val.c_str(), val.length());   
    iarchive iarc(istrm);
    iarc >> owner;
    // try to deserialize vdata if exists
    istrm.peek();
    if (!istrm.eof()) iarc >> vdata;
    return true;
  }

  bool disk_atom::get_edge_data(vertex_id_type src, vertex_id_type target, std::string &edata) {
    edata = "";
    std::string key = "e"+id_to_str(src)+"_"+id_to_str(target);
    if (db.get(key, &edata) == false) {
      return false; 
    }
    return true;
  }



  std::map<uint16_t, uint32_t> disk_atom::enumerate_adjacent_atoms() {
    std::vector<disk_atom::vertex_id_type> vids = enumerate_vertices();
    std::map<uint16_t, uint32_t> ret;
    for (size_t i = 0;i < vids.size(); ++i) {
      uint16_t owner;
      get_vertex(vids[i], owner);
      if (owner != atomid) ret[owner]++;
    }
    return ret;
  }

  disk_atom::vertex_color_type disk_atom::max_color() {
    disk_atom::vertex_color_type mcolor = 0;
    
    std::vector<disk_atom::vertex_id_type> vids = enumerate_vertices();
    std::map<uint16_t, uint32_t> ret;
    for (size_t i = 0;i < vids.size(); ++i) {
      disk_atom::vertex_color_type c = get_color(vids[i]);
      if (c != disk_atom::vertex_color_type(-1)) {
        mcolor = std::max(mcolor, c);
      }
    }
    return mcolor;
  }


  std::vector<disk_atom::vertex_id_type> disk_atom::get_in_vertices(disk_atom::vertex_id_type vid) {
    std::vector<disk_atom::vertex_id_type> ret;
    std::string val;
    std::string key = "i"+id_to_str(vid);
    
    if (db.get(key, &val) == false) {
      return ret;
    }
    
    const uint64_t* v = reinterpret_cast<const uint64_t*>(val.c_str());
    ASSERT_TRUE(val.length() % 8 == 0);
    size_t numel = val.length() / 8;
    ret.resize(numel);
    for (size_t i = 0;i < numel; ++i) ret[i] = v[i];

    return ret;    
  }
   
   

  std::vector<disk_atom::vertex_id_type> disk_atom::get_out_vertices(disk_atom::vertex_id_type vid) {
    std::vector<disk_atom::vertex_id_type> ret;
    std::string val;
    std::string key = "o"+id_to_str(vid);
    
    if (db.get(key, &val) == false) {
      return ret;
    }

    const uint64_t* v = reinterpret_cast<const uint64_t*>(val.c_str());
    size_t numel = val.length() / 8;
    ASSERT_TRUE(val.length() % 8 == 0);
    ret.resize(numel);
    for (size_t i = 0;i < numel; ++i) ret[i] = v[i];

    return ret;    
  }



  disk_atom::vertex_color_type 
  disk_atom::get_color(disk_atom::vertex_id_type vid) {
    std::string key = "c" + id_to_str(vid);
    disk_atom::vertex_color_type  ret;
    if (db.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) == -1)  return disk_atom::vertex_color_type(-1);
    return ret;
  }


  void disk_atom::set_color(disk_atom::vertex_id_type vid, 
                            disk_atom::vertex_color_type color) {
    std::string key = "c" + id_to_str(vid);
    db.set(key.c_str(), key.length(), (char*)&color, sizeof(color));
  }


  uint16_t disk_atom::get_owner(disk_atom::vertex_id_type vid) {
    std::string key = "h" + id_to_str(vid);
    uint16_t ret;
    if (db.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) == -1) return (uint16_t)(-1); 
    
    return ret;
  }


  void disk_atom::set_owner(disk_atom::vertex_id_type vid, uint16_t owner) {
    std::string key = "h" + id_to_str(vid);
    db.set(key.c_str(), key.length(), (char*)&owner, sizeof(owner));
  }


  void disk_atom::clear() {
    numv.value = 0;
    nume.value = 0;
    db.clear();
  }


  disk_atom::vertex_id_type disk_atom::vertex_key_to_id(std::string s) {
    vertex_id_type vid;
    decompress_int<vertex_id_type>(s.c_str() + 1, vid);
    return vid;
  }


  void disk_atom::build_memory_atom(std::string fname) {
    memory_atom matom(fname, atomid);
    std::vector<vertex_id_type> vertices = disk_atom::enumerate_vertices();
    // add all the vertices
    for (size_t i = 0; i < vertices.size(); ++i) {
      uint16_t owner;
      std::string vdata;
      get_vertex_data(vertices[i], owner, vdata);
      if (vdata.length() == 0) matom.add_vertex(vertices[i], owner);
      else matom.add_vertex_with_data(vertices[i], owner, vdata);
      
      vertex_color_type c = get_color(vertices[i]);
      if (c != vertex_color_type(-1)) matom.set_color(vertices[i], c);
    }
    
    // add all the edges
    for (size_t i = 0; i < vertices.size(); ++i) {
      std::vector<vertex_id_type> inv = disk_atom::get_in_vertices(vertices[i]); 
      for (size_t j = 0;j < inv.size(); ++j) {
        std::string edata;
        disk_atom::get_edge_data(inv[j], vertices[i], edata);
        matom.add_edge_with_data(inv[j], vertices[i], edata);
      }
    }
    ASSERT_EQ(matom.num_edges(), num_edges());
    ASSERT_EQ(matom.num_vertices(), num_vertices());
    
    // finally, add all the owner hashes
    // open a cursor
    std::string key, val;
    storage_type::Cursor* cur = get_db().cursor();
    cur->jump();

    // begin iteration
    while (cur->get_key(&key, false)) {
      if (key[0] == 'h') {
        vertex_id_type vid = vertex_key_to_id(key);       
        cur->get_value(&val, true);
        uint16_t owner= *reinterpret_cast<const uint16_t*>(val.c_str());
        matom.set_owner(vid, owner);
      }
      else {
        cur->step();
      }
    }
    delete cur;
  }

} // namespace graphlab

