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
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/graph/disk_atom.hpp>
#include <kchashdb.h>


namespace graphlab {

  disk_atom::disk_atom(std::string filename, 
                       uint16_t atomid, bool constant_in_memory):
                          atomid(atomid),filename(filename),const_in_mem(constant_in_memory) {
    if (const_in_mem) {
      // try to load the fast file
      std::ifstream fin((filename+".fast").c_str(), std::ios::binary);
      if (fin.is_open() && fin.good()) {
        iarchive iarc(fin);
        iarc >> cache;
        if (cache_get("numv", 4, (char*)&numv.value, sizeof(numv.value)) == -1) numv = 0;
        if (cache_get("nume", 4, (char*)&nume.value, sizeof(nume.value)) == -1) nume = 0;
        if (cache_get("numlocalv", 9, (char*)&numlocalv.value, sizeof(numlocalv.value)) == -1) numlocalv = 0;
        if (cache_get("numlocale", 9, (char*)&numlocale.value, sizeof(numlocale.value)) == -1) numlocale = 0;

      }
      else {
        open_db();
        build_cache();
      }
    }
    else {
      open_db();
    }
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
    if (db.get("numlocalv", 9, (char*)&numlocalv.value, sizeof(numlocalv.value)) == -1) numlocalv = 0;
    if (db.get("numlocale", 9, (char*)&numlocale.value, sizeof(numlocale.value)) == -1) numlocale = 0;
  }
  

  disk_atom::~disk_atom() {
    synchronize();
    db.close();
  }

  void disk_atom::synchronize() {
    // update the linked list
    db.set("numv", 4, (char*)&numv.value, sizeof(numv.value));
    db.set("nume", 4, (char*)&nume.value, sizeof(nume.value));
    db.set("numlocalv", 9, (char*)&numlocalv.value, sizeof(numlocalv.value));
    db.set("numlocale", 9, (char*)&numlocale.value, sizeof(numlocale.value));
    db.synchronize();
  }

  void disk_atom::add_vertex(disk_atom::vertex_id_type vid, uint16_t owner) {
    ASSERT_FALSE(const_in_mem);
    if (!add_vertex_skip(vid, owner)) {
      std::stringstream strm;
      oarchive oarc(strm);    
      oarc << owner;
      strm.flush();
      db.set("v"+id_to_str(vid), strm.str());
    }
  }


  bool disk_atom::add_vertex_skip(disk_atom::vertex_id_type vid, uint16_t owner) {
    ASSERT_FALSE(const_in_mem);
    mut[vid % 511].lock();
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner;
    strm.flush();
    if (db.add("v"+id_to_str(vid), strm.str())) {
      uint64_t v64 = (uint64_t)vid;
      db.append("_vidlist", 8, (char*)&v64, sizeof(v64));
      numv.inc();
      if (owner == atomid) numlocalv.inc();
      mut[vid % 511].unlock();
      return true;
    }
    mut[vid % 511].unlock();
    return false;
  }


  void disk_atom::add_edge(disk_atom::vertex_id_type src, disk_atom::vertex_id_type target) {
    ASSERT_FALSE(const_in_mem);
    if (!add_edge_skip(src, target)) {
      db.set("e"+id_to_str(src)+"_"+id_to_str(target), std::string(""));
    }
  }


  bool disk_atom::add_edge_skip(disk_atom::vertex_id_type src, disk_atom::vertex_id_type target) {
    ASSERT_FALSE(const_in_mem);
    mut[(src ^ target) % 511].lock();
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
      mut[(src ^ target) % 511].unlock();
      return true;
    }
    mut[(src ^ target) % 511].unlock();
    return false;
  }

  std::vector<disk_atom::vertex_id_type> disk_atom::enumerate_vertices() {
    std::vector<disk_atom::vertex_id_type> ret;
    // read the entire vertex list
    std::string vidlist;
    if (const_in_mem && cache_get(std::string("_vidlist"), &vidlist) == false) {
      return ret;
    }
    else if (db.get(std::string("_vidlist"), &vidlist) == false) {
      return ret;
    }
    
    ASSERT_EQ(vidlist.size() % sizeof(uint64_t), 0);
    ret.resize(vidlist.size() / sizeof(uint64_t));
    const uint64_t* arr = reinterpret_cast<const uint64_t*>(vidlist.c_str());
    std::copy(arr, arr + ret.size(), ret.begin());
    return ret;
  }


  bool disk_atom::get_vertex(disk_atom::vertex_id_type vid, uint16_t &owner) {
    std::string val;
    std::string key = "v"+id_to_str(vid);
    
    if (const_in_mem && cache_get(key, &val) == false) {
      return false;
    }
    else if (db.get("v"+id_to_str(vid), &val) == false) {
      return false;
    }

    boost::iostreams::stream<boost::iostreams::array_source> 
      istrm(val.c_str(), val.length());   
    iarchive iarc(istrm);
    iarc >> owner;
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
    
    if ((const_in_mem && cache_get(key, &val)) || db.get(key, &val)) {
      const uint64_t* v = reinterpret_cast<const uint64_t*>(val.c_str());
      ASSERT_TRUE(val.length() % 8 == 0);
      size_t numel = val.length() / 8;
      ret.resize(numel);
      for (size_t i = 0;i < numel; ++i) ret[i] = v[i];
    }
    return ret;    
  }
   
   

  std::vector<disk_atom::vertex_id_type> disk_atom::get_out_vertices(disk_atom::vertex_id_type vid) {
    std::vector<disk_atom::vertex_id_type> ret;
    std::string val;
    std::string key = "o"+id_to_str(vid);
    if ((const_in_mem && cache_get(key, &val)) || db.get(key, &val)) {
      const uint64_t* v = reinterpret_cast<const uint64_t*>(val.c_str());
      size_t numel = val.length() / 8;
      ASSERT_TRUE(val.length() % 8 == 0);
      ret.resize(numel);
      for (size_t i = 0;i < numel; ++i) ret[i] = v[i];
    }
    return ret;    
  }



  disk_atom::vertex_color_type 
  disk_atom::get_color(disk_atom::vertex_id_type vid) {
    std::string key = "c" + id_to_str(vid);
    disk_atom::vertex_color_type  ret;
    if (const_in_mem && 
        cache_get(key.c_str(), key.length(), 
                  (char*)&ret, sizeof(ret)) != -1) return ret;

    if (db.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) == -1) 
      ret = disk_atom::vertex_color_type(-1);
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
    if (const_in_mem && 
        cache_get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) != -1) return ret;
  
    if (db.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) == -1) ret = (uint16_t)(-1);
    return ret;
  }


  void disk_atom::set_owner(disk_atom::vertex_id_type vid, uint16_t owner) {
    std::string key = "h" + id_to_str(vid);
    db.set(key.c_str(), key.length(), (char*)&owner, sizeof(owner));
  }


  void disk_atom::clear() {
    numv.value = 0;
    nume.value = 0;
    numlocalv.value = 0;
    numlocale.value = 0;
    db.clear();
  }

  void disk_atom::build_cache() {
    cache.clear();
    storage_type::Cursor* cur = db.cursor();
    cur->jump();
    std::string key, val;
    while (cur->get(&key, &val, true)) {
      cache[key] = val;
    }
  }
  
  void disk_atom::build_fast_finalized_atom() {
    build_cache();
    std::string finalizedatom = filename + ".fast";  

    std::ofstream fout(finalizedatom.c_str(), std::ios::binary);
    oarchive oarc(fout);
    oarc << cache;
    fout.close();
  }

} // namespace graphlab

