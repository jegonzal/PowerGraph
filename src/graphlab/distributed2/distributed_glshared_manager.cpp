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


#include <graphlab/distributed2/distributed_glshared.hpp>
#include <graphlab/distributed2/distributed_glshared_manager.hpp>



namespace graphlab {

distributed_glshared_manager::distributed_glshared_manager(distributed_control &dc):
                  rmi(dc, this),
                  glsharedobjs(distgl_impl::get_global_dist_glshared_registry()),
                  dht(dc){
  dht.attach_modification_trigger(boost::bind(&distributed_glshared_manager::invalidate,
                                              this, _1, _2, _3));
  for (size_t i = 0; i < glsharedobjs.size(); ++i) {
    logstream(LOG_INFO) << "registered entry " << i << " with type " 
                        << glsharedobjs[i]->type_name() << std::endl;
    if (glsharedobjs[i]->manager != NULL) {
      logger(LOG_WARNING, "glshared objects are still attached to a previous manager!");
    }
    glsharedobjs[i]->manager = this;
    glsharedobjs[i]->id = i;
    objrevmap[glsharedobjs[i]] = i;
    if (dht.owning_machine(i) == rmi.procid()) {
      std::stringstream strm;
      oarchive oarc(strm);
      glsharedobjs[i]->save(oarc);
      dht.set(i, strm.str());
    }
  }
  // perform the sets
}

distributed_glshared_manager::~distributed_glshared_manager() {
  // detach from the dht
  dht.detach_modification_trigger();
  // detach from the glshared objects
  for (size_t i = 0; i < glsharedobjs.size(); ++i) {
    glsharedobjs[i]->manager = NULL;
  }
}
/*
completes an atomic exchange of an entry
*/
std::string distributed_glshared_manager::exchange(size_t entry, 
                                                   const std::string &val) {
  std::string ret;
  if (dht.owning_machine(entry) == rmi.procid()) {
    std::string& valref = dht.begin_critical_section(entry);
    ret = valref;
    valref = val;
    dht.end_critical_section(entry);
    dht.push_changes(entry, false, procid_t(-1));
  }
  else {
    ret = rmi.remote_request(dht.owning_machine(entry),
                             &distributed_glshared_manager::exchange,
                             entry,
                             val);
  }
  dht.invalidate(entry);
  return ret;
}


void distributed_glshared_manager::invalidate(size_t entry, const std::string &value, bool incache) {
  if (incache) {
//    read_synchronize(entry, false);
    std::stringstream strm(value);
    iarchive iarc(strm);
    glsharedobjs[entry]->load(iarc);
  }
  else {
    glsharedobjs[entry]->invalidated = true;
  }
}


/**
Synchronize variable with index i.
Call
*/
void distributed_glshared_manager::write_synchronize(size_t entry, bool async) {
//  logstream(LOG_DEBUG) << rmi.procid() << ": " << "write synchronize on " << entry << " async = " << async << std::endl;
  std::stringstream strm;
  oarchive oarc(strm);
  glsharedobjs[entry]->save(oarc);
  glsharedobjs[entry]->invalidated = false;
  if (async) {
    dht.set(entry, strm.str());
  }
  else {
    dht.set_synchronous(entry, strm.str());
  }
}

void distributed_glshared_manager::write_synchronize(distributed_glshared_base* obj, 
                                                    bool async) {
  write_synchronize(objrevmap[obj], async);
}


void distributed_glshared_manager::read_synchronize(size_t entry, 
                                                    bool async) {
//  logstream(LOG_DEBUG) << rmi.procid() << ": " << "read synchronize on " << entry << " async = " << async << std::endl;
  if (async) {
    bool ret = dht.asynchronous_get(entry);
    if (ret == false) return;
    // if ret = true, it is cache and we pass through 
    // to the synchronous read
  }
  
  std::pair<bool, std::string> ret = dht.get(entry);
  assert(ret.first);
  std::stringstream strm(ret.second);
  iarchive iarc(strm);
  glsharedobjs[entry]->load(iarc);
  glsharedobjs[entry]->invalidated = false;
}

void distributed_glshared_manager::read_synchronize(distributed_glshared_base* obj, 
                                                    bool async) {
  read_synchronize(objrevmap[obj], async);
}

} // namespace graphlab

