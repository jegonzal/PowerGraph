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


// #ifndef FROM_DISTRIBUTED_GRAPH_INCLUDE
// #warning "distributed_graph_synchronization.hpp should not be included directly."
// #warning "You should include only distributed_graph.hpp"
// #warning "I will fix this for you now, but don't do it again!"

// #include <graphlab/distributed2/graph/distributed_graph.hpp>
// #include <boost/range/iterator_range.hpp>
// #else

// //#define DGRAPH_DEBUG

// /**
//  * synchronize the data on vertex with global id vid
//  * vid must be a ghost
//  */
// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::synchronize_vertex(vertex_id_t vid, bool async) {
//   vertex_id_t localvid = global2localvid[vid];
//   if (is_ghost(vid)) {
//     vertex_conditional_store out;
//     out.hasdata = localstore.vertex_modified(localvid);
//     if (out.hasdata) {
//       localstore.set_vertex_modified(localvid, false);
//       out.data.first = localstore.vertex_data(localvid);
//     }
//     if (async == false) {
//       vertex_conditional_store v;
//       v = rmi.remote_request(localvid2owner[localvid],
//                             &distributed_graph<VertexData, EdgeData>::get_vertex_if_version_less_than,
//                             vid,
//                             localstore.vertex_version(localvid),
//                             out);
//      if (v.hasdata) {
//         update_vertex_data_and_version(vid, v);
//       }
//     }
//     else {
//       pending_async_updates.flag.inc();
//       rmi.remote_call(localvid2owner[localvid],
//                        &distributed_graph<VertexData, EdgeData>::async_get_vertex_if_version_less_than,
//                        rmi.procid(),
//                        vid,
//                        localstore.vertex_version(localvid),
//                        out);
//     }
//   }
// }


// /**
//  * synchronize the data on edge with global id eid
//  * target of edge must be a ghost
//  */
// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::synchronize_edge(edge_id_t eid, bool async) {
//   vertex_id_t localtargetvid = localstore.target(eid);
//   vertex_id_t targetvid = local2globalvid[localtargetvid];
//   vertex_id_t localsourcevid = localstore.source(eid);
//   vertex_id_t sourcevid = local2globalvid[localsourcevid];

//   if (is_ghost(targetvid)) {

//     edge_conditional_store out;
//     out.hasdata = localstore.edge_modified(eid);
//     if (out.hasdata) {
//       localstore.set_edge_modified(eid, false);
//       out.data.first = localstore.edge_data(eid);
//     }
//     if (async == false) {
//       edge_conditional_store e = rmi.remote_request(localvid2owner[localtargetvid],
//                                                     &distributed_graph<VertexData, EdgeData>::get_edge_if_version_less_than2,
//                                                     sourcevid,
//                                                     targetvid,
//                                                     localstore.edge_version(eid),
//                                                     out);
//       if (e.hasdata) {
//         update_edge_data_and_version2(sourcevid, targetvid, e);
//       }
//     }
//     else {
//       pending_async_updates.flag.inc();
//       rmi.remote_call(localvid2owner[localtargetvid],
//                         &distributed_graph<VertexData, EdgeData>::async_get_edge_if_version_less_than2,
//                         rmi.procid(),
//                         sourcevid,
//                         targetvid,
//                         localstore.edge_version(eid),
//                         out);
//     }
//   }
// }

// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::synchronize_scope(vertex_id_t vid, bool async) {
//   // construct he requests
//   std::map<procid_t, std::pair<block_synchronize_request2, std::vector<vertex_id_t>::iterator> > requests;
//   synchronize_scope_construct_req(vid, requests);
  
//   if (async) {
//     // if asynchronous, the reply goes to pending_async_updates
//     typename std::map<procid_t, 
//           std::pair<block_synchronize_request2, std::vector<vertex_id_t>::iterator> >::iterator iter;
//     iter = requests.begin();
//     size_t replytarget = reinterpret_cast<size_t>(&pending_async_updates);
//     pending_async_updates.flag.inc(requests.size());

//     while(iter != requests.end()) {
//       rmi.remote_call(iter->first,
//                       &distributed_graph<VertexData, EdgeData>::async_get_alot2,
//                       rmi.procid(),
//                       iter->second.first,
//                       replytarget,
//                       0);
//       ++iter;
//     }
//   }
//   else {
//     // otherwise we collect it into a local reply ret tye
//     dc_impl::reply_ret_type reply(true, 0);
//     typename std::map<procid_t, 
//           std::pair<block_synchronize_request2, std::vector<vertex_id_t>::iterator> >::iterator iter;
//     iter = requests.begin();
//     size_t replytarget = reinterpret_cast<size_t>(&reply);
//     reply.flag.inc(requests.size());

//     while(iter != requests.end()) {
//       rmi.remote_call(iter->first,
//                       &distributed_graph<VertexData, EdgeData>::async_get_alot2,
//                       rmi.procid(),
//                       iter->second.first,
//                       replytarget,
//                       0);
//       ++iter;
//     }

//     reply.wait();
//   }
  
  
  
  
  
// }



// /**
//  * Constructs the request for synchronizing the scope for vertex vid
//  * vid must be owned by the current machine. 
//  */
// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::synchronize_scope_construct_req(vertex_id_t vid, 
//       std::map<procid_t, std::pair<block_synchronize_request2, std::vector<vertex_id_t>::iterator> > &requests) {
//   ASSERT_FALSE(is_ghost(vid));
//   if (boundaryscopesset.find(vid) == boundaryscopesset.end()) return;

//   // now this is very annoying. A significant amount of code is identical here.
//   // whether with edge canonical numbering on or not. But I cannot refactor it
//   // easily because the types are different and I don't really want to introduce
//   // templates here for something so trivial.
  
//     vertex_id_t localvid = global2localvid[vid];
//     requests.clear();
//     // I should have all the in-edges. but I need the in vertices.
//     // need to track the vertices added so I don't add duplicate vertices
//     // if the vertex has both in-out edges to this vertex.
//     // trick! vids are ordered!

//   foreach(edge_id_t localineid, localstore.in_edge_ids(localvid)) {
//     vertex_id_t localsourcevid = localstore.source(localineid);
//     if (localvid_is_ghost(localsourcevid)) {
//       // need to synchronize incoming vertex
//       procid_t sourceowner = localvid2owner[localsourcevid];
//       block_synchronize_request2 &req = requests[sourceowner].first;
//       req.vid.push_back(local2globalvid[localsourcevid]);
//       req.vidversion.push_back(localstore.vertex_version(localsourcevid));
//       vertex_conditional_store vs;
//       vs.hasdata = localstore.vertex_modified(localsourcevid);
//       if (vs.hasdata) {
//         localstore.set_vertex_modified(localsourcevid, false);
//         vs.data.first = localstore.vertex_data(localsourcevid);
//         vs.data.second = localstore.vertex_version(localsourcevid);
//       }
//       requests[sourceowner].second=req.vid.end();
//       req.vstore.push_back(vs);
//     }
//   }
//   // now for the out edges
//   foreach(edge_id_t localouteid, localstore.out_edge_ids(localvid)) {
//     vertex_id_t localtargetvid = localstore.target(localouteid);
//     procid_t targetowner = localvid2owner[localstore.target(localouteid)];

//     if (localvid_is_ghost(localtargetvid)) {
//       block_synchronize_request2 &req = requests[targetowner].first;
//       // need to synchronize outgoing vertex and outgoing edge
//       // do outgoing vertex first
//       if (std::binary_search(req.vid.begin(), requests[targetowner].second, local2globalvid[localtargetvid]) == false) {
//         req.vid.push_back(local2globalvid[localtargetvid]);
//         req.vidversion.push_back(localstore.vertex_version(localtargetvid));
//         vertex_conditional_store vs;
//         vs.hasdata = localstore.vertex_modified(localtargetvid);
//         if (vs.hasdata) {
//           localstore.set_vertex_modified(localtargetvid, false);
//           vs.data.first = localstore.vertex_data(localtargetvid);
//           vs.data.second = localstore.vertex_version(localtargetvid);
//         }
//         req.vstore.push_back(vs);
//       }
//       // now for the outgoing edge
      
//       req.srcdest.push_back(std::pair<vertex_id_t, vertex_id_t>(vid, local2globalvid[localtargetvid]));
//       req.edgeversion.push_back(localstore.edge_version(localouteid));
//       edge_conditional_store es;
//       es.hasdata = localstore.edge_modified(localouteid);
//       if (es.hasdata) {
//         localstore.set_edge_modified(localouteid, false);
//         es.data.first = localstore.edge_data(localouteid);
//         es.data.second = localstore.edge_version(localouteid);
//       }
//       req.estore.push_back(es);
//     }
//   }
// }



// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::async_synchronize_scope_callback(vertex_id_t vid, 
//                                                 boost::function<void (void)> callback){
//   vertex_id_t localvid = global2localvid[vid];
//   ASSERT_TRUE(scope_callbacks[localvid].callback == NULL);
//   ASSERT_TRUE(boundary_scopes_set().find(vid) != boundary_scopes_set().end());
//   // register the callback
//   scope_callbacks[localvid].callback = callback;
//   // construct the requests
//   std::map<procid_t, std::pair<block_synchronize_request2, std::vector<vertex_id_t>::iterator> > requests;
//   synchronize_scope_construct_req(vid, requests);
  
//   //send the stuff
//   typename std::map<procid_t, 
//         std::pair<block_synchronize_request2, std::vector<vertex_id_t>::iterator> >::iterator iter;
//   iter = requests.begin();
  
//   ASSERT_EQ(scope_callbacks[localvid].counter.value, 0);
//   scope_callbacks[localvid].counter.inc(requests.size());
//   while(iter != requests.end()) {

//     // the reply target is 0. see reply_alot2
//     rmi.remote_call(iter->first,
//                     &distributed_graph<VertexData, EdgeData>::async_get_alot2,
//                     rmi.procid(),
//                     iter->second.first,
//                     0,
//                     localvid);
//     ++iter;
//   }


// }

// /**
//  * Waits for all asynchronous data synchronizations to complete
//  */
// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::wait_for_all_async_syncs() {
//   pending_async_updates.wait();
// }



// template <typename VertexData, typename EdgeData> 
// typename distributed_graph<VertexData, EdgeData>::vertex_conditional_store 
// distributed_graph<VertexData, EdgeData>::get_vertex_if_version_less_than(vertex_id_t vid,
//                                                          uint64_t  vertexversion,
//                                                          vertex_conditional_store &vdata) {
//   vertex_conditional_store ret;
//   size_t localvid = global2localvid[vid];
//   uint64_t local_vertex_version = localstore.vertex_version(localvid);
//   // Now I must the the owner of this vertex
//   ASSERT_EQ(localvid2owner[localvid], rmi.procid());
//   #ifdef DGRAPH_DEBUG
//   logstream(LOG_DEBUG) << "get vertex: " << vid << ":" << vertexversion << " vs " << local_vertex_version << ". " << vdata.hasdata << std::endl;
//   #endif
//   ret.hasdata = false;
  
//   if (local_vertex_version  > vertexversion) {
//     // send if I have a later version
//     ret.hasdata = true;
//     ret.data.first = localstore.vertex_data(localvid);
//     ret.data.second = local_vertex_version;
//   }
//   else if (local_vertex_version == vertexversion) {
//     // if version is the same and there is data, store and increment the version    
//     if (vdata.hasdata) {
//       localstore.increment_and_update_vertex(localvid, vdata.data.first);
//     }
//   }
//   else {
//     logstream(LOG_FATAL) << "Remote node attempted to update vertex " 
//                         << vid << " with a newer version than the owner" << std::endl;
//   }
//   return ret;
// }


// template <typename VertexData, typename EdgeData> 
// typename distributed_graph<VertexData, EdgeData>::edge_conditional_store 
// distributed_graph<VertexData, EdgeData>::get_edge_if_version_less_than2(vertex_id_t source,
//                                                       vertex_id_t target,
//                                                       uint64_t  edgeversion,
//                                                       edge_conditional_store &edata) {
//   edge_conditional_store ret;
//   size_t localsource = global2localvid[source];
//   size_t localtarget = global2localvid[target];


//   std::pair<bool, edge_id_t> findret = localstore.find(localsource, localtarget);
//   assert(findret.first);
//   edge_id_t localeid = findret.second;
  
//   uint64_t  local_edge_version = localstore.edge_version(localeid);
//   ret.hasdata = false;
//   #ifdef DGRAPH_DEBUG
  
//   logstream(LOG_DEBUG) << "get edge2: " << "(" << source << "->" << target << ")" << ":" << edgeversion << " vs " << local_edge_version << ". " << edata.hasdata << std::endl;
//   #endif
//   if (local_edge_version > edgeversion) {
//     ret.hasdata = true;
//     ret.data.first = localstore.edge_data(localeid);
//     ret.data.second = local_edge_version;
//   }
//   else if (local_edge_version == edgeversion) {
//     if (edata.hasdata) {
//       localstore.increment_and_update_edge(localeid, edata.data.first);
//     }
//   }
//   else {
//     logstream(LOG_FATAL) << "Remote node attempted to update edge (" 
//                         << source <<  ", " << target << ") with a newer version than the owner" << std::endl;
//   }
//   return ret;
// }


// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::async_get_vertex_if_version_less_than(
//                                                               procid_t srcproc, 
//                                                               vertex_id_t vid, 
//                                                               uint64_t vertexversion,
//                                                               vertex_conditional_store &vdata) {
//   rmi.remote_call(srcproc,
//                   &distributed_graph<VertexData, EdgeData>::reply_vertex_data_and_version,
//                   vid,
//                   get_vertex_if_version_less_than(vid, vertexversion, vdata));
// }

// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::async_get_edge_if_version_less_than2(
//                                               procid_t srcproc, 
//                                               vertex_id_t source, 
//                                               vertex_id_t target, 
//                                               uint64_t edgeversion,
//                                               edge_conditional_store &edata) {
//   rmi.remote_call(srcproc,
//                   &distributed_graph<VertexData, EdgeData>::reply_edge_data_and_version2,
//                   source,
//                   target,
//                   get_edge_if_version_less_than2(source, target, edgeversion, edata));
// }


// template <typename VertexData, typename EdgeData> 
// typename distributed_graph<VertexData, EdgeData>::block_synchronize_request2& 
// distributed_graph<VertexData, EdgeData>::get_alot2(
//                           distributed_graph<VertexData, EdgeData>::block_synchronize_request2 &request) {
//   std::vector<vertex_conditional_store> vresponse(request.vid.size());
//   std::vector<edge_conditional_store> eresponse(request.srcdest.size());
//   for (size_t i = 0;i < request.vid.size(); ++i) {
//     request.vstore[i] = get_vertex_if_version_less_than(request.vid[i], 
//                                                         request.vidversion[i], 
//                                                         request.vstore[i]);
//   }
//   for (size_t i = 0;i < request.srcdest.size(); ++i) {
//     request.estore[i] = get_edge_if_version_less_than2(request.srcdest[i].first, 
//                                                        request.srcdest[i].second, 
//                                                        request.edgeversion[i], 
//                                                        request.estore[i]);
//   }
//   request.vidversion.clear();
//   request.edgeversion.clear();
//   return request;
// }

// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::async_get_alot2(
//                     procid_t srcproc,
//                     distributed_graph<VertexData, EdgeData>::block_synchronize_request2 &request,
//                     size_t replytarget,
//                     size_t tag) {
//   get_alot2(request);
//   rmi.remote_call(srcproc,
//                   &distributed_graph<VertexData, EdgeData>::reply_alot2,
//                   request,
//                   replytarget, 
//                   tag);
// }

// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::reply_vertex_data_and_version(
//                         vertex_id_t vid, 
//                         distributed_graph<VertexData, EdgeData>::vertex_conditional_store &vstore) {
//   if (vstore.hasdata) update_vertex_data_and_version(vid, vstore);
//   // the dc and procid are meaningless. Just pass something
//   reply_increment_counter(rmi.dc(), 0, 
//                           reinterpret_cast<size_t>(&pending_async_updates), dc_impl::blob());
// }

// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::reply_edge_data_and_version2(
//                   vertex_id_t source, 
//                   vertex_id_t target, 
//                   distributed_graph<VertexData, EdgeData>::edge_conditional_store &estore) {
//   if (estore.hasdata) update_edge_data_and_version2(source, target, estore);
//   reply_increment_counter(rmi.dc(), 0, 
//                           reinterpret_cast<size_t>(&pending_async_updates), dc_impl::blob());

// }



// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::update_vertex_data_and_version(
//                         vertex_id_t vid, 
//                         distributed_graph<VertexData, EdgeData>::vertex_conditional_store &vstore) {
//   vertex_id_t localvid = global2localvid[vid];
//   // this must be a ghost
//   ASSERT_NE(localvid2owner[localvid], rmi.procid());
// #ifdef DGRAPH_DEBUG
//   logstream(LOG_DEBUG) << "Receiving vertex " << vid << "(" << localvid << ")  "  << ". "
//                        << vstore.data.second << " vs " << localstore.vertex_version(localvid) << std::endl;
// #endif
//   localstore.conditional_update_vertex(localvid, vstore.data.first, vstore.data.second);
// }

// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::update_edge_data_and_version2(
//                           vertex_id_t source, 
//                           vertex_id_t target, 
//                           distributed_graph<VertexData, EdgeData>::edge_conditional_store &estore) {
//   if (estore.hasdata) {
//     vertex_id_t localsourcevid = global2localvid[source];
//     vertex_id_t localtargetvid = global2localvid[target];
//     ASSERT_NE(localvid2owner[localtargetvid], rmi.procid());
//     std::pair<bool, edge_id_t> findret = localstore.find(localsourcevid, localtargetvid);
    
//     assert(findret.first);
//     #ifdef DGRAPH_DEBUG
//     logstream(LOG_DEBUG) << "Receiving edge (" << source << ","<<target << ")  " << ". "
//                        << estore.data.second << " vs " << localstore.edge_version(findret.second) << std::endl;
//     #endif

//     localstore.conditional_update_edge(findret.second, estore.data.first, estore.data.second);
//   }
// }


// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::update_alot2(
//                   distributed_graph<VertexData, EdgeData>::block_synchronize_request2 &request) {
//   for (size_t i = 0;i < request.vid.size(); ++i) {
//     update_vertex_data_and_version(request.vid[i], request.vstore[i]);
//   }

//   for (size_t i = 0;i < request.srcdest.size(); ++i) {
//     update_edge_data_and_version2(request.srcdest[i].first, request.srcdest[i].second, request.estore[i]);
//   }
// }


// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::reply_alot2(
//                     distributed_graph<VertexData, EdgeData>::block_synchronize_request2 &request,
//                     size_t replytarget,
//                     size_t tag) {
//   update_alot2(request);
  
//   // special handling for callbacks
//   if (replytarget != 0)  {
//     reply_increment_counter(rmi.dc(), 0, 
//                           replytarget, dc_impl::blob());  
//   }
//   else {
//     // tag is local vid
//     vertex_id_t localvid = tag;

//     ASSERT_TRUE(scope_callbacks[localvid].callback != NULL);
//     if (scope_callbacks[localvid].counter.dec() == 0) {
//       // make a copy of it and clear the callbacks entry.
//       boost::function<void (void)> tmp = scope_callbacks[localvid].callback;
//       scope_callbacks[localvid].callback = NULL;
//       tmp();
//     }
//   } 
// }


// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::synchronize_all_vertices(bool async) {
//   foreach(vertex_id_t vid, ghostvertices) {
//     synchronize_vertex(vid, async);
//   }
// }

// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::synchronize_all_edges(bool async) {
//   foreach(vertex_id_t vid, ghostvertices) {
//     foreach(edge_id_t eid, localstore.in_edge_ids(global2localvid[vid])) {
//       synchronize_edge(eid, async);
//     }
//   }
// }

// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::synchronize_all_scopes(bool async) {
//   foreach(vertex_id_t vid, boundaryscopes) {
//     synchronize_scope(vid, async);
//   }
// }



// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::update_vertex_data_and_version_and_reply(
//                         vertex_id_t vid, 
//                         distributed_graph<VertexData, EdgeData>::vertex_conditional_store &vstore,
//                         procid_t srcproc,
//                         size_t reply) {
                        
//   update_vertex_data_and_version(vid, vstore);

//   if (srcproc != procid_t(-1)) {
//     rmi.dc().remote_call(srcproc, reply_increment_counter, reply, dc_impl::blob());
//   }
// }

// template <typename VertexData, typename EdgeData> 
// void distributed_graph<VertexData, EdgeData>::update_edge_data_and_version_and_reply2(
//                           vertex_id_t source, 
//                           vertex_id_t target, 
//                           distributed_graph<VertexData, EdgeData>::edge_conditional_store &estore,
//                           procid_t srcproc, size_t reply) {
//   update_edge_data_and_version2(source, target, estore);
//   if (srcproc != procid_t(-1)) {
//     rmi.dc().remote_call(srcproc, reply_increment_counter, reply, dc_impl::blob());
//   }
// }

// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::push_owned_vertex_to_replicas(vertex_id_t vid, bool async, bool untracked) {
//   vertex_id_t localvid = global2localvid[vid];
//   // get the replicas
//   const std::vector<procid_t>& replicas = localvid_to_replicas(localvid);
//   // owner is a replica too. if there are no other replicas quit
//   if (replicas.size() <= 1) return;
  
//   dc_impl::reply_ret_type ret(true, replicas.size() - 1);
  
//   // if async, set the return reply to go to the global pending push updates
//   size_t retptr;
//   if (async == false) {
//     retptr = reinterpret_cast<size_t>(&ret);
//   }
//   else {
//     retptr = reinterpret_cast<size_t>(&pending_push_updates);
//     pending_push_updates.flag.inc(replicas.size() - 1);
//   }
  
//   /**
//   If untracked, set the reply procid to -1. That will mean that
//   no reply will be returned at all
//   */
//   procid_t srcprocid;
//   if (untracked == false) srcprocid = rmi.procid();
//   else srcprocid = procid_t(-1);
//   // build the store
//   vertex_conditional_store vstore;
//   vstore.hasdata = true;
//   vstore.data.first = localstore.vertex_data(localvid);
//   vstore.data.second = localstore.vertex_version(localvid);
  
//  foreach(procid_t proc, replicas) {
//     if (proc != rmi.procid()) {
// #ifdef DGRAPH_DEBUG
//       logger(LOG_DEBUG, "Pushing vertex %d to proc %d", vid, proc);
// #endif
//       rmi.remote_call(proc,
//                       &distributed_graph<VertexData, EdgeData>::update_vertex_data_and_version_and_reply,
//                       vid,
//                       vstore,
//                       srcprocid,
//                       retptr);
//     }
//   }
//   if (async == false && untracked == false)  ret.wait();
// }


// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::push_owned_edge_to_replicas(edge_id_t eid, bool async, bool untracked) {
//   // get the replicas
//   vertex_id_t globalsource = source(eid);
//   vertex_id_t globaltarget = target(eid);
//   // firstly I must own this edge
//   if (!is_owned(globaltarget)) return;
  
//   // Now, there are 2 cases. 
//   // Case 1, I own both source and target.
//   //         in that case there is nothing to sync. Any other
//   //         machine at most has ghosts of source and target, but will
//   //         not have the end itself
//   // Case 2, I own the target, but not the source
//   //         then the only replica is the owner of the source

//   procid_t sendto;
  
//   if (is_owned(globalsource)) {
//     return;
//   }
//   else {
//     // otherwise, it is the owner of the source
//     sendto = localvid2owner[localstore.source(eid)];
//   }

//   dc_impl::reply_ret_type ret(true, 1);  
//   // if async, set the return reply to go to the global pending push updates
//   size_t retptr;
//   if (async == false) {
//     retptr = reinterpret_cast<size_t>(&ret);
//   }
//   else {
//     retptr = reinterpret_cast<size_t>(&pending_push_updates);
//   }
  
//   /**
//   If untracked, set the reply procid to -1. That will mean that
//   no reply will be returned at all
//   */
//   procid_t srcprocid;
//   if (untracked == false) srcprocid = rmi.procid();
//   else srcprocid = procid_t(-1);
//   // build the store
//   edge_conditional_store estore;
//   estore.hasdata = true;
//   estore.data.first = localstore.edge_data(eid);
//   estore.data.second = localstore.edge_version(eid);
  
// #ifdef DGRAPH_DEBUG
//   logstream(LOG_DEBUG) << "Pushing edge (" << globalsource << ", " << globaltarget << ") to proc " << sendto << std::endl;
// #endif
//   rmi.remote_call(sendto,
//                   &distributed_graph<VertexData, EdgeData>::update_edge_data_and_version_and_reply2,
//                   globalsource,
//                   globaltarget,
//                   estore,
//                   srcprocid,
//                   retptr);

//   if (async == false && untracked == false)  ret.wait();
// }


// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::push_all_owned_vertices_to_replicas() {
//   // perform block collective
//   std::vector<block_synchronize_request2> blockpushes(rmi.numprocs()); 

//   foreach(vertex_id_t vid, boundary_scopes()) {

//     vertex_id_t localvid = global2localvid[vid];
//     // get the replicas
//     const std::vector<procid_t>& replicas = localvid_to_replicas(localvid);
//     // owner is a replica too. if there are no other replicas quit
//     if (replicas.size() <= 1) continue;
  
 
//    // build the store
//     vertex_conditional_store vstore;
//     vstore.hasdata = true;
//     vstore.data.first = localstore.vertex_data(localvid);
//     vstore.data.second = localstore.vertex_version(localvid);
  
//     foreach(procid_t proc, replicas) {
//       if (proc != rmi.procid()) {
//         blockpushes[proc].vid.push_back(vid);
//         blockpushes[proc].vidversion.push_back(localstore.vertex_version(localvid));
//         blockpushes[proc].vstore.push_back(vstore);
//         if (blockpushes[proc].vid.size() >= 1024*1024/sizeof(VertexData)) {
//           rmi.remote_call(proc,
//                           &distributed_graph<VertexData, EdgeData>::update_alot2,
//                           blockpushes[proc]);
//           blockpushes[proc].clear();
//         }
//       }
//     }
//   }

//   for(size_t proc = 0; proc < rmi.numprocs(); ++proc) {
//     if (blockpushes[proc].vid.size() > 0) {
//       assert(proc != rmi.procid());
//       rmi.remote_call(proc,
//                       &distributed_graph<VertexData, EdgeData>::update_alot2,
//                       blockpushes[proc]);
//       blockpushes[proc].clear();
//     }
//   }
// }

// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::push_all_owned_edges_to_replicas() {
//   std::vector<std::vector<block_synchronize_request2> > blockpushes(omp_get_max_threads()); 
//   for (size_t i = 0;i < blockpushes.size(); ++i) blockpushes[i].resize(rmi.numprocs());
  
  
//   #pragma omp parallel for
//   for (long i = 0;i < (long)ghostvertices.size(); ++i) {
//     int thrid = omp_get_thread_num();
    
//     vertex_id_t vid = ghost_vertices()[i];
//     vertex_id_t localvid = global2localvid[vid];
//     procid_t proc = localvid2owner[localvid];
//     foreach(edge_id_t localeid, localstore.out_edge_ids(localvid)) {
//       vertex_id_t targetvid = local2globalvid[localstore.target(localeid)];
//       edge_conditional_store estore;
//       estore.hasdata = true;
//       estore.data.first = localstore.edge_data(localeid);
//       estore.data.second = localstore.edge_version(localeid);
  

//       blockpushes[thrid][proc].srcdest.push_back(std::make_pair<vertex_id_t, vertex_id_t>(vid, targetvid));
//       blockpushes[thrid][proc].edgeversion.push_back(localstore.edge_version(localeid));
//       blockpushes[thrid][proc].estore.push_back(estore);
//       if (blockpushes[thrid][proc].srcdest.size() >= 1*1024*1024/sizeof(EdgeData)) {
//         rmi.remote_call(proc,
//                         &distributed_graph<VertexData, EdgeData>::update_alot2,
//                         blockpushes[thrid][proc]);
//         blockpushes[thrid][proc].clear();
//       }
//     }
//   }
//   for (size_t i = 0;i < blockpushes.size(); ++i) {
//     for(size_t proc = 0; proc < rmi.numprocs(); ++proc) {
//       if (blockpushes[i][proc].srcdest.size() > 0) {
//         assert(proc != rmi.procid());
//         rmi.remote_call(proc,
//                         &distributed_graph<VertexData, EdgeData>::update_alot2,
//                         blockpushes[i][proc]);
//         blockpushes[i][proc].clear();
//       }
//     }
//   }
// }


// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::push_owned_scope_to_replicas(vertex_id_t vid, 
//                                                                             bool onlymodified, 
//                                                                             bool clearmodified, 
//                                                                             bool async,
//                                                                             bool untracked) {
//   // fast exit if this is not on a boundary
//   if (boundaryscopesset.find(vid) == boundaryscopesset.end()) return;
//   if (0) {
//    if (is_owned(vid)) {
//       vertex_id_t localvid = global2localvid[vid];
//       if (localstore.vertex_modified(localvid)) {
//         localstore.set_vertex_modified(localvid, false);
//         push_owned_vertex_to_replicas(vid, async, untracked);
        
//       }
//       foreach(edge_id_t eid, in_edge_ids(vid)) {
//         if (localstore.edge_modified(eid)) {
//           localstore.set_edge_modified(eid, false);
//           push_owned_edge_to_replicas(eid, async, untracked);
//         }
//       }
//     }
//   }
//   else {
//     if (is_owned(vid)) {
//       push_owned_vertex_to_replicas(vid, async, untracked);
//       foreach(edge_id_t eid, in_edge_ids(vid)) {
//         push_owned_edge_to_replicas(eid, async, untracked);
//       }
//     }
//   }
// }

// template <typename VertexData, typename EdgeData>
// void distributed_graph<VertexData, EdgeData>::wait_for_all_async_pushes() {
//   pending_push_updates.wait();
// }

// #endif


