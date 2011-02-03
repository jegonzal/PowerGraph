
#include <graphlab/distributed2/graph/graph_fragment.hpp>

#include <iostream>
#include <iomanip>
#include <sstream>


#include <graphlab/macros_def.hpp>

using namespace graphlab;
using namespace graphlab::graph_fragment;


const std::string graphlab::graph_fragment::
structure_suffix(".structure");
const std::string graphlab::graph_fragment::
vdata_suffix(".vdata");
const std::string graphlab::graph_fragment::
edata_suffix(".edata");

file_description::file_description() { }

file_description::    
file_description(const std::string& base,
                     vertex_id_t id, 
                     vertex_id_t nfragments,
                     vertex_id_t nverts,
                     vertex_id_t nedges) :
  id(id), 
  nfragments(nfragments), 
  nverts(nverts), 
  nedges(nedges),
  num_local_verts(0),
  num_local_edges(0) {
  
  // Make the file name
  {
    std::stringstream strm;
    strm << base
         << std::setw(3) << std::setfill('0')
         << id
         << structure_suffix;
    structure_filename = strm.str();
  }
  {
    std::stringstream strm;
    strm << base
         << std::setw(3) << std::setfill('0')
         << id
         << vdata_suffix;
    vertex_data_filename = strm.str();
  }
  {
    std::stringstream strm;
    strm << base
         << std::setw(3) << std::setfill('0')
         << id
         << edata_suffix;
    edge_data_filename = strm.str();
  }

  // Compute the vertex range
  fragment_size = nverts / nfragments;
  assert(fragment_size >= 1);
  fragment_remainder = nverts % nfragments;
  begin_vertex = fragment_size * id + std::min(id, fragment_remainder);
  end_vertex = fragment_size * (id + 1) + std::min(id+1, fragment_remainder);

  assert(end_vertex >= begin_vertex);
  num_local_verts = end_vertex - begin_vertex;

  assert(begin_vertex < end_vertex);
  // std::cout << "id(" << id << ") : " 
  //           << begin_vertex << " to " << end_vertex 
  //           << " with " << num_local_verts() << " local objects."
  //           << std::endl;
  
  
}



bool file_description::
is_local(vertex_id_t vid) const {
  return vid >= begin_vertex && vid < end_vertex;
}

vertex_id_t file_description::
local_id(vertex_id_t gvid) const {
  assert(is_local(gvid));
  return gvid - begin_vertex;  
}

vertex_id_t file_description::
owning_fragment(vertex_id_t vid) const {   
  if(vid < (fragment_size+1) * fragment_remainder) 
    return vid / (fragment_size+1);
  else 
    return (vid - ((fragment_size+1) * fragment_remainder)) / fragment_size 
      + fragment_remainder;
}



void file_description::
vids_to_fragmentids(const std::vector<vertex_id_t>& vids,
                    std::vector<vertex_id_t>& fragmentids) const {

  fragmentids.clear();
  fragmentids.resize(vids.size());
  for(vertex_id_t i = 0; i < vids.size(); ++i) {
    fragmentids[i] = owning_fragment(vids[i]);
  }

}



void file_description::
save(graphlab::oarchive& oarc) const {
  oarc << structure_filename
       << edge_data_filename
       << vertex_data_filename
       << id
       << nfragments
       << nverts
       << nedges
       << fragment_size
       << fragment_remainder
       << begin_vertex
       << end_vertex
       << num_local_verts
       << num_local_edges;
}



void file_description::
load(graphlab::iarchive& iarc) {
  iarc >> structure_filename
       >> edge_data_filename
       >> vertex_data_filename
       >> id
       >> nfragments
       >> nverts
       >> nedges
       >> fragment_size
       >> fragment_remainder
       >> begin_vertex
       >> end_vertex
       >> num_local_verts
       >> num_local_edges;
}



std::ostream& 
operator<<(std::ostream& out, 
           const graphlab::graph_fragment::file_description& desc) {
  return out << "structure filename:  " << desc.structure_filename << std::endl
             << "edge data filename:  " << desc.edge_data_filename << std::endl
             << "vdata filename:      " << desc.vertex_data_filename << std::endl
             << "id:                  " << desc.id << std::endl
             << "nblocks:             " << desc.nfragments << std::endl
             << "nverts:              " << desc.nverts << std::endl
             << "nedges:              " << desc.nedges << std::endl
             << "fragment_size:       " << desc.fragment_size << std::endl
             << "fragment_remainder:  " << desc.fragment_remainder << std::endl
             << "begin_vertex:        " << desc.begin_vertex << std::endl
             << "end_vertex:          " << desc.end_vertex << std::endl
             << "local verts:         " << desc.num_local_verts << std::endl
             << "local edges:         " << desc.num_local_edges << std::endl;  
}


structure_description::structure_description() { }

structure_description::
structure_description(const file_description& desc) :
  desc(desc), 
  neighbor_ids(desc.num_local_verts),
  in_neighbor_ids(desc.num_local_verts),
  in_edge_ids(desc.num_local_verts) { }

bool structure_description::
add_edge(const vertex_id_t source, const vertex_id_t target) {
  assert(desc.is_local(source) || desc.is_local(target));
  bool is_inbound(desc.is_local(target));
  if(is_inbound) {
    is_inbound = true;
    vertex_id_t local_target(desc.local_id(target));
    neighbor_ids[local_target].push_back(source);
    in_neighbor_ids[local_target].push_back(source);
    edge_id_t local_eid(desc.num_local_edges++);
    in_edge_ids[local_target].push_back(local_eid);

  } 
  if(desc.is_local(source)) {
    assert(desc.is_local(source));
    vertex_id_t local_source(desc.local_id(source));
    neighbor_ids[local_source].push_back(target);
  } 
  return is_inbound;
} // end of add edge



void structure_description::
save(graphlab::oarchive& oarc) const {
  oarc << desc 
       << neighbor_ids
       << in_neighbor_ids 
       << in_edge_ids;

}

void structure_description::
load(graphlab::iarchive& iarc)  {
  iarc >> desc
       >> neighbor_ids
       >> in_neighbor_ids 
       >> in_edge_ids;
}














#include <graphlab/macros_undef.hpp>


