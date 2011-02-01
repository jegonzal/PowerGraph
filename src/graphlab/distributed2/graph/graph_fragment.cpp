#include <graphlab/distributed2/graph/graph_fragment.hpp>

#include <iostream>
#include <iomanip>
#include <sstream>


#include <graphlab/macros_def.hpp>

using namespace graphlab;
using namespace graphlab::graph_fragment;




description::description() { }

description::    
description(const std::string& base,
                     vertex_id_t id, 
                     vertex_id_t nfragments,
                     vertex_id_t nverts,
                     vertex_id_t nedges) :
  id(id), 
  nfragments(nfragments), 
  nverts(nverts), 
  nedges(nedges) {
  
  // Make the file name
  std::stringstream strm;
  strm << base
       << std::setw(3) << std::setfill('0')
       << id
       << ".bin";
  name = strm.str();

  // Compute the vertex range
  fragment_size = nverts / nfragments;
  assert(fragment_size >= 1);
  fragment_remainder = nverts % nfragments;
  begin_vertex = fragment_size * id + std::min(id, fragment_remainder);
  end_vertex = fragment_size * (id + 1) + std::min(id+1, fragment_remainder);

  assert(begin_vertex < end_vertex);
  // std::cout << "id(" << id << ") : " 
  //           << begin_vertex << " to " << end_vertex 
  //           << " with " << num_local_verts() << " local objects."
  //           << std::endl;
  
  
  }



bool description::
is_local(vertex_id_t vid) const {
  return vid >= begin_vertex && vid < end_vertex;
}

vertex_id_t description::
owning_fragment(vertex_id_t vid) const {   
  if(vid < (fragment_size+1) * fragment_remainder) 
    return vid / (fragment_size+1);
  else 
    return (vid - ((fragment_size+1) * fragment_remainder)) / fragment_size 
      + fragment_remainder;
}

vertex_id_t description::
num_local_verts() const {
  return end_vertex - begin_vertex;
}


void description::
save(graphlab::oarchive& oarc) const {
  oarc << name
       << id
       << nfragments
       << nverts
       << nedges
       << fragment_size
       << fragment_remainder
       << begin_vertex
       << end_vertex;
}

void description::
load(graphlab::iarchive& iarc) {
  iarc >> name
       >> id
       >> nfragments
       >> nverts
       >> nedges
       >> fragment_size
       >> fragment_remainder
       >> begin_vertex
       >> end_vertex;
}

void description::
vids_to_fragmentids(const std::vector<vertex_id_t>& vids,
                    std::vector<vertex_id_t>& fragmentids) const {

  fragmentids.clear();
  fragmentids.resize(vids.size());
  for(vertex_id_t i = 0; i < vids.size(); ++i) {
    fragmentids[i] = owning_fragment(vids[i]);
  }

}

std::ostream& 
operator<<(std::ostream& out, 
           const graphlab::graph_fragment::description& desc) {
  return out << "name:                " << desc.name << std::endl
             << "id:                  " << desc.id << std::endl
             << "nblocks:             " << desc.nfragments << std::endl
             << "nverts:              " << desc.nverts << std::endl
             << "nedges:              " << desc.nedges << std::endl
             << "fragment_size:       " << desc.fragment_size << std::endl
             << "fragment_remainder:  " << desc.fragment_remainder << std::endl
             << "begin_vertex:        " << desc.begin_vertex << std::endl
             << "end_vertex:          " << desc.end_vertex << std::endl
             << "local verts:         " << desc.num_local_verts() << std::endl;
  
}










#include <graphlab/macros_undef.hpp>


