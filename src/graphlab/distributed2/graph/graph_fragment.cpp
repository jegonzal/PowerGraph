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
                     vertex_id_t nblocks,
                     vertex_id_t nverts,
                     vertex_id_t nedges) :
  id(id), 
  nblocks(nblocks), 
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
  block_size = nverts / nblocks;
  assert(block_size >= 1);
  block_remainder = nverts % nblocks;
  begin_vertex = block_size * id + std::min(id, block_remainder);
  end_vertex = block_size * (id + 1) + std::min(id+1, block_remainder);

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
owning_block(vertex_id_t vid) const {   
  if(vid < (block_size+1) * block_remainder) 
    return vid / (block_size+1);
  else 
    return (vid - ((block_size+1) * block_remainder)) / block_size 
      + block_remainder;
}

vertex_id_t description::
num_local_verts() const {
  return end_vertex - begin_vertex;
}


void description::
save(graphlab::oarchive& oarc) const {
  oarc << name
       << id
       << nblocks
       << nverts
       << nedges
       << block_size
       << block_remainder
       << begin_vertex
       << end_vertex;
}

void description::
load(graphlab::iarchive& iarc) {
  iarc >> name
       >> id
       >> nblocks
       >> nverts
       >> nedges
       >> block_size
       >> block_remainder
       >> begin_vertex
       >> end_vertex;
}








#include <graphlab/macros_undef.hpp>


