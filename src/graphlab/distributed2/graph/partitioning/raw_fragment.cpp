#include <boost/filesystem.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <graphlab/distributed2/graph/partitioning/raw_fragment.hpp>

#include <graphlab/macros_def.hpp>

using namespace graphlab;


const std::string raw_fragment::structure_suffix = ".structure";
const std::string raw_fragment::vdata_suffix = ".vdata";
const std::string raw_fragment::edata_suffix = ".edata";


raw_fragment::raw_fragment(const std::string& base,
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

  neighbor_ids.resize(num_local_verts);
  in_neighbor_ids.resize(num_local_verts);


}



bool raw_fragment::
is_local(vertex_id_t vid) const {
  return vid >= begin_vertex && vid < end_vertex;
}

vertex_id_t raw_fragment::
local_id(vertex_id_t gvid) const {
  assert(is_local(gvid));
  return gvid - begin_vertex;  
}

vertex_id_t raw_fragment::
owning_fragment(vertex_id_t vid) const {   
  if(vid < (fragment_size+1) * fragment_remainder) 
    return vid / (fragment_size+1);
  else 
    return (vid - ((fragment_size+1) * fragment_remainder)) / fragment_size 
      + fragment_remainder;
}



void raw_fragment::
vids_to_fragmentids(const std::vector<vertex_id_t>& vids,
                    std::vector<vertex_id_t>& fragmentids) const {

  fragmentids.clear();
  fragmentids.resize(vids.size());
  for(vertex_id_t i = 0; i < vids.size(); ++i) {
    fragmentids[i] = owning_fragment(vids[i]);
  }

}








bool raw_fragment::
add_edge(const vertex_id_t source, const vertex_id_t target) {
  assert(is_local(source) || is_local(target));
  bool is_inbound(is_local(target));
  if(is_inbound) {
    is_inbound = true;
    vertex_id_t local_target(local_id(target));
    neighbor_ids[local_target].push_back(source);
    in_neighbor_ids[local_target].push_back(source);
  } 
  if(is_local(source)) {
    assert(is_local(source));
    vertex_id_t local_source(local_id(source));
    neighbor_ids[local_source].push_back(target);
  } 
  return is_inbound;
} // end of add edge








void raw_fragment::finalize_adjacency_lists() {
  for(size_t i = 0; i < neighbor_ids.size(); ++i) {
    std::vector<vertex_id_t>& vec(neighbor_ids[i]);
    std::sort(vec.begin(), vec.end());
    typedef std::vector<vertex_id_t>::iterator iterator;
    iterator new_end = std::unique(vec.begin(), vec.end());
    vec.resize(std::distance(vec.begin(), new_end));
  }
}





void raw_fragment::
list_structure_files(const std::string& pathname, 
                     std::vector<std::string>& files) {
  namespace fs = boost::filesystem;
  fs::path path(pathname);
  assert(fs::exists(path));
  for(fs::directory_iterator iter( path ), end_iter; 
      iter != end_iter; ++iter) {
    if( ! fs::is_directory(iter->status()) ) {
      std::string filename(iter->path().filename());
      size_t period_loc = filename.rfind('.');
      if(period_loc != std::string::npos) {
        std::string ending(filename.substr(period_loc));
        if(ending == raw_fragment::structure_suffix) {
          files.push_back(iter->path().filename());
        }
      }
    }
  }
  std::sort(files.begin(), files.end());
} // end of list_structure_files






void raw_fragment::
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
       << num_local_edges
    // Actual structural payload
       << neighbor_ids
       << in_neighbor_ids;
}






void raw_fragment::
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
       >> num_local_edges
       >> neighbor_ids
       >> in_neighbor_ids;
}


void raw_fragment::
partial_load(graphlab::iarchive& iarc) {
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
graphlab::operator<<(std::ostream& out,  const raw_fragment& frag) {
  return out << "structure filename:  " << frag.structure_filename << std::endl
             << "edge data filename:  " << frag.edge_data_filename << std::endl
             << "vdata filename:      " << frag.vertex_data_filename << std::endl
             << "id:                  " << frag.id << std::endl
             << "nblocks:             " << frag.nfragments << std::endl
             << "nverts:              " << frag.nverts << std::endl
             << "nedges:              " << frag.nedges << std::endl
             << "fragment_size:       " << frag.fragment_size << std::endl
             << "fragment_remainder:  " << frag.fragment_remainder << std::endl
             << "begin_vertex:        " << frag.begin_vertex << std::endl
             << "end_vertex:          " << frag.end_vertex << std::endl
             << "local verts:         " << frag.num_local_verts << std::endl
             << "local edges:         " << frag.num_local_edges << std::endl;  
}




















#include <graphlab/macros_undef.hpp>


