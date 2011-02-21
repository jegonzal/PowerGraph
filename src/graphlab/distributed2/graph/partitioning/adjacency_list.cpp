#include <boost/filesystem.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <graphlab/distributed2/graph/partitioning/adjacency_list.hpp>

#include <graphlab/macros_def.hpp>

using namespace graphlab;


const std::string adjacency_list::elist_suffix = ".elist";
const std::string adjacency_list::vlist_suffix = ".vlist";
const std::string adjacency_list::vdata_suffix = ".vdata";
const std::string adjacency_list::edata_suffix = ".edata";





vertex_id_t adjacency_list::
add_vertex(const vertex_id_t& vid) {
  typedef global2local_type::const_iterator iterator;
  iterator iter = global2local.find(vid);
  // If the vertex has already been added just return the local id
  if(iter != global2local.end()) {
    return iter->second;
  } else {
    vertex_id_t localvid = local_vertices.size();
    local_vertices.push_back(vid);
    in_neighbor_ids.resize(in_neighbor_ids.size() + 1);
    global2local[vid] = localvid;
    return localvid;
  }  
} // end of add vertex


vertex_id_t adjacency_list::
get_local_vid(const vertex_id_t& vid) const {
  typedef global2local_type::const_iterator iterator;
  iterator iter = global2local.find(vid);
  assert(iter != global2local.end());
  return iter->second;
} // end of get local vid


void adjacency_list::
add_edge(const vertex_id_t& source, const vertex_id_t& target,
         const bool require_target_ownership) {
  vertex_id_t target_localvid(0);
  if(require_target_ownership) 
    target_localvid = get_local_vid(target);
  else 
    target_localvid = add_vertex(target);
  in_neighbor_ids.at(target_localvid).push_back(source);
} // end of add edge


void adjacency_list::
load(const std::string& fname) {
  // be sure to have the vlist fname
  std::string vlist_fname = change_suffix(fname, vlist_suffix);
  {
    std::ifstream fin(vlist_fname.c_str());
    while(fin.good()) {
      size_t vid(0);
      fin >> vid;
      if(fin.good())
        add_vertex(vid);
    }
    fin.close();
  }
  // be sure to have the elist fname
  std::string elist_fname = change_suffix(fname, elist_suffix);
  const bool REQUIRE_HAS_TARGET(true);
  {
    std::ifstream fin(elist_fname.c_str());
    while(fin.good()) {
      size_t source(0);
      size_t target(0);
      fin >> source >> target;
      if(fin.good()) 
        add_edge(source, target, REQUIRE_HAS_TARGET);  
    }
    fin.close();
  }
  assert(local_vertices.size() == in_neighbor_ids.size());
  assert(local_vertices.size() == global2local.size());
} // end of load file


void adjacency_list::
save(const std::string& fname, const size_t& id) const {
  assert(local_vertices.size() == in_neighbor_ids.size());
  assert(local_vertices.size() == global2local.size());
  {
    // be sure to have the vlist fname
    std::string vlist_fname = make_fname(fname, id, vlist_suffix);
    std::ofstream fout(vlist_fname.c_str());
    assert(fout.good());
    for(size_t i = 0; i < local_vertices.size(); ++i) {
      fout << local_vertices[i] << '\n';
      assert(fout.good());
    }
    fout.close();
  }
  {
    // be sure to have the elist fname
    std::string elist_fname = make_fname(fname, id, elist_suffix);
    std::ofstream fout(elist_fname.c_str());
    assert(fout.good());
    for(size_t i = 0; i < in_neighbor_ids.size(); ++i) {
      vertex_id_t target( local_vertices[i] );
      for(size_t j = 0; j < in_neighbor_ids[i].size(); ++j) {
        vertex_id_t source(in_neighbor_ids[i][j]);
        fout << source << '\t' << target << '\n';
        assert(fout.good());
      }
    }
    fout.close();
  }
} // end of save file


void adjacency_list::
operator+=(const adjacency_list& other) {
  for(size_t i = 0; i < other.in_neighbor_ids.size(); ++i) {    
    vertex_id_t target = other.local_vertices.at(i);
    for(size_t j = 0; j < other.in_neighbor_ids[i].size(); ++j) {
      vertex_id_t source = other.in_neighbor_ids[i][j];
      add_edge(source, target);
    }
  }
}





 // adjacency_list::adjacency_list(const std::string& base,
//                            vertex_id_t id, 
//                            vertex_id_t nfragments,
//                            vertex_id_t nverts,
//                            vertex_id_t nedges) :
//   id(id), 
//   nfragments(nfragments), 
//   nverts(nverts), 
//   nedges(nedges) {
  
//   // Make the file name
//   {
//     std::stringstream strm;
//     strm << base
//          << std::setw(3) << std::setfill('0')
//          << id
//          << structure_suffix;
//     structure_filename = strm.str();
//   }
//   {
//     std::stringstream strm;
//     strm << base
//          << std::setw(3) << std::setfill('0')
//          << id
//          << vdata_suffix;
//     vertex_data_filename = strm.str();
//   }
//   {
//     std::stringstream strm;
//     strm << base
//          << std::setw(3) << std::setfill('0')
//          << id
//          << edata_suffix;
//     edge_data_filename = strm.str();
//   }
// }


std::string adjacency_list::
make_fname(const std::string& base,
                     const size_t& id,
                     const std::string& suffix) {
  std::stringstream strm;
  strm << base << "_"
       << std::setw(3) << std::setfill('0')
       << id
       << suffix;
  return strm.str();
}



std::string adjacency_list::
change_suffix(const std::string& base,
              const std::string& new_suffix) {
 
  size_t pos = base.rfind('.');
  assert(pos != std::string::npos); 
  std::string new_base = base.substr(0, pos);
  return new_base + new_suffix;

}


void adjacency_list::
list_vlist_files(const std::string& pathname, 
                 std::vector<std::string>& files) {
  namespace fs = boost::filesystem;
  fs::path path(pathname);
  assert(fs::exists(path));
  for(fs::directory_iterator iter( path ), end_iter; 
      iter != end_iter; ++iter) {
    if( ! fs::is_directory(iter->status()) ) {
      std::string filename(iter->path().filename());
      size_t pos = 
        filename.size() >= adjacency_list::vlist_suffix.size()?
        filename.size() - adjacency_list::vlist_suffix.size() : 0;
      std::string ending(filename.substr(pos));
      if(ending == adjacency_list::vlist_suffix) {
        files.push_back(iter->path().filename());
      }
      // size_t period_loc = filename.rfind('.');
      // if(period_loc != std::string::npos) {
      //   std::string ending(filename.substr(period_loc));
      //   if(ending == adjacency_list::alist_suffix) {
      //     files.push_back(iter->path().filename());
      //   }
      // }
    }
  }
  std::sort(files.begin(), files.end());
} // end of list_structure_files





void adjacency_list::
check_local_structures(size_t nverts) const {
  assert(in_neighbor_ids.size() < nverts);
  assert(local_vertices.size() == in_neighbor_ids.size());
  { // Check local vertices is a set
    std::set<vertex_id_t> lvids;
    lvids.insert(local_vertices.begin(), local_vertices.end());
    assert(lvids.size() == local_vertices.size());
  }

  for(size_t i = 0; i < local_vertices.size(); ++i) {    
    vertex_id_t vid(local_vertices[i]);
    assert(vid < nverts);
    {    // Check Global to local map is valid
      typedef global2local_type::const_iterator iterator;
      iterator iter = global2local.find(vid);
      assert(iter != global2local.end());
      assert(iter->second == i);
    }
    {    // Check that neighbors have valid ids
      for(size_t j = 0; j < in_neighbor_ids[i].size(); ++j)
        assert(in_neighbor_ids[i][j] < nverts);
    }
    {     // check that neighbors list is a set
      std::set<vertex_id_t> neighbors;
      neighbors.insert(in_neighbor_ids[i].begin(), 
                       in_neighbor_ids[i].end());
      assert(neighbors.size() == in_neighbor_ids[i].size());
    }
  } 
} // end of check local structures










#include <graphlab/macros_undef.hpp>


