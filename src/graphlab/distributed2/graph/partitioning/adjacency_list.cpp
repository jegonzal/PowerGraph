#include <boost/filesystem.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <graphlab/distributed2/graph/partitioning/adjacency_list.hpp>

#include <graphlab/macros_def.hpp>

using namespace graphlab;


const std::string adjacency_list::alist_suffix = ".structure";
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


void adjacency_list::
add_edge(const vertex_id_t& source, const vertex_id_t& target) {
  vertex_id_t target_localvid = add_vertex(target);
  in_neighbor_ids.at(target_localvid).push_back(source);
} // end of add edge


void adjacency_list::
load(const std::string& fname) {
  std::ifstream fin(fname.c_str());
  assert(fin.good());
  while(fin.good()) {
    size_t source(0);
    size_t target(0);
    fin >> source >> target;
    if(fin.good()) add_edge(source, target);  
  }
  fin.close();
} // end of load file


void adjacency_list::
save(const std::string& fname) const {
  std::ofstream fout(fname.c_str());
  for(size_t i = 0; i < in_neighbor_ids.size(); ++i) {    
    vertex_id_t target = local_vertices.at(i);
    for(size_t j = 0; j < in_neighbor_ids[i].size(); ++j) {
      vertex_id_t source = in_neighbor_ids[i][j];
      fout << source << '\t' << target << '\n';
    }
    assert(fout.good());
  }
  fout.close();  
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
make_alist_fname(const std::string& base,
                 const size_t& id) {
  return make_fname(base, id, alist_suffix);
}

std::string adjacency_list::
make_vdata_fname(const std::string& base,
                 const size_t& id) {
  return make_fname(base, id, vdata_suffix);
}

std::string adjacency_list::
make_edata_fname(const std::string& base,
                 const size_t& id) {
  return make_fname(base, id, edata_suffix);
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
list_adjacency_files(const std::string& pathname, 
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
        if(ending == adjacency_list::alist_suffix) {
          files.push_back(iter->path().filename());
        }
      }
    }
  }
  std::sort(files.begin(), files.end());
} // end of list_structure_files

















#include <graphlab/macros_undef.hpp>


