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


#include <boost/functional/hash.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


#include <graphlab/graph/disk_graph.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/graph/mr_disk_graph_construction.hpp>

#include <graphlab/util/fs_util.hpp>



#include "dist_pagerank.hpp"



graphlab::atomic<size_t> ctr;

typedef graphlab::igraph_constructor<vertex_data, edge_data>
igraph_constructor_type;

class graph_constructor : 
  public igraph_constructor_type {

public: 
  typedef igraph_constructor_type base_type;
  typedef base_type::vertex_id_type vertex_id_type;

private:

  size_t procid;
  size_t nprocs;
  std::vector<std::string> links_fnames;
  std::vector<std::string> urls_fnames;


  boost::hash<vertex_id_type> hash_function;
 
public:
  graph_constructor(const std::vector<std::string>& links_fnames, 
                    const std::vector<std::string>& urls_fnames) : 
    links_fnames(links_fnames), urls_fnames(urls_fnames) { }
  

  iteration_method begin(size_t procid, size_t nprocs) {
    // Chunk the file into blocks
    this->procid = procid;
    this->nprocs = nprocs;
    return CallBack;
  }
  
  uint16_t vertex_to_atomid(vertex_id_type vid, uint16_t numatoms)  {
    return (hash_function(vid) * 179424673) % numatoms;
  }


  void load_vertices() {
    size_t ctr = 0;
    namespace bios = boost::iostreams;
    for(size_t fileid = procid; fileid < urls_fnames.size(); fileid += nprocs) {
      std::cout << "Opening url file: " << urls_fnames[fileid] << std::endl;
      // open the file
      std::ifstream in_file(urls_fnames[fileid].c_str(), 
                            std::ios_base::in | std::ios_base::binary);
      assert(in_file.good());
      bios::filtering_stream<bios::input> fin;  
      fin.push(bios::gzip_decompressor());
      fin.push(in_file);
      assert(fin.good());
      // Loop over the contents
      std::string url;
      while(fin.good()) {
        // Load a vertex
        size_t source_id(-1);
        try { fin >> source_id >> url; } catch ( ... ) { 
          std::cout << "Error: on reading url." << std::endl;
          break;
        }
        if(!fin.good()) break;
        if(source_id == size_t(-1)) {
          std::cout << "Read NULL LINK" << std::endl;  break;
        }
        const vertex_data vdata;
        const vertex_color_type color = 0;
        add_vertex(source_id, vdata, color);
        ++ctr;
        if (ctr % 1000000 == 0) {
          std::cout << source_id << ": " << url  << std::endl;
          std::cout << "Added vdata for " << ctr << " vertices" << std::endl;
        }
      } // end of loop over file
      std::cout << "Finished file: " << urls_fnames[fileid] << std::endl;
    } // end of loop over files
  } // end load vertices

  void load_edges() {
    namespace bios = boost::iostreams;
    size_t ctr = 0;
    for(size_t fileid = procid; fileid < links_fnames.size(); fileid += nprocs) {
      std::cout << "Opening link file: " << links_fnames[fileid] << std::endl;
      // open the file
      std::ifstream in_file(links_fnames[fileid].c_str(), 
                            std::ios_base::in | std::ios_base::binary);
      assert(in_file.good());
      bios::filtering_stream<bios::input> fin;  
      fin.push(bios::gzip_decompressor());
      fin.push(in_file);
      assert(fin.good());
      // Loop over the contents
      while(fin.good()) {
        // Load a vertex
        size_t source_id(-1);
        size_t nlinks = 0;
        try { fin >> source_id >> nlinks; } catch ( ... ) { 
          std::cout << "Error: on read sourceid nlinks." << std::endl;
          break;
        }
        if(!fin.good()) break;
        if(source_id == size_t(-1)) {
          std::cout << "Read NULL LINK" << std::endl;  break;
        }
        // add neighbors
        for(size_t i = 0; i < nlinks; ++i) {
          size_t target_id(-1);
          try { fin >> target_id; } catch ( ... ) {
            std::cout << "Error reading neighbor" << std::endl;
            break;
          }
          if(target_id == size_t(-1)) break;
          // add edge
          const edge_data edata(1/double(nlinks));
          add_edge(std::make_pair(source_id, target_id), edata);  
        } // end of loop over neighbors
        ++ctr;
        if (ctr % 1000000 == 0) {
          std::cout << "Added edata for " << ctr << " vertices: " 
                    << source_id << std::endl;
        }
      } // end of loop over file
      std::cout << "Finished file: " << links_fnames[fileid] << std::endl;
    } // end of loop over files
  } // end of generate callback

  void generate_callback() {
    load_vertices();
    load_edges();
  } // end of generate callback
};





  
int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  // create distributed control
  graphlab::distributed_control dc(param);

  const std::string links_directory = std::string(argv[1]) + "/";
  const std::string urls_directory = std::string(argv[2]) + "/";
  const std::string tmp_directory = std::string(argv[3]) + "/";
  const std::string target_directory = std::string(argv[4]) + "/";
 
  std::cout << "links_directory:  " << links_directory << std::endl 
            << "urls directory:   " << urls_directory << std::endl
            << "tmp directory:    " << tmp_directory << std::endl
            << "target directory: " << target_directory << std::endl;

  // Load the links
  std::vector<std::string> links_fnames;
  std::vector<std::string> urls_fnames;
  graphlab::fs_util::list_files_with_suffix(links_directory, ".gz", links_fnames);
  for (size_t i = 0;i < links_fnames.size(); ++i) { 
    links_fnames[i] = links_directory + links_fnames[i];
    std::cout << "links file: " << links_fnames[i] << std::endl;
  }
  graphlab::fs_util::list_files_with_suffix(urls_directory, ".gz", urls_fnames);
  for (size_t i = 0;i < urls_fnames.size(); ++i) 
    urls_fnames[i] = urls_directory + urls_fnames[i]; 
 
  //-----------------------------------------------------------------------
  // create an instance of a graph cosntructor
  graph_constructor constructor(links_fnames, urls_fnames);
  std::cout << "Starting MR construction..." << std::endl;
  graphlab::timer ti;
  ti.start();
  // call the mr_disk_graph_construction function
  const size_t natoms = 512;
  const size_t max_per_node = 1;
  const std::string output_basename = "altavista";
  graphlab::
    mr_disk_graph_construction<graph_constructor, vertex_data, edge_data>
    (dc, 
     constructor, 
     max_per_node,
     output_basename, 
     natoms,
     graphlab::disk_graph_atom_type::WRITE_ONLY_ATOM,
     tmp_directory,
     target_directory);
  std::cout << "Completed in " << ti.current_time() << " s" << std::endl;
  //-----------------------------------------------------------------------
  graphlab::mpi_tools::finalize();
}
