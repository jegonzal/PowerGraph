/* Copyright (c) 2009 Carnegie Mellon University. 
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
 *  Code written by Danny Bickson, CMU
 *  Any changes to the code must include this original license notice in full.
 */


#ifndef MATRIX_LOADER_HPP
#define MATRIX_LOADER_HPP

#include <omp.h>
#include <vector>


#include "matrixmarket/mmio.h"

#include <graphlab.hpp>


#include <graphlab/macros_def.hpp>

struct matrix_descriptor {
  int rows, cols, nonzeros;
  matrix_descriptor() : rows(0), cols(0), nonzeros(0) { }
}; // end of matrix descriptor

template<typename Graph>
struct matrix_entry {
  typedef Graph graph_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  vertex_id_type source, target;
  edge_data_type edata;
  matrix_entry() : source(0), target(0) { }
  matrix_entry(const vertex_id_type& source, 
               const vertex_id_type& target,
               const edge_data_type& edata) :
    source(source), target(target), edata(edata) { }
}; // end of matrix_entry



template<typename Graph>
bool load_matrixmarket(const std::string& fname,
                             matrix_descriptor& desc,
                             std::vector< matrix_entry<Graph> >& test_set) {
  typedef Graph graph_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  typedef matrix_entry<graph_type> matrix_entry_type;

  // Open the file 
  FILE* fptr = fopen(fname.c_str(), "r");
  if(fptr == NULL) {
    logstream(LOG_ERROR) << "Unable to open file " << fname << std::endl;
    return false;
  }
  // read Matrix market header
  MM_typecode matcode;
  if(mm_read_banner(fptr, &matcode)) {
    logstream(LOG_ERROR) << "Unable to read banner" << std::endl;
    return false;
  }
  // Screen header type
  if (mm_is_complex(matcode) || !mm_is_matrix(matcode)) {
    logstream(LOG_ERROR) 
      << "Sorry, this application does not support matrixmarket type: "
      <<  mm_typecode_to_str(matcode) << std::endl;
    return false;
  }
  // load the matrix descriptor
  if(mm_read_mtx_crd_size(fptr, &desc.rows, &desc.cols, &desc.nonzeros)) {
    logstream(LOG_ERROR) << "Error reading dimensions" << std::endl;
  }
  std::cout << "Rows:      " << desc.rows << std::endl
            << "Cols:      " << desc.cols << std::endl
            << "Nonzeros:  " << desc.nonzeros << std::endl;
  std::cout << "Constructing all vertices." << std::endl;
   std::cout << "Adding edges." << std::endl;
  for(size_t i = 0; i < size_t(desc.nonzeros); ++i) {    
    int row = 0, col = 0;  double val = 0;
    if(fscanf(fptr, "%d %d %lg\n", &row, &col, &val) != 3) {
      logstream(LOG_ERROR) 
        << "Error reading file on line: " << i << std::endl;
      return false;
    } --row; --col;
    ASSERT_LT(row, desc.rows);
    ASSERT_LT(col, desc.cols);
    const vertex_id_type source = row;
    const vertex_id_type target = col + desc.rows;
    const edge_data_type edata(val);
    test_set.push_back(matrix_entry_type(source, target, edata));
  } // end of for loop  
  return true;
} // end of load matrixmarket graph


template<typename Graph>
bool load_matrixmarket_graph(const std::string& fname,
                             matrix_descriptor& desc,
                             Graph& graph) {
  typedef Graph graph_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  typedef matrix_entry<graph_type> matrix_entry_type;

  // Open the file 
  FILE* fptr = fopen(fname.c_str(), "r");
  if(fptr == NULL) {
    logstream(LOG_ERROR) << "Unable to open file " << fname << std::endl;
    return false;
  }
  // read Matrix market header
  MM_typecode matcode;
  if(mm_read_banner(fptr, &matcode)) {
    logstream(LOG_ERROR) << "Unable to read banner" << std::endl;
    return false;
  }
  // Screen header type
  if (mm_is_complex(matcode) || !mm_is_matrix(matcode)) {
    logstream(LOG_ERROR) 
      << "Sorry, this application does not support matrixmarket type: "
      <<  mm_typecode_to_str(matcode) << std::endl;
    return false;
  }
  // load the matrix descriptor
  if(mm_read_mtx_crd_size(fptr, &desc.rows, &desc.cols, &desc.nonzeros)) {
    logstream(LOG_ERROR) << "Error reading dimensions" << std::endl;
  }
  std::cout << "Rows:      " << desc.rows << std::endl
            << "Cols:      " << desc.cols << std::endl
            << "Nonzeros:  " << desc.nonzeros << std::endl;
  std::cout << "Constructing all vertices." << std::endl;
  graph.resize(desc.rows + desc.cols);
  std::cout << "Adding edges." << std::endl;
  for(size_t i = 0; i < size_t(desc.nonzeros); ++i) {    
    int row = 0, col = 0;  double val = 0;
    if(fscanf(fptr, "%d %d %lg\n", &row, &col, &val) != 3) {
      logstream(LOG_ERROR) 
        << "Error reading file on line: " << i << std::endl;
      return false;
    } --row; --col;
    ASSERT_LT(row, desc.rows);
    ASSERT_LT(col, desc.cols);
    const vertex_id_type source = row;
    const vertex_id_type target = col + desc.rows;
    const edge_data_type edata(val);
    graph.add_edge(source, target, edata);
  } // end of for loop  
  std::cout << "Graph size:    " << graph.num_edges() << std::endl;
  graph.finalize();
  return true;
} // end of load matrixmarket graph


template<typename Graph>
bool load_tsv_graph(const std::string& fname,
                    matrix_descriptor& desc,
                    Graph& graph) {  
  return false;
} // end of laod tsv graph


template<typename Graph>
bool load_graph(const std::string& fname,
                const std::string& format,
                matrix_descriptor& desc,
                Graph& graph) {
  if(format == "matrixmarket") 
    return load_matrixmarket_graph(fname, desc, graph);
  else if(format == "tsv")
    return load_tsv_graph(fname, desc, graph);
  else std::cout << "Invalid file format!" << std::endl;
  return false;
} // end of load graph

extern bool debug;

template<typename Graph>
void initialize_vertex_data(const size_t nlatent, Graph& graph) {
  typedef typename Graph::vertex_id_type vertex_id_type;
  typedef typename Graph::vertex_data_type vertex_data_type;
#pragma omp parallel for
  for(ssize_t vid = 0; vid < ssize_t(graph.num_vertices()); ++vid) {
    // Randomly initialize the vertex data
    vertex_data_type& vdata = graph.vertex_data(vid);
    vdata.latent.resize(nlatent);
    for(size_t i = 0; i < nlatent; ++i) 
      vdata.latent(i) = debug? 0.1: graphlab::random::gaussian();
  } 
} // end of initialize vertex data









#include <graphlab/macros_undef.hpp>
#endif // end of matrix_loader_hpp
