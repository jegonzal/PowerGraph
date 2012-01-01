/* Copyright (c) 2009 Carnegie Mellon University. 
 * ad
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


#ifndef IO_HPP
#define IO_HPP



#include <omp.h>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>



#include "mmio.h"
#include "mathlayer.hpp"
#include "types.hpp"
#include <graphlab.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include <graphlab/macros_def.hpp>


extern bool debug;


enum matrix_market_parser{
   MATRIX_MARKET_3 = 1,
   MATRIX_MARKET_4 = 2,
   MATRIX_MARKET_5 = 3,
   MATRIX_MARKET_6 = 4
};



/*
 * open a file and verify open success
 */
FILE * open_file(const char * name, const char * mode, bool optional = false){
  FILE * f = fopen(name, mode);
  if (f == NULL && !optional){
      perror("fopen failed");
      logstream(LOG_FATAL) <<" Failed to open file" << name << std::endl;
   }
  return f;
}


/*
 * list all existing files inside a specificed directory
 */
std::vector<std::string> list_all_files_in_dir(const std::string & dir, const std::string &prefix){
  namespace fs = boost::filesystem;  
  std::vector<std::string> ret;

  fs::path dir_path(dir);
  fs::directory_iterator end_iter;
  if ( fs::exists(dir_path) && fs::is_directory(dir_path)) {
    for( fs::directory_iterator dir_iter(dir_path) ; dir_iter != end_iter ; ++dir_iter) {
      if (fs::is_regular_file(dir_iter->status()) ) {
        std::string filename;
#define BOOST_FILESYSTEM_VERSION 2
#if BOOST_FILESYSTEM_VERSION == 2 
        filename = dir_iter->leaf();
#else
        filename = dir_iter->path().filename().string();
#endif
        if (prefix.size() > 0 && !boost::starts_with(filename,prefix)) 
          continue;
        ret.push_back(filename);
      }
    }
  }
  return ret;
}

/*
 * extract the output from node data ito a vector of values
 */
template<typename graph_type>
vec  fill_output(graph_type * g, bipartite_graph_descriptor & matrix_info, int field_type){
  typedef typename graph_type::vertex_data_type vertex_data_type;

  vec out = zeros(matrix_info.num_nodes(false));
  for (int i = matrix_info.get_start_node(false); i < matrix_info.get_end_node(false); i++){
    out[i] = g->vertex_data(i).get_output(field_type);
  }
  return out;
}

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
                       bipartite_graph_descriptor& desc,
                       std::vector< matrix_entry<Graph> >& test_set) {
  typedef Graph graph_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  typedef matrix_entry<graph_type> matrix_entry_type;

  // Open the file 
  FILE* fptr = open_file(fname.c_str(), "r");
  
  // read Matrix market header
  MM_typecode matcode;
  if(mm_read_banner(fptr, &matcode)) {
    logstream(LOG_FATAL) << "Unable to read banner" << std::endl;
  }
  // Screen header type
  if (mm_is_complex(matcode) || !mm_is_matrix(matcode)) {
    logstream(LOG_FATAL) 
      << "Sorry, this application does not support matrixmarket type: "
      <<  mm_typecode_to_str(matcode) << std::endl;
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
      logstream(LOG_FATAL) 
        << "Error reading file on line: " << i << std::endl;
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
bool load_matrixmarket_cpp_graph(const std::string& fname,
                             bipartite_graph_descriptor& desc,
                             Graph& graph,
		             bool gzip = false,
			     int parse_type = MATRIX_MARKET_3){ 
  typedef Graph graph_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  typedef matrix_entry<graph_type> matrix_entry_type;

  // Open the file 
  logstream(LOG_INFO) << "Reading matrix market file: " << fname << std::endl;
  //FILE* fptr = open_file(fname.c_str(), "r");
  std::ifstream in_file(fname.c_str(), std::ios::binary);
  boost::iostreams::filtering_stream<boost::iostreams::input> fin;
  if (gzip)
    fin.push(boost::iostreams::gzip_decompressor());
  fin.push(in_file); 
    
  MM_typecode matcode;
  int tmprows, tmpcols, tmpnnz;

    // read Matrix market header
    if(mm_read_cpp_banner(fin, &matcode)) {
      logstream(LOG_FATAL) << "Unable to read banner" << std::endl;
    }
    // Screen header type
    if (mm_is_complex(matcode) || !mm_is_matrix(matcode)) {
      logstream(LOG_FATAL) 
      << "Sorry, this application does not support matrixmarket type: "
      <<  mm_typecode_to_str(matcode) << std::endl;
     return false;
    }
    // load the matrix descriptor
    if(mm_read_cpp_mtx_crd_size(fin, &tmprows, &tmpcols, &tmpnnz)) {
      logstream(LOG_FATAL) << "Error reading dimensions" << std::endl;
    }

  //update dataset size from file only with empty structure
  if (desc.rows == 0 && desc.cols == 0 && desc.nonzeros == 0){
     desc.rows = tmprows; desc.cols = tmpcols; desc.nonzeros = tmpnnz;
  }

  std::cout << "Rows:      " << desc.rows << std::endl
            << "Cols:      " << desc.cols << std::endl
            << "Nonzeros:  " << desc.nonzeros << std::endl;
  std::cout << "Constructing all vertices." << std::endl;

  if ((int)graph.num_vertices() < desc.total())
    graph.resize(desc.total());
  bool is_square = desc.is_square();

  char line[MM_MAX_LINE_LENGTH];

  std::cout << "Adding edges." << std::endl;
  for(size_t i = 0; i < size_t(desc.nonzeros); ++i) {    
    int row = 0, col = 0, inttime = 0, intdate = 0;  
    double val = 0;
    unsigned long long time = 0;
	

    if (fin.eof()){
       if (i < (size_t)(desc.nonzeros - 1))
            logstream(LOG_WARNING) <<"File " << fname << " ended after " << i << " expected lines: " << desc.nonzeros << std::endl;
       break;
    }
    fin.getline(line, MM_MAX_LINE_LENGTH);
    

    //regular matrix market format. [from] [to] [val]
    if (parse_type == MATRIX_MARKET_3){ 
        if(sscanf(line, "%d %d %lg\n", &row, &col, &val) != 3) {
        logstream(LOG_ERROR) 
          << "Error reading file on line: " << i << std::endl;
        return false;
      }
    }
     //extended matrix market format. [from] [to] [val from->to] [val to->from] [ignored] [ignored]
    else if (parse_type == MATRIX_MARKET_4){ 
        if(sscanf(line, "%d %d %llu %lg\n", &row, &col, &time, &val) != 4) {
        logstream(LOG_ERROR) 
          << "Error reading file on line: " << i << std::endl;
        return false;
      }
    }
     //extended matrix market format. [from] [to] [val from->to] [val to->from] [ignored] [ignored]
    else if (parse_type == MATRIX_MARKET_5){ 
        if(sscanf(line, "%d %d %d %d %lg\n", &row, &col, &intdate, &inttime, &val) != 5) {
        logstream(LOG_ERROR) 
          << "Error reading file on line: " << i << std::endl;
        return false;
      }
     //extended matrix market format. [from] [to] [val from->to] [val to->from] [ignored] [ignored]
    }  else if (parse_type == MATRIX_MARKET_6){
      double val2, zero, zero1;
      if(sscanf(line, "%d %d %lg %lg %lg %lg\n", &row, &col, &val, &val2, &zero, &zero1) != 6) {
        logstream(LOG_FATAL) 
          << "Error reading file " << fname << " on line: " << i << std::endl;
        return false;
      }
      val += val2; //sum up to values to have a single undirected link
    }
    else assert(false);

    --row; --col;

    ASSERT_LT(row, desc.rows);
    ASSERT_LT(col, desc.cols);
    ASSERT_GE(row, 0);
    ASSERT_GE(col, 0);
    const vertex_id_type source = row;
    const vertex_id_type target = col + (is_square ? 0 : desc.rows);
    const edge_data_type edata(val);

    if (debug && desc.nonzeros < 100)
      logstream(LOG_INFO)<<"Adding an edge: " << source << "->" << target << " with val: " << std::endl;

    if(is_square && source == target) 
      graph.vertex_data(source).add_self_edge(val);
    else {
     graph.add_edge(source, target, edata); 
      if (mm_is_symmetric(matcode))
        graph.add_edge(target, source, edata);
    }
  } // end of for loop  
  std::cout << "Graph size:    " << graph.num_edges() << std::endl;
  //graph.finalize();
  
   if (gzip)
     fin.pop();
   fin.pop();
   in_file.close();
   return true;
} // end of load matrixmarket graph


template<typename Graph>
bool load_matrixmarket_graph(const std::string& fname,
                             bipartite_graph_descriptor& desc,
                             Graph& graph,
			     int parse_type = MATRIX_MARKET_3){ 
  typedef Graph graph_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  typedef matrix_entry<graph_type> matrix_entry_type;

  // Open the file 
  logstream(LOG_INFO) << "Reading matrix market file: " << fname << std::endl;
  FILE* fptr = open_file(fname.c_str(), "r");
  
  // read Matrix market header
  MM_typecode matcode;
  if(mm_read_banner(fptr, &matcode)) {
    logstream(LOG_FATAL) << "Unable to read banner" << std::endl;
  }
  // Screen header type
  if (mm_is_complex(matcode) || !mm_is_matrix(matcode)) {
    logstream(LOG_FATAL) 
      << "Sorry, this application does not support matrixmarket type: "
      <<  mm_typecode_to_str(matcode) << std::endl;
    return false;
  }
  // load the matrix descriptor
  if(mm_read_mtx_crd_size(fptr, &desc.rows, &desc.cols, &desc.nonzeros)) {
    logstream(LOG_FATAL) << "Error reading dimensions" << std::endl;
  }
  std::cout << "Rows:      " << desc.rows << std::endl
            << "Cols:      " << desc.cols << std::endl
            << "Nonzeros:  " << desc.nonzeros << std::endl;
  std::cout << "Constructing all vertices." << std::endl;
  graph.resize(desc.total());
  bool is_square = desc.is_square();

  std::cout << "Adding edges." << std::endl;
  for(size_t i = 0; i < size_t(desc.nonzeros); ++i) {    
    int row = 0, col = 0;  
    double val = 0;

    //regular matrix market format. [from] [to] [val]
    if (parse_type == MATRIX_MARKET_3){ 
      if(fscanf(fptr, "%d %d %lg\n", &row, &col, &val) != 3) {
        logstream(LOG_ERROR) 
          << "Error reading file on line: " << i << std::endl;
        return false;
      }
     //extended matrix market format. [from] [to] [val from->to] [val to->from] [ignored] [ignored]
    } else if (parse_type == MATRIX_MARKET_6){
      double val2, zero, zero1;
      if(fscanf(fptr, "%d %d %lg %lg %lg %lg\n", &row, &col, &val, &val2, &zero, &zero1) != 6) {
        logstream(LOG_FATAL) 
          << "Error reading file " << fname << " on line: " << i << std::endl;
        return false;
      }
      val += val2; //sum up to values to have a single undirected link
    }
    else assert(false);

    --row; --col;

    ASSERT_LT(row, desc.rows);
    ASSERT_LT(col, desc.cols);
    ASSERT_GE(row, 0);
    ASSERT_GE(col, 0);
    const vertex_id_type source = row;
    const vertex_id_type target = col + (is_square ? 0 : desc.rows);
    const edge_data_type edata(val);

    if (debug && desc.nonzeros < 100)
      logstream(LOG_INFO)<<"Adding an edge: " << source << "->" << target << " with val: " << std::endl;

    if(is_square && source == target) 
      graph.vertex_data(source).add_self_edge(val);
    else {
     graph.add_edge(source, target, edata); 
      if (mm_is_symmetric(matcode))
        graph.add_edge(target, source, edata);
    }
  } // end of for loop  
  std::cout << "Graph size:    " << graph.num_edges() << std::endl;
  //graph.finalize();
  return true;
} // end of load matrixmarket graph


template<typename Graph>
bool load_tsv_graph(const std::string& fname,
                    bipartite_graph_descriptor& desc,
                    Graph& graph) {  
  return false;
} // end of laod tsv graph


template<typename Graph>
bool load_graph(const std::string& fname,
                const std::string& format,
                bipartite_graph_descriptor& desc,
                Graph& graph, 
	        int format_type = MATRIX_MARKET_3) {

  if(format == "matrixmarket") 
    return load_matrixmarket_graph(fname, desc, graph, format_type);
  else if(format == "tsv")
    return load_tsv_graph(fname, desc, graph);
  else std::cout << "Invalid file format!" << std::endl;
  return false;
} // end of load graph

template<typename Graph>
bool load_cpp_graph(const std::string& fname,
                const std::string& format,
                bipartite_graph_descriptor& desc,
                Graph& graph, 
	        bool gzip = false,
	        int format_type = MATRIX_MARKET_3) {

  if(format == "matrixmarket") 
    return load_matrixmarket_cpp_graph(fname, desc, graph, gzip, format_type);
  else std::cout << "Invalid file format!" << std::endl;
  return false;
} // end of load graph


template <typename graph_type>
void load_matrix_market_vector(const std::string & filename, const bipartite_graph_descriptor & desc, graph_type & g, int type, bool optional_field)
{
    typedef typename graph_type::vertex_data_type vertex_data;
    
    int ret_code;
    MM_typecode matcode;
    int M, N, nz;   
    int i;

    logstream(LOG_INFO) <<"Going to read matrix market vector from input file: " << filename << std::endl;
  
    FILE * f = open_file(filename.c_str(), "r", optional_field);
    //if optional file not found return
    if (f== NULL && optional_field){
       return;
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        logstream(LOG_ERROR) << "Could not process Matrix Market banner." << std::endl;
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        logstream(LOG_ERROR) << "sorry, this application does not support " << std::endl << 
          "Market Market type: " << mm_typecode_to_str(matcode) << std::endl;
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0){
       logstream(LOG_ERROR) << "failed to read matrix market cardinality size " << std::endl; 
       exit(1);
    }


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    int row,col; 
    double val;

    for (i=0; i<nz; i++)
    {
        int rc = fscanf(f, "%d %d %lg\n", &row, &col, &val);
        if (rc != 3){
	  logstream(LOG_FATAL) << "Failed reading input file: " << filename << "Problm at data row " << i << " (not including header and comment lines)" << std::endl;

        }
        row--;  /* adjust from 1-based to 0-based */
        col--;
        //some users have gibrish in text file - better check both I and J are >=0 as well
        assert(row >=0 && row< M);
        assert(col == 0);
        //set observation value
        vertex_data & vdata = g.vertex_data(row);
        vdata.set_val(val, type);
    }
    fclose(f);

}


template<typename Graph>
void load_vector(const std::string& fname,
                   const std::string& format,
                   const bipartite_graph_descriptor& desc,
                   Graph& graph, 
		   int type,
		   bool optional_field) {

  if (format == "matrixmarket"){
     load_matrix_market_vector(fname, desc, graph, type, optional_field);
     return;
  }
  else assert(false); //TODO other formats

}

inline void write_row(int row, int col, double val, FILE * f, bool issparse){
    if (issparse)
      fprintf(f, "%d %d %10.3g\n", row, col, val);
    else fprintf(f, "%10.3g ", val);
}

inline void write_row(int row, int col, int val, FILE * f, bool issparse){
    if (issparse)
      fprintf(f, "%d %d %d\n", row, col, val);
    else fprintf(f, "%d ", val);
}

template<typename T>
inline void set_typecode(MM_typecode & matcore);

template<>
inline void set_typecode<vec>(MM_typecode & matcode){
   mm_set_real(&matcode);
}

template<>
inline void set_typecode<ivec>(MM_typecode & matcode){
  mm_set_integer(&matcode);
}


template<typename vec>
void save_matrix_market_format_vector(const std::string datafile, const vec & output, bool issparse)
{
    MM_typecode matcode;                        
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);

    if (issparse)
       mm_set_sparse(&matcode);
    else mm_set_dense(&matcode);

    set_typecode<vec>(matcode);

    FILE * f = fopen(datafile.c_str(),"w");
    assert(f != NULL);
    mm_write_banner(f, matcode); 
    mm_write_mtx_crd_size(f, output.size(), 1, output.size());

    for (int j=0; j<(int)output.size(); j++){
      write_row(j+1, 1, output[j], f, issparse);
      if (!issparse) 
        fprintf(f, "\n");
    }

    fclose(f);
}

template<typename mat>
void save_matrix_market_format_matrix(const std::string datafile, const mat & output, bool issparse)
{
    MM_typecode matcode;                        
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);

    if (issparse)
      mm_set_sparse(&matcode);
    else mm_set_dense(&matcode);

    set_typecode<vec>(matcode);

    FILE * f = fopen(datafile.c_str(),"w");
    assert(f != NULL);
    mm_write_banner(f, matcode); 
    mm_write_mtx_crd_size(f, output.size(), 1, output.size());

    for (int j=0; j<(int)output.rows(); j++)
      for (int i=0; i< (int)output.cols(); i++){
         write_row(j+1, i+1, get_val(output, i, j), f, issparse);
         if (!issparse && (i == (int)output.cols() - 1))
           fprintf(f, "\n");
      }

    fclose(f);
}



//read a vector from file and return an array
inline double * read_vec(FILE * f, size_t len){
  double * vec = new double[len];
  assert(vec != NULL);
  int rc = fread(vec, sizeof(double), len, f);
  assert(rc == (int)len);
  return vec;
}


//write an output vector to file
template<typename T>
inline void write_vec(const FILE * f, const int len, const T * array){
  int total = 0;
  assert(f != NULL && array != NULL);
  while(total < len){
     int rc = fwrite(array, sizeof(T), len, (FILE*)f);
    if (rc <= 0){
      perror("fwrite");
      logstream(LOG_FATAL) << "Failed writing array" << std::endl; 
    }
    total += rc;
  }
  assert(total == len);
}

//write an output vector to file
inline void write_vec(const FILE * f, const int len, const int * array){
  assert(f != NULL && array != NULL);
  int rc = fwrite(array, sizeof(int), len, (FILE*)f);
  assert(rc == len);
}

inline void write_output_vector_binary(const std::string & datafile, const uint* output, int size){

   FILE * f = open_file(datafile.c_str(), "w");
   std::cout<<"Writing result to file: "<<datafile<<std::endl;
   write_vec(f, size, output);
   fclose(f);
}



template<typename vec>
inline void write_output_vector_binary(const std::string & datafile, const vec& output){

   FILE * f = open_file(datafile.c_str(), "w");
   std::cout<<"Writing result to file: "<<datafile<<std::endl;
   std::cout<<"You can read the file in Matlab using the load_c_gl.m matlab script"<<std::endl;
   write_vec(f, output.size(), &output[0]);
   fclose(f);
}

template<typename mat>
inline void write_output_matrix_binary(const std::string & datafile, const mat& output){

   FILE * f = open_file(datafile.c_str(), "w");
   std::cout<<"Writing result to file: "<<datafile<<std::endl;
   std::cout<<"You can read the file in Matlab using the load_c_gl.m matlab script"<<std::endl;
   write_vec(f, output.size(), data(output));
   fclose(f);
}


template<typename vec>
inline void write_output_vector(const std::string & datafile, const std::string & format, const vec& output, bool issparse){

  if (format == "binary")
    write_output_vector_binary(datafile, output);
  else if (format == "matrixmarket")
    save_matrix_market_format_vector(datafile, output,issparse); 
  else assert(false);
}

template<typename mat>
inline void write_output_matrix(const std::string & datafile, const std::string & format, const mat& output, bool isparse){

  if (format == "binary")
    write_output_matrix_binary(datafile, output);
  else if (format == "matrixmarket")
    save_matrix_market_format_matrix(datafile, output); 
  else assert(false);
}


//read matrix size from a binary file
FILE * load_matrix_metadata(const char * filename, bipartite_graph_descriptor & desc){
   printf("Loading %s\n", filename);
   FILE * f = open_file(filename, "r", false);

   int rc = fread(&desc.rows, sizeof(desc.rows), 1, f);
   assert(rc == 1);
   rc = fread(&desc.cols, sizeof(desc.cols), 1, f);
   assert(rc == 1);
   return f;
}


template<typename Graph>
bool load_binary_graph(const std::string& fname,
                             bipartite_graph_descriptor& desc,
                             Graph& graph) {
  typedef Graph graph_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  typedef matrix_entry<graph_type> matrix_entry_type;

  // Open the file 
  logstream(LOG_INFO) << "Reading matrix market file: " << fname << std::endl;
  FILE* fptr = open_file(fname.c_str(), "r");
  
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
  graph.resize(desc.total());
  bool is_square = desc.is_square();

  std::cout << "Adding edges." << std::endl;
  for(size_t i = 0; i < size_t(desc.nonzeros); ++i) {    
    int row = 0, col = 0;  
    double val = 0;
    if(fscanf(fptr, "%d %d %lg\n", &row, &col, &val) != 3) {
      logstream(LOG_FATAL) 
        << "Error reading file " << fname << " on line: " << i << std::endl;
    } --row; --col;
    ASSERT_LT(row, desc.rows);
    ASSERT_LT(col, desc.cols);
    ASSERT_GE(row, 0);
    ASSERT_GE(col, 0);
    const vertex_id_type source = row;
    const vertex_id_type target = col + (is_square ? 0 : desc.rows);
    const edge_data_type edata(val);

    if (debug && desc.nonzeros < 100)
      logstream(LOG_INFO)<<"Adding an edge: " << source << "->" << target << " with val: " << std::endl;

    if(is_square && source == target) 
      graph.vertex_data(source).add_self_edge(val);
    else
      graph.add_edge(source, target, edata);
  } // end of for loop  
  std::cout << "Graph size:    " << graph.num_edges() << std::endl;
  //graph.finalize();
  return true;
} // end of load matrixmarket graph

uint array_from_file(std::string filename, uint *& array){
          struct stat sb;
          int fd = open (filename.c_str(), O_RDONLY);
          if (fd == -1) {
                  perror ("open");
                  logstream(LOG_FATAL) << "Failed to open input file: " << filename << std::endl;
          }

          if (fstat (fd, &sb) == -1) {
                  perror ("fstat");
                  logstream(LOG_FATAL) << "Failed to get size of  input file: " << filename << std::endl;
          }

          if (!S_ISREG (sb.st_mode)) {
                  logstream(LOG_FATAL) << "Input file: " << filename 
              << " is not a regular file and can not be mapped " << std::endl;
          }
	  close(fd);
 
	  int toread = sb.st_size/4; 
          array = new uint[toread];
          int total = 0;
	  FILE * f = fopen(filename.c_str(), "r");
          if (f == NULL){
	     perror("fopen");
             logstream(LOG_FATAL) << "Failed to open input file: " << filename << std::endl;
          }
         
          while(total < toread){
	     int rc = fread(array+total, sizeof(uint), toread-total,f);
	     if (rc < 0 ){
	       perror("fread");
               logstream(LOG_FATAL) << "Failed to read from input file: " << filename << std::endl;
	     }
	     total += rc; 
          }
          return sb.st_size;
}


uint mmap_from_file(std::string filename, uint *& array){
          struct stat sb;
          int fd = open (filename.c_str(), O_RDONLY);
          if (fd == -1) {
                  perror ("open");
                  logstream(LOG_FATAL) << "Failed to open input file: " << filename << std::endl;
          }

          if (fstat (fd, &sb) == -1) {
                  perror ("fstat");
                  logstream(LOG_FATAL) << "Failed to get size of  input file: " << filename << std::endl;
          }

          if (!S_ISREG (sb.st_mode)) {
                  logstream(LOG_FATAL) << "Input file: " << filename 
              << " is not a regular file and can not be mapped " << std::endl;
          }

          array = (uint*)mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
          if (array == MAP_FAILED) {
                  perror ("mmap");
                  logstream(LOG_FATAL) << "Failed to map input file: " << filename << std::endl;
          }
          return sb.st_size;
}


#include <graphlab/macros_undef.hpp>
#endif // end of matrix_loader_hpp
