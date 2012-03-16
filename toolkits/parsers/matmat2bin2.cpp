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


#include <cstdio>
#include <map>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/unordered_map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/unordered_map.hpp>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;
using namespace std;


bool debug = false;
bool quick = true;
std::string datafile;
graphlab::timer mytime;
unsigned long long total_lines = 0;
unsigned long long self_edges = 0;
bool gzip = false;
bool reverse_edges = false;

struct vertex_data {
  string filename;
  vertex_data(){ };
  vertex_data(std::string _filename) : filename(_filename) { }
  void add_self_edge(double val){ };
}; // end of vertex_data

struct edge_data {
   edge_data(double val){ };
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;



struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {
    
   std::string dir = context.get_global<std::string>("PATH");
   std::string outdir = context.get_global<std::string>("OUTPATH");
    //open file
    vertex_data& vdata = context.vertex_data();
    std::ifstream in_file((dir + vdata.filename).c_str(), std::ios::binary);
    logstream(LOG_INFO)<<"Opening input file: " << dir << vdata.filename << std::endl;
    if (!in_file.good()){
       logstream(LOG_ERROR)<<"Failed to open file!"<< std::endl;
       return;
    }
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;
    if (gzip)
      fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);  

    bipartite_graph_descriptor info;
    graph_type graph;
    //get matrix market header size
    load_graph(vdata.filename, "matrixmarket", info, graph, MATRIX_MARKET_3, false, true);

    FILE * deg_file = open_file(vdata.filename + ".nodes", "w");
    FILE * edge_file = open_file(vdata.filename + ".edges", "w");
    FILE * weight_file = open_file(vdata.filename + ".weights", "w");

    //int deg_ptr = 0;
    uint edges_so_far = 0;
    fwrite(&edges_so_far, sizeof(int), 1, deg_file);
    int num_degree_written = 1;  
    
    if (reverse_edges){
        for (uint k=0; k < info.rows-1; k++){
             fwrite(&edges_so_far, sizeof(int), 1, deg_file);
             num_degree_written++;
           }
     }
    
    int last_row = 0, last_col = 0;
    char linebuf[256], buf1[256], buf2[256], buf3[256];
    char saveptr[1024];
    uint line = 1;
    int lines = context.get_global<int>("LINES");
    bool header = true;
 
    while(true){
      fin.getline(linebuf, 128);
      if (fin.eof())
        break;
      if (linebuf[0] == '%'){ //skip matrix market header
         continue;
      }
      else if (header){ //skip matrix market size
          header = false;
          continue;
      }
      char *pch = strtok_r(linebuf," \r\n\t,",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl << " line is: " << linebuf << std::endl;
        return;
       }
      strncpy(buf1, pch, 20);
      pch = strtok_r(NULL, " \r\n\t;,",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
       strncpy(buf2, pch, 20);
        pch = strtok_r(NULL, "\r\n\t ;,",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }

       strncpy(buf3, pch, 40);
        
        uint from, to;
        double val;
        from = atoi(buf1);
        to = atoi(buf2);
        if (reverse_edges) //file is sorted by 2nd column
           assert(to > (uint)last_col);
        else assert(from > (uint)last_row); //file is sorted by 1st column

        val = atof(buf3);
        //matrix market size line- skip
        if (from == (uint)info.rows && to == (uint)info.cols && (size_t)val == info.nonzeros)
           continue;

        assert(from >= 1 && from <= (uint)info.rows);
        assert(to >=1 && to <= (uint)info.cols);
  	from--; to--;

        fwrite(&val, sizeof(double), 1, weight_file);

        if (reverse_edges){
           fwrite(&from, sizeof(int), 1, edge_file);
        }
        else {
           uint abs_dst = to + info.rows;
           fwrite(&abs_dst, sizeof(int), 1, edge_file);
        }

        if (reverse_edges && to > (uint)last_col){
           for (uint k=last_col; k < to; k++){
             fwrite(&edges_so_far, sizeof(int), 1, deg_file);
             num_degree_written++;
           }
        }
        else if (!reverse_edges && from > (uint)last_row){
	   for (uint k=last_row; k < from; k++){
             fwrite(&edges_so_far, sizeof(int), 1, deg_file);
             num_degree_written++;
           }
        }
        edges_so_far++;

        if (debug && line <= 10)
            cout<<"Read line: " << line << " From: " << from << " To: " << to << " val: " << val << " total edges so far: " << edges_so_far << endl;
       
      line++;
      total_lines++;

      last_row = from; last_col = to;
      if (lines && line>=lines)
	 break;

      if (debug && (line % 1000000 == 0))
        logstream(LOG_INFO) << "Parsed line: " << line << endl;
    } 

   logstream(LOG_INFO) <<"Finished parsing total of " << line << " lines in file " << vdata.filename << endl;
   if (reverse_edges) {
     for (int k=last_col; k < info.cols; k++){
       fwrite(&edges_so_far, sizeof(int), 1, deg_file);
       num_degree_written++;
     }
   }
   else {
     for (int k=last_row; k < info.rows+info.cols-1; k++){
       fwrite(&edges_so_far, sizeof(int), 1, deg_file);
       num_degree_written++;
     }
   }
   
    fwrite(&edges_so_far, sizeof(int), 1, deg_file);
    num_degree_written++;
    assert(num_degree_written == info.rows+info.cols+1); 
    assert((size_t)edges_so_far == info.nonzeros);
    
    fclose(deg_file);
    fclose(weight_file);
    fclose(edge_file);

    // close file
    fin.pop(); if (gzip) fin.pop();
    //fout.pop(); fout.pop();
    in_file.close();
    //out_file.close();
  }


};





int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string format = "plain";
  std::string dir = "/mnt/bigbrofs/usr3/bickson/phone_calls/";
  std::string outdir = "/mnt/bigbrofs/usr3/bickson/out_phone_calls/";
  int unittest = 0;
  int lines = 0;
  bool load = false;
  bool save_to_text = false;
  std::string filter = "";

  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("lines", &lines, lines, "limit number of read lines to XX");
  clopts.attach_option("quick", &quick, quick, "quick mode");
  clopts.attach_option("dir", &dir, dir, "path to files");
  clopts.attach_option("outdir", &outdir, outdir, "output directory");
  clopts.attach_option("load", &load, load, "load map from file");
  clopts.attach_option("save_to_text", & save_to_text, save_to_text, 
                       "save output map in text file");
  clopts.attach_option("filter", &filter, filter, "Filter input files starting with prefix.. ");
  clopts.attach_option("gzip", &gzip, gzip, "gzipped input file");
  clopts.attach_option("reverse_edges", &reverse_edges, reverse_edges, "true = matrix market file is sorted by 2nd col. False - matrix market file is sorted by 1st column");
  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab Linear solver library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Gaussian Belief Propagation, "
    << "Jacobi method, Conjugate Gradient" << std::endl;



  // Create a core
  graphlab::core<graph_type, stringzipparser_update> core;
  core.set_options(clopts); // Set the engine options
  core.set_scope_type("vertex");
  mytime.start();

  std::vector<std::string> in_files;
  if (datafile.size() > 0)
     in_files.push_back(datafile);
  else in_files = list_all_files_in_dir(dir, filter);
  assert(in_files.size() >= 1);
  for (int i=0; i< (int)in_files.size(); i++){
      vertex_data data(in_files[i]);
      core.graph().add_vertex(data);
  }

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(stringzipparser_update());
 
  //accumulator acum;
  //core.add_sync("sync", acum, sync_interval);
  core.add_global("LINES", lines); 
  core.add_global("PATH", dir);
  core.add_global("OUTPATH", outdir);


  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;
  std::cout << "Total number of edges: " << self_edges << std::endl;
   return EXIT_SUCCESS;
}



