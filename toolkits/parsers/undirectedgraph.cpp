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
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;
using namespace std;


bool debug = false;
bool quick = true;
boost::unordered_map<string,uint> hash2nodeid;
std::string datafile;
//atomic<unsigned int> conseq_id;
uint conseq_id;
graphlab::mutex mymutex;
graphlab::timer mytime;
unsigned long long total_lines = 0;
unsigned long long self_edges = 0;

struct vertex_data {
  string filename;
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data

std::set<uint> * out_edges;
std::set<uint> * in_edges;

struct edge_data {
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
unsigned long int datestr2uint64(const std::string & data, int & dateret, int & timeret);



void assign_id(uint & outval, const string &name){

  boost::unordered_map<string,uint>::iterator it = hash2nodeid.find(name);
  if (it != hash2nodeid.end()){
     outval = it->second;
     return;
  }
  logstream(LOG_ERROR)<<"Did not find map entry: " << name << std::endl;
  assert(false);
}



void find_ids(uint & from, uint & to, const string &buf1, const string& buf2){

   assign_id(from, buf1);
   assign_id(to, buf2);
   assert(from > 0 && to > 0);
}


/***
* Line format is: PnLaCsEnqei atslBvPNusB 050803 235959 590 
*/

struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {
    
   std::string dir = context.get_global<std::string>("PATH");
   std::string outdir = context.get_global<std::string>("OUTPATH");

    //open file
    vertex_data& vdata = context.vertex_data();
    std::ifstream in_file((dir + vdata.filename).c_str(), std::ios::binary);
    logstream(LOG_INFO)<<"Opening input file: " << dir << vdata.filename << std::endl;
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;
    fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);  

    std::ofstream out_file(std::string(outdir + vdata.filename + ".out.gz").c_str(), std::ios::binary);
    logstream(LOG_INFO)<<"Opening output file " << outdir << vdata.filename << ".out.gz" << std::endl;
    boost::iostreams::filtering_stream<boost::iostreams::output> fout;
    fout.push(boost::iostreams::gzip_compressor());
    fout.push(out_file);
   
    MM_typecode out_typecode;
    mm_clear_typecode(&out_typecode);
    mm_set_integer(&out_typecode); 
    mm_set_dense(&out_typecode); 
    mm_set_matrix(&out_typecode);
    mm_write_cpp_banner(fout, out_typecode);
    mm_write_cpp_mtx_crd_size(fout, 987654321, 987654321, 987654322);


    char linebuf[256], buf1[256], buf2[256], buf3[256];
    char saveptr[1024];
    int duration;
    int line = 1;
    int lines = context.get_global<int>("LINES");
    int dateret, timeret;
   
    while(true){
      fin.getline(linebuf, 128);
      if (fin.eof())
        break;

       char *pch = strtok_r(linebuf," \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
        return;
       }
      strncpy(buf1, pch, 18);
      pch = strtok_r(NULL, " \r\n\t;",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
       strncpy(buf2, pch, 18);
      if (!quick){
        pch = strtok_r(NULL, " ",(char**)&saveptr);
        strncpy(buf3, pch, 6);
        buf3[6] = ' ';
        pch = strtok_r(NULL, " ",(char**)&saveptr);
        strncpy(buf3+7,pch,6);
        pch = strtok_r(NULL, " \r\n\t",(char**)&saveptr);
        duration = atoi(pch);
        datestr2uint64(std::string(buf3), dateret, timeret);
        uint from, to;
        find_ids(from, to, buf1, buf2);
        if (debug && line <= 10)
            cout<<"Read line: " << line << " From: " << from << " To: " << to << " timeret: " << dateret << " time: " << timeret << " val: " << duration << endl;
         fout << from << " " << to << " " << dateret << " " << timeret << " " << duration << endl;
      }
      else {
        uint from,to;
        find_ids(from, to, buf1, buf2);
        if (from != to){
          in_edges[from].insert(to);
          out_edges[to].insert(from);
        }
        else self_edges++;
      }

      //fin.read(buf1,1); //go over \n
      line++;
      total_lines++;
      if (lines && line>=lines)
	 break;

      if (debug && (line % 50000 == 0))
        logstream(LOG_INFO) << "Parsed line: " << line << " map size is: " << hash2nodeid.size() << std::endl;
          if (hash2nodeid.size() % 500000 == 0)
        logstream(LOG_INFO) << "Hash map size: " << hash2nodeid.size() << " at time: " << mytime.current_time() << " edges: " << total_lines << std::endl;
    } 

   logstream(LOG_INFO) <<"Finished parsing total of " << line << " lines in file " << vdata.filename << endl <<
	                 "total map size: " << hash2nodeid.size() << endl;

    // close file
    fin.pop(); fin.pop();
    fout.pop(); fout.pop();
    in_file.close();
    out_file.close();
  }


};


/*
class accumulator :
  public graphlab::iaccumulator<graph_type, stringparser_update, accumulator> {
private:
  real_type real_norm, relative_norm;
public:
  accumulator() {}
  void operator()(icontext_type& context) {
  }
  void operator+=(const accumulator& other) { 
  }
  void finalize(iglobal_context_type& context) {
  }
}; // end of  accumulator
*/




int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string format = "plain";
  std::string dir = "/mnt/bigbrofs/usr0/bickson/phone_calls/";
  std::string outdir = "/mnt/bigbrofs/usr0/bickson/out_phone_calls/";
  int unittest = 0;
  int lines = 0;
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
  //unit testing
  if (unittest == 1){
  }

  
  std::vector<std::string> in_files;
  if (datafile.size() > 0)
     in_files.push_back(datafile); 
  else in_files = list_all_files_in_dir(dir);
  assert(in_files.size() >= 1);
  for (int i=0; i< in_files.size(); i++){
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

   logstream(LOG_INFO)<<"Reading hash map from file" << std::endl;
    std::ifstream in_file((outdir + ".map.gz").c_str(), std::ios::binary);
    logstream(LOG_INFO)<<"Opening input file: " << outdir << ".map.gz" << std::endl;
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;
    fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);  

   int line = 0;
    char linebuf[128];
    char saveptr[128], buf1[128], buf2[128];
    while(true){
      fin.getline(linebuf, 128);
      if (fin.eof()){
        logstream(LOG_INFO) << "File ended after " << line << " lines " << std::endl;
        break;
     }
    line++;
      char *pch = strtok_r(linebuf," ",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing imap file: " << ":" << line <<std::endl;
        return EXIT_FAILURE;
       }
      strncpy(buf1, pch, strlen(pch)+1);
      pch = strtok_r(NULL, " \n",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: "  << ":" << line <<std::endl;
         return EXIT_FAILURE;
       }
       strncpy(buf2, pch, strlen(pch)+1);

      hash2nodeid[std::string(buf1)] = boost::lexical_cast<uint>(buf2);
     } 
   logstream(LOG_INFO)<<"Read total of " << hash2nodeid.size() << " entries" << std::endl;

   
   in_edges = new std::set<uint>[hash2nodeid.size()+2];
   out_edges = new std::set<uint>[hash2nodeid.size()+2];


   double runtime= core.start();
  std::cout << "Finished in " << runtime << std::endl;

     std::ofstream out_file(std::string(outdir + ".out.gz").c_str(), std::ios::binary);
    logstream(LOG_INFO)<<"Opening output file " << outdir << ".out.gz" << std::endl;
    boost::iostreams::filtering_stream<boost::iostreams::output> fout;
    fout.push(boost::iostreams::gzip_compressor());
    fout.push(out_file);
   
  {
    fout << in_edges;
    fout << out_edges;
   }

  logstream(LOG_INFO)<<"Wrote total edges: " << total_lines << std::endl;
 

  //vec ret = fill_output(&core.graph(), matrix_info, JACOBI_X);

  //write_output_vector(datafile + "x.out", format, ret);


   
   return EXIT_SUCCESS;
}



