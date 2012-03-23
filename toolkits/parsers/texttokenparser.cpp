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
bool gzip = false;
boost::unordered_map<string,uint> string2nodeid;
boost::unordered_map<uint,string> nodeid2hash;
std::string datafile;
uint conseq_id;
graphlab::mutex mymutex;
graphlab::timer mytime;
unsigned long long total_lines = 0;
unsigned long long self_edges = 0;

struct vertex_data {
  string filename;
  vertex_data() { }
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data

struct edge_data {
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
unsigned long int datestr2uint64(const std::string & data, int & dateret, int & timeret, int thread_id);



void assign_id(uint & outval, const string &name, const int line, const string &filename){

  boost::unordered_map<string,uint>::iterator it = string2nodeid.find(name);
  if (it != string2nodeid.end()){
     outval = it->second;
     return;
  }
  mymutex.lock();
  outval = string2nodeid[name];
  if (outval == 0){
     string2nodeid[name] = ++conseq_id;
     outval = conseq_id;
     nodeid2hash[outval] = name;
  }
  mymutex.unlock();
}



void find_ids(uint & to, const string& buf2, const int line, const string &filename){
   assign_id(to, buf2, line, filename);
}


/***
* Line format is:
1::gift card
2::
3::
4::might have art deco roots but the exaggerated basket weave pattern */
struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {
    
   std::string dir = context.get_global<std::string>("PATH");
   std::string outdir = context.get_global<std::string>("OUTPATH");
   //int mythreadid = thread::thread_id();
    //open file
   vertex_data& vdata = context.vertex_data();
   gzip_in_file fin((dir + vdata.filename), gzip);
   gzip_out_file fout((outdir + vdata.filename + ".out"), gzip);
    /*MM_typecode out_typecode;
    mm_clear_typecode(&out_typecode);
    mm_set_integer(&out_typecode); 
    mm_set_dense(&out_typecode); 
    mm_set_matrix(&out_typecode);
    mm_write_cpp_banner(fout, out_typecode);
    mm_write_cpp_mtx_crd_size(fout, 987654321, 987654321, 987654322);
*/

    char linebuf[24000], buf1[256], buf2[256];
    char saveptr[1024];
    uint line = 1;
    uint lines = context.get_global<uint>("LINES");
   
    while(true){
      fin.get_sp().getline(linebuf, 24000);
      if (fin.get_sp().eof())
        break;

      char *pch = strtok_r(linebuf," \r\n\t:",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line << "[" << linebuf << "]" << std::endl;
        return;
       }
      strncpy(buf1, pch, 20);
      int from = atoi(buf1);
      if (from <= 0){
         logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line << " document ID is zero or less: " << from << std::endl;
         return;
      }
  
      while(pch != NULL){
         pch = strtok_r(NULL, " \r\n\t:",(char**)&saveptr);
         if (pch != NULL && strlen(pch) > 0){ 
           strncpy(buf2, pch, 20);
           uint to; 
           find_ids(to, buf2, line, vdata.filename);
           if (debug && line <= 10)
              cout<<"Read line: " << line << " From: " << from << " To: " << to << " string: " << buf2 << endl;

           fout.get_sp() << from << " " << to << " 1" << endl;
          }
      }

      line++;
      total_lines++;
      if (lines && line>=lines)
	 break;

      if (debug && (line % 50000 == 0))
        logstream(LOG_INFO) << "Parsed line: " << line << " map size is: " << string2nodeid.size() << std::endl;
      if (string2nodeid.size() % 500000 == 0)
        logstream(LOG_INFO) << "Hash map size: " << string2nodeid.size() << " at time: " << mytime.current_time() << " edges: " << total_lines << std::endl;
    } 

   logstream(LOG_INFO) <<"Finished parsing total of " << line << " lines in file " << vdata.filename << endl <<
	                 "total map size: " << string2nodeid.size() << endl;

  }


};



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string format = "plain";
  std::string dir = "/mnt/bigbrofs/usr3/bickson/phone_calls/";
  std::string outdir = "/mnt/bigbrofs/usr3/bickson/out_phone_calls/";
  // int unittest = 0;
  uint lines = 0;
  bool load = false;
  bool save_to_text = false;
  std::string filter = "";

  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("lines", &lines, lines, "limit number of read lines to XX");
  clopts.attach_option("dir", &dir, dir, "path to files");
  clopts.attach_option("outdir", &outdir, outdir, "output directory");
  clopts.attach_option("load", &load, load, "load map from file");
  clopts.attach_option("save_to_text", & save_to_text, save_to_text, 
                       "save output map in text file");
  clopts.attach_option("filter", &filter, filter, "Filter input files starting with prefix.. ");
  clopts.attach_option("gzip", &gzip, gzip, "Gzipped input file?");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab parsers library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented parsers are: Call data records, document tokens "<< std::endl;

  // Create a core
  graphlab::core<graph_type, stringzipparser_update> core;
  core.set_options(clopts); // Set the engine options
  core.set_scope_type("none");
  mytime.start();

  std::vector<std::string> in_files;
  if (datafile.size() > 0)
     in_files.push_back(datafile);
  else in_files = list_all_files_in_dir(dir, filter);
  assert(in_files.size() >= 1);
  for (int i=0; i< (int)in_files.size(); i++){
    vertex_data data(in_files[i]);
    core.graph().add_vertex(vertex_id_type(i), data);
  }

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(stringzipparser_update());
 
  core.add_global("LINES", lines); 
  core.add_global("PATH", dir);
  core.add_global("OUTPATH", outdir);

  if (load){
    mytime.start();
    logstream(LOG_INFO)<<"Opening input file " << outdir << ".map" << std::endl;
    load_map_from_file(string2nodeid, outdir + ".map");
    logstream(LOG_INFO)<<"Finished reading input file in " << mytime.current_time() << std::endl;
  }

  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;
  std::cout << "Total number of edges: " << self_edges << std::endl;

    if (save_to_text){
     save_map_to_text_file(string2nodeid, outdir + ".map.text.gz", gzip);
     save_map_to_text_file(nodeid2hash, outdir + ".reverse.map.text.gz", gzip);
   }

    save_map_to_file(string2nodeid, outdir + ".map");
    save_map_to_file(nodeid2hash, outdir + ".reverse.map");
   return EXIT_SUCCESS;
}



