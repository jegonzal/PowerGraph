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
bool quick = false;
boost::unordered_map<string,uint> hash2nodeid;
std::string datafile;

uint conseq_id;
graphlab::mutex mymutex;
graphlab::timer mytime;
unsigned long long total_lines = 0;
unsigned long long self_edges = 0;
unsigned long long filtered_out = 0;
int min_time = 0, max_time = 24*3600; //filter out records based on time range
bool negate_time = false;
int total_graphs_done = 0;

struct vertex_data {
  string filename;
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data


struct edge_data {
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
unsigned long int datestr2uint64(const std::string & data, int & dateret, int & timeret, int threadid);



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
   int mythreadid = thread::thread_id();

    //open file
    vertex_data& vdata = context.vertex_data();
    gzip_in_file fin((dir + vdata.filename));
    std::string out_filename = outdir + vdata.filename + boost::lexical_cast<std::string>(min_time) + "-" + 
       boost::lexical_cast<std::string>(max_time) + (negate_time? "neg" : "") + ".out.gz"; 

    gzip_out_file fout(out_filename);
   
    MM_typecode out_typecode;
    mm_clear_typecode(&out_typecode);
    mm_set_integer(&out_typecode); 
    mm_set_dense(&out_typecode); 
    mm_set_matrix(&out_typecode);
    mm_write_cpp_banner(fout.get_sp(), out_typecode);
    mm_write_cpp_mtx_crd_size(fout.get_sp(), 987654321, 987654321, 987654322);


    char linebuf[256], buf1[256], buf2[256], buf3[256];
    char saveptr[1024];
    int duration;
    int line = 1;
    int lines = context.get_global<int>("LINES");
    int dateret, timeret;
   
    while(true){
      fin.get_sp().getline(linebuf, 128);
      if (fin.get_sp().eof())
        break;

       char *pch = strtok_r(linebuf," \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
        return;
       }
      strncpy(buf1, pch, 20);
      pch = strtok_r(NULL, " \r\n\t;",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
       strncpy(buf2, pch, 20);
      if (!quick){
        pch = strtok_r(NULL, " ",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
         strncpy(buf3, pch, 6);
        buf3[6] = ' ';
        pch = strtok_r(NULL, " ",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
         strncpy(buf3+7,pch,6);
        pch = strtok_r(NULL, " \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
         duration = atoi(pch);
      if (duration < 0){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
         datestr2uint64(std::string(buf3), timeret, dateret, mythreadid);
        uint from, to;
        find_ids(from, to, buf1, buf2);
        if (debug && line <= 10)
            cout<<"Read line: " << line << " From: " << from << " To: " << to << " timeret: " << timeret << " date: " << dateret << " val: " << duration << endl;
         
        if (from == to)
          self_edges++;
        else if (!negate_time && ((timeret < min_time*3600) || (timeret > max_time*3600)))
            filtered_out++;
         else if (negate_time && ((timeret >= min_time*3600) && (timeret <= max_time*3600)))
            filtered_out++;
         else // fout << from << " " << to << " " << dateret << " " << timeret << " " << duration << endl;
            { 
            fout.get_sp() << from << " " << to << endl;
            total_lines++;
         }
      }
      else {
        uint from,to;
        find_ids(from, to, buf1, buf2);
        if (from != to){
          fout.get_sp() << from << " " << to << endl;
          //fout << to << " " << from << endl;
        }
        else self_edges++;
      }

      //fin.read(buf1,1); //go over \n
      line++;
      if (lines && line>=lines)
	 break;

      if (line % 5000000 == 0)
        logstream(LOG_INFO) << mytime.current_time() << ") Parsed line: " << line << " chosen lines " << total_lines <<  " filtered out: " << filtered_out << " slef edges: " << self_edges << std::endl;
    } 

   total_graphs_done++;
   logstream(LOG_INFO) << mytime.current_time() << ") Finished parsing total of " << line << " lines in file " << vdata.filename <<
	                 "total lines " << total_lines << " total graphs done: " << total_graphs_done << endl;

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
  std::string outdir = "/usr2/bickson/filtered.hours/";
  std::string filter;
  int unittest = 0;
  int lines = 0;
  int nodes = 121408373;

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
  clopts.attach_option("min_time", &min_time, min_time, "filter out records < min_time");
  clopts.attach_option("max_time", &max_time, max_time, "filter out records > max_time");
  clopts.attach_option("filter", &filter, filter, "select files starting with prefi [filter]");
  clopts.attach_option("negate_time", &negate_time, negate_time, "invert time selection");
  clopts.attach_option("outdir", &outdir, outdir, "output directory");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;

  assert(min_time < max_time);

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
      if (in_files[i].find(".gz") != string::npos){
        vertex_data data(in_files[i]);
        core.graph().add_vertex(data);
     }
  }

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(stringzipparser_update());
 
  //accumulator acum;
  //core.add_sync("sync", acum, sync_interval);
  core.add_global("LINES", lines); 
  core.add_global("PATH", dir);
  core.add_global("OUTPATH", outdir);
    mytime.start();
    //hash2nodeid.rehash(nodes);
   load_map_from_file(hash2nodeid, outdir + ".map"); 

   double runtime= core.start();
  std::cout << "Finished in " << runtime << std::endl;


  logstream(LOG_INFO)<<"Wrote total edges: " << total_lines << " in time: " << mytime.current_time() << std::endl;
  logstream(LOG_INFO)<<"Total edges filtered out: " << filtered_out << std::endl; 

   return EXIT_SUCCESS;
}



