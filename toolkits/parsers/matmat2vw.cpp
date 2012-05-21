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
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;
using namespace std;

bool debug = false;
string datafile;
bool info = false;
size_t total_lines = 0;
bool gzip = false;

int get_day(time_t pt);

struct vertex_data {
  string filename;
  vertex_data() { }
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data

struct edge_data {
};



struct singlerating{
   int user;
   int item;
   double rating;
   //int time;
   singlerating(int _user, int _item, double _rating/*, int _time*/) : user(_user), item(_item), rating(_rating)/*, time(_time)*/{ };
};


typedef graphlab::graph<vertex_data, edge_data> graph_type;

struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {

    
   vertex_data& vdata = context.vertex_data();
   gzip_in_file fin(vdata.filename +  (info ? ".data" : ""), gzip);
   //gzip_out_file fout(vdata.filename + ".out", gzip); 
   FILE * pfile = open_file((vdata.filename + ".out").c_str(), "w");
    char linebuf[24000];
    char saveptr[1024];
    int last_from = -1, last_to = -1;
    std::vector<singlerating> multiple_ratings;
    int rows, cols;
    size_t nz;
 
    MM_typecode matcode;
    if (!info){
      int rc = mm_read_cpp_banner(fin.get_sp(), &matcode);
      if (rc)
        logstream(LOG_FATAL)<<"Failed to read mm banner" << endl;
      rc = mm_read_cpp_mtx_crd_size(fin.get_sp(), &rows, &cols, &nz);
      if (rc)
        logstream(LOG_FATAL)<<"Failed to read mm banner" << endl;
    }
    else {
     gzip_in_file info_file(vdata.filename + ".info", gzip);
      int rc = mm_read_cpp_banner(info_file.get_sp(), &matcode);
      if (rc)
        logstream(LOG_FATAL)<<"Failed to read mm banner" << endl;
      rc = mm_read_cpp_mtx_crd_size(info_file.get_sp(), &rows, &cols, &nz);
      if (rc)
        logstream(LOG_FATAL)<<"Failed to read mm banner" << endl;
 
    }

    while(true){
      fin.get_sp().getline(linebuf, 24000);
      if (fin.get_sp().eof())
        break;

      char *pch = strtok_r(linebuf," \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << "[" << linebuf << "]" << std::endl;
        return;
       }
      int from = atoi(pch);
      if (from <= 0 || from > rows){
         logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << " first column is out of range " << from << std::endl;
         return;
      }
      pch = strtok_r(NULL," \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << "[" << linebuf << "]" << std::endl;
        return;
       }
      int to = atoi(pch);
      if (to <= 0 || to > cols){
         logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << " second column is out of range" << to << std::endl;
         return;
      }
      pch = strtok_r(NULL," \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << "[" << linebuf << "]" << std::endl;
        return;
       }
      double val = atof(pch);
      total_lines++;
      if (debug && (total_lines % 50000 == 0))
        logstream(LOG_INFO) << "Parsed line: " << total_lines << std::endl;

      singlerating thisrating(from, to, val);
      if ((last_from > -1) && (last_from != from)){
        assert(multiple_ratings.size() > 0 && multiple_ratings.size() < 1000000);
        //fout.get_sp() << "0 | ";
        fprintf(pfile, "0 | ");
        for (uint i=0; i< multiple_ratings.size(); i++){
          //fout.get_sp()<<multiple_ratings[i].item<<":"<<multiple_ratings[i].rating<<" ";
          fprintf(pfile, "%d:%lg ", multiple_ratings[i].item, multiple_ratings[i].rating);
        } 
        //fout.get_sp()<<endl;
        fprintf(pfile, "\n");
        //fout.get_sp().strict_sync();
        multiple_ratings.clear();
        multiple_ratings.push_back(thisrating);
      }
      else multiple_ratings.push_back(thisrating);
   
      last_from = from;
      last_to = to;
    } 

   logstream(LOG_INFO) <<"Finished parsing total of " << total_lines << " lines in file " << vdata.filename << endl;
   if (total_lines != nz)
     logstream(LOG_FATAL)<<"Expected a total of " << nz << " lines, while in practice tehre where: " << total_lines << endl;
   assert(multiple_ratings.size() > 0);
        fprintf(pfile, "0 | ");
        for (uint i=0; i< multiple_ratings.size(); i++){
          //fout.get_sp()<<multiple_ratings[i].item<<":"<<multiple_ratings[i].rating<<" ";
          fprintf(pfile, "%d:%lg ", multiple_ratings[i].item, multiple_ratings[i].rating);
        } 
        //fout.get_sp()<<endl;
        fprintf(pfile, "\n");
        //fout.get_sp().strict_sync();
        multiple_ratings.clear();
        

  }


};



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Parsers Library");
  clopts.attach_option("data", &datafile, datafile,
                       "training data input file");
  clopts.add_positional("data");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("gzip", &gzip, gzip, "Gzipped input file?");
  clopts.attach_option("info", &info, info, "matrix market file is given in .info file");
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

  vertex_data data(datafile);
  core.graph().add_vertex(0, data);
  core.schedule_all(stringzipparser_update());
  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;
}



