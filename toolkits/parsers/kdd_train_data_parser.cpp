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

int nodes = 2421057;
int num_edges = 39752419;
int split_time = 1321891200; //time to split test data
int split_training_time = 1320595199;
int ratingA = 19349608;
int ratingB = 15561329;
int split_day_of_year = 310;
int level = 2;
bool debug = false;
string datafile;
string example_submission; 
unsigned long long total_lines = 0;
bool gzip = false;

int get_day(time_t pt);

struct vertex_data {
  string filename;
  vertex_data() { }
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data

struct edge_data {
};

struct vertex_data2 {
  double A_ii;
  double value;
  vertex_data2(): A_ii(1) { } 
  void add_self_edge(double value) { A_ii = value; }
  void set_val(double value, int field_type) { 
  }  
  double get_output(int field_type){ return -1; }
}; // end of vertex_data


struct edge_data2 {
  int rating;
  int time;
  edge_data2(int _rating, int _time) : rating(_rating),time(_time) { }
  int get_field(int pos){ return pos == 0 ? rating : time; }
  void set_field(int pos, int val){ if (pos == 0) rating = val; else time = val; }
};

struct singlerating{
   int user;
   int item;
   int rating;
   int time;

   singlerating(int _user, int _item, int _rating, int _time) : user(_user), item(_item), rating(_rating), time(_time){ };
};


typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::graph<vertex_data2, edge_data2> graph_type2;
typedef graph_type2::edge_type edge_type;

void add_single_edge(int from, int to, edge_data2 & edge, graph_type2 & pgraph, int & examples_added){

  edge_type found = pgraph.find(from - 1, to+nodes-1);
  if (!found.empty())
    return;
  pgraph.add_edge(from - 1, to+nodes-1, edge); 
  examples_added++;
}

struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {

typedef graphlab::graph<vertex_data2, edge_data2>::edge_list_type edge_list;
    
   vertex_data& vdata = context.vertex_data();
   gzip_in_file fin(vdata.filename, gzip);
   gzip_in_file fin_example(example_submission, gzip);
    
    char linebuf[24000];
    char saveptr[1024];
    int added = 0;
    int last_from = -1, last_to = -1, last_rating = 0, last_time = 0;
    int ignore_last_to = -1, ignore_last_from = -1;
    int positive_examples = 0, negative_examples = 0;
    int maxdayofyear = 0; 
    int mindayofyear = 10000000;
    graph_type2 out_graph;
    out_graph.resize(2*nodes);
    graph_type2 out_graph_validation;

    std::vector<singlerating> multiple_ratings; 
 
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
      if (from <= 0){
         logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << " document ID is zero or less: " << from << std::endl;
         return;
      }
      pch = strtok_r(NULL," \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << "[" << linebuf << "]" << std::endl;
        return;
       }
      int to = atoi(pch);
      if (to <= 0){
         logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << " document ID is zero or less: " << from << std::endl;
         return;
      }
      pch = strtok_r(NULL," \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << "[" << linebuf << "]" << std::endl;
        return;
       }
      int rating = atoi(pch);
      if (rating != -1 && rating != 1){
         logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << " invalid rating, not -1 or 1 " << from << std::endl;
         return;
      }
      pch = strtok_r(NULL," \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << "[" << linebuf << "]" << std::endl;
        return;
       }
      int time = atoi(pch);
      //131834878,1321027199
      if (time < 131834878 || time > 1321027199){
         logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << " invalid time " << from << std::endl;
         return;
      }

      total_lines++;
      singlerating thisrating(from, to, rating, time);
      if (debug && (total_lines % 50000 == 0))
        logstream(LOG_INFO) << "Parsed line: " << total_lines << " selected lines: " << added << std::endl;

      //duplicate entry in different time -  store to choose at random one of rating later
      if (last_from == from && last_to == to && last_rating == rating){
         multiple_ratings.push_back(thisrating);
      }
      //item with opposite ratings, ignore
      else if (last_from == from && last_to == to && last_rating != rating){
         ignore_last_from = from;
         ignore_last_to = to;
         multiple_ratings.clear();
      }
      //a different rating encountered
      else if ((last_from > -1) && (last_from != from || last_to != to)){
         //if we marked this entry to ignore it, don't add it.
         if (last_from == ignore_last_from && last_to == ignore_last_to){
            multiple_ratings.clear();
         }
         else { 
            int dayofyear = get_day(last_time);
            if (dayofyear < mindayofyear)
                mindayofyear = dayofyear;
            if (dayofyear > maxdayofyear)
                maxdayofyear = dayofyear;
            assert(multiple_ratings.size() > 0);
            int randint = rand() % multiple_ratings.size();
            singlerating chosen = multiple_ratings[randint];
            edge_data2 edge(chosen.rating, get_day(chosen.time)-283);
            //if (dayofyear >=  split_day_of_year)
            if (chosen.time > split_training_time)
              out_graph_validation.add_edge(chosen.user - 1, chosen.item+nodes-1, edge);
            else out_graph.add_edge(chosen.user - 1, chosen.item+nodes-1, edge);
            added++;
            if (chosen.rating == -1)
	            negative_examples++;
            else positive_examples++;
            multiple_ratings.clear();
            //cout<<"day of year: " << get_day(last_time) << endl;
         }
          multiple_ratings.push_back(thisrating); 
      }
      else multiple_ratings.push_back(thisrating);
   
      last_from = from;
      last_to = to;
      last_rating = rating;
      last_time = time;
    } 

   logstream(LOG_INFO) <<"Finished parsing total of " << total_lines << " lines in file " << vdata.filename << endl;
   logstream(LOG_INFO) <<"Min day of year " << mindayofyear << " maxday " << maxdayofyear << endl;
    
   int dayofyear = get_day(last_time);
   if (dayofyear < mindayofyear)
                mindayofyear = dayofyear;
   if (dayofyear > maxdayofyear)
                maxdayofyear = dayofyear;
   edge_data2 last_edge(last_rating, dayofyear-283/*last_time*/);
   //if (dayofyear >= split_day_of_year)
   if (last_time > split_training_time)
      out_graph_validation.add_edge(last_from - 1, last_to+nodes-1, last_edge);
   else out_graph.add_edge(last_from - 1, last_to+nodes-1, last_edge);
   if (last_rating == -1)
	   negative_examples++;
   else positive_examples++;

    /* handle addition of synthetit time example taken from the test data */
    int examples_added = 0;
    bool first = true; 
    while(true){
      fin_example.get_sp().getline(linebuf, 24000);
      if (fin_example.get_sp().eof())
        break;
      if (first){ //skip the line: id,clicks
        first = false;
        continue;
      }

      char *pch = strtok_r(linebuf," ,\r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing example file: " << endl; 
        return;
       }
      int from = atoi(pch);
      if (from <= 0){
         logstream(LOG_ERROR) << "Error when parsing file: " << example_submission << " line was: " << linebuf <<endl;
         return;
      }
      pch = strtok_r(NULL," ,\r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << "[" << linebuf << "]" << std::endl;
        return;
       }
      int to = atoi(pch), to2 = 0, to3 = 0;
      if (to <= 0){
         logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << total_lines << " document ID is zero or less: " << from << std::endl;
         return;
      }
      pch = strtok_r(NULL," \r\n\t",(char**)&saveptr);
      if (pch){
         to2 = atoi(pch);
         pch = strtok_r(NULL," \r\n\t",(char**)&saveptr);
         if (pch)
            to3 = atoi(pch);
      }
      edge_data2 edge(1,0/* 131834878*/);
      add_single_edge(from, to, edge, out_graph, examples_added);
      if (level >=2 && to2 > 0 ){
        add_single_edge(from, to2, edge, out_graph, examples_added);
       if (level>= 3 && to3> 0){
         add_single_edge(from, to3, edge, out_graph, examples_added);
       }
      } 

   }
   logstream(LOG_ERROR)<<"Added a total of " << examples_added << " from test data" << endl;
   positive_examples+= examples_added;
   out_graph.finalize();
   out_graph_validation.finalize();
    
    bipartite_graph_descriptor training_info, validation_info;
    training_info.rows = training_info.cols = nodes;
    validation_info.rows = validation_info.cols = nodes;
    save_matrixmarket_graph(vdata.filename +".out", training_info, out_graph, MATRIX_MARKET_4);
    logstream(LOG_INFO)<<"Finished exporting a total of " << out_graph.num_edges() << " ratings to training graph " << std::endl << " Positive ratings: " << positive_examples << " Negative ratings: " << negative_examples << std::endl;

  }


};



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Parsers Library");
  clopts.attach_option("data", &datafile, datafile,
                       "training data input file");
  clopts.add_positional("data");
  clopts.attach_option("example_submission", &example_submission, example_submission, "example 2 column submission file location");
  clopts.add_positional("example_submission");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("gzip", &gzip, gzip, "Gzipped input file?");
  clopts.attach_option("split_day_of_year", &split_day_of_year, split_day_of_year, "split training set to validation set, for days >= split_day_of_year");
  clopts.attach_option("level", &level, level, "take XX examples for each user from test data and add to training. 0 means no examples.");
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



