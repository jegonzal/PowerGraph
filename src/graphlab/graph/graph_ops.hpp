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


/**
 * \file graph_io.hpp
 *
 * This file supports basic graph io operations to simplify reading
 * and writing adjacency structures from files.
 *
 */

#ifndef GRAPHLAB_GRAPH_UTILITY_HPP
#define GRAPHLAB_GRAPH_UTILITY_HPP

#include <iostream>
#include <fstream>
#include <string>



namespace graphlab {
  namespace graph_utility {


    /**
     * builds a topological_sort of the graph returning it in topsort. 
     * 
     * \param[out] topsort Resultant topological sort of the graph vertices.
     *
     * function will return false if graph is not acyclic.
     */
    template<typename Graph>
    bool topological_sort(const Graph& graph, 
                          std::vector<vertex_id_type>& topsort) const {
      typedef typename Graph::edge_type edge_type;
      typedef typename Graph::edge_list_type edge_list_type;
      topsort.clear();
      topsort.reserve(num_vertices());
      std::vector<size_t> indeg;
      indeg.resize(graph.num_vertices());
      std::queue<vertex_id_type> q;
      for (size_t i = 0;i < graph.num_vertices(); ++i) {
        indeg[i] = graph.get_in_edges(i).size();
        if (indeg[i] == 0) {
          q.push(i);
        }
      }
    
      while (!q.empty()) {
        vertex_id_type v = q.front();
        q.pop();
        topsort.push_back(v);
        foreach(edge_type edge, graph.get_out_edges(v)) {
          vertex_id_type destv = edge.target();
          --indeg[destv];
          if (indeg[destv] == 0) {
            q.push(destv);
          }
        }
      }
      if (q.empty() && topsort.size() != graph.num_vertices()) 
        return false;
      return true;
    } // end of topological sort


    template<typename Graph>
    size_t unique_neighbors(const Graph& graph) {
      
    } // end of unique_neighbors
    

    template<typename Graph>
    bool load_snap(const std::string& filename,
                   Graph& graph) {
      std::ifstream fin(filename.c_str());
      if(!fin.good()) return false;
      // Loop through file reading each line
      size_t self_edges = 0;
      while(fin.good() && !fin.eof()) {
        if(fin.peek() == '#') {
          std::string str;
          std::getline(fin, str);
          std::cout << str << std::endl;
          continue;
        }
        size_t source = 0;
        size_t target = 0;
        fin >> source;
        if(!fin.good()) break;
        fin >> target; assert(fin.good());
        // Ensure that the number of vertices is correct
        if(source >= graph.num_vertices() || target >= graph.num_vertices())
          graph.resize(std::max(source, target) + 1);
        if(source != target) graph.add_edge(source, target);
        else if(self_edges++ == 0) 
          logstream(LOG_WARNING) 
            << "Self edge encountered but not supported!" << std::endl
            << "\t Further warnings will be surpressed." << std::endl;
      } // end of while loop       
      fin.close();
      logstream(LOG_INFO) 
        << "Finished loading graph with: " << std::endl
        << "\t Vertices: " << graph.num_vertices() << std::endl
        << "\t Edges:  " << graph.num_edges() << std::endl;
      if(self_edges > 0) 
        logstream(LOG_INFO) << "\t Dropped self edges: " << self_edges 
                            << std::endl;
      return true;
    } // end of load SNAP



    template<typename Graph>
    bool load_edge_list(const std::string& filename,
                        graph_type& graph) {
      typedef typename Graph::vertex_id_type vertex_id_type;
      std::ifstream fin(filename.c_str());
      if(!fin.good()) return false;
      size_t self_edges = 0;
      // Loop through file reading each line
      while(fin.good() && !fin.eof()) {
        vertex_id_type source = 0, target = 0;
        fin >> source >> target;
        if(!fin.good()) break;
        // Ensure that the number of vertices is correct
        if(source >= graph.num_vertices() ||
           target >= graph.num_vertices())
          graph.resize(std::max(source, target) + 1);
        if(source != target) graph.add_edge(source, target);
        else if(self_edges++ == 0) 
          logstream(LOG_WARNING) 
            << "Self edge encountered but not supported!" << std::endl
            << "\t Further warnings will be surpressed." << std::endl
      }            
      fin.close();
      logstream(LOG_INFO) 
        << "Finished loading graph with: " << std::endl
        << "\t Vertices: " << graph.num_vertices() << std::endl
        << "\t Edges: " << graph.num_edges() << std::endl;        
      if(self_edges > 0) 
        logstream(LOG_INFO) << "\t Dropped self edges: " << self_edges 
                            << std::endl;  
      return true;
    } // end of load edge list

    
    inline void skip_newline(std::ifstream& fin) {
      char next_char = ' ';
      fin.get(next_char);
      ASSERT_EQ(next_char, '\n');  
    }

    template<typename Graph>
    bool load_metis(const std::string& filename,
                    Graph& graph) { 
      typedef typename Graph::vertex_id_type vertex_id_type;
      std::ifstream fin(filename.c_str());
      if(!fin.good()) return false;
      size_t nverts = 0, nedges = 0;
      fin >> nverts >> nedges;
      logstream(LOG_INFO) 
        << "Loading graph with: " << std::endl
        << "\t Vertices: " << nverts << std::endl
        << "\t Edges: " << nedges << std::endl;          
      skip_newline(fin);
      graph.resize(nverts);
      for(vertex_id_type source = 0; source < nverts; ++source) {
        while(fin.peek() != '\n') {
          ASSERT_TRUE(fin.good());
          vertex_id_type target = 0;
          fin >> target; 
          ASSERT_GT(target, 0);
          // decrement the value since starting value is 1 not zero
          target--; 
          ASSERT_LT(target, graph.num_vertices());     
          if(source != target) graph.add_edge(source, target);
                  else if(self_edges++ == 0) 
          logstream(LOG_WARNING) 
            << "Self edge encountered but not supported!" << std::endl
            << "\t Further warnings will be surpressed." << std::endl;
          
        }
        skip_newline(fin);
      }
      fin.close();
      logstream(LOG_INFO) 
        << "Finished loading graph with: " << std::endl
        << "\t Vertices: " << graph.num_vertices() << std::endl
        << "\t Edges: " << graph.num_edges() << std::endl;      
      if(self_edges > 0) 
        logstream(LOG_INFO) << "\t Dropped self edges: " << self_edges 
                            << std::endl;
      return true;
    } // end of load metis



    template<typename Graph>
    bool load(const std::string& fname,
              const std::string& format,
              Graph& graph) {
      if (format == "metis") return load_metis(fname, graph);
      else if (format == "snap") return load_snap(fname, graph);
      else if (format == "tsv") return load_edge_list(fname, graph);
      else {
        logstream(LOG_WARNING)
          << "Invalid format \"" << format << "\".  "
          << "Unable to load file \"" << fname << "\"!" << std::endl;     
      }
      return false;
    }



    template<typename Graph>
    bool load(const std::string& fname,
              Graph& graph) {
      const size_t pos = fname.rfind('.');
      if(pos == std::string::npos || pos + 1 >= fname.size()) {
        logstream(LOG_WARNING) 
          << "Filename \"" << fname 
          << "\" does not have a suffix." << std::endl
          << "Unable to infer file format!" << std::endl;
        return false;
      }      
      const std::string format(fname.substr(pos+1, std::string::npos));
      logstream(LOG_INFO) << "File format: " << format << std::endl;
      return load(fname, format, graph);
    } // end of load













    template<typename Graph>
    bool load_snap(const std::string& filename,
                   Graph& graph) {
      std::ifstream fin(filename.c_str());
      if(!fin.good()) return false;
      // Loop through file reading each line
      size_t self_edges = 0;
      while(fin.good() && !fin.eof()) {
        if(fin.peek() == '#') {
          std::string str;
          std::getline(fin, str);
          std::cout << str << std::endl;
          continue;
        }
        size_t source = 0;
        size_t target = 0;
        fin >> source;
        if(!fin.good()) break;
        fin >> target; assert(fin.good());
        // Ensure that the number of vertices is correct
        if(source >= graph.num_vertices() || target >= graph.num_vertices())
          graph.resize(std::max(source, target) + 1);
        if(source != target) graph.add_edge(source, target);
        else if(self_edges++ == 0) 
          logstream(LOG_WARNING) 
            << "Self edge encountered but not supported!" << std::endl
            << "\t Further warnings will be surpressed." << std::endl;
      } // end of while loop       
      fin.close();
      logstream(LOG_INFO) 
        << "Finished loading graph with: " << std::endl
        << "\t Vertices: " << graph.num_vertices() << std::endl
        << "\t Edges:  " << graph.num_edges() << std::endl;
      if(self_edges > 0) 
        logstream(LOG_INFO) << "\t Dropped self edges: " << self_edges 
                            << std::endl;
      return true;
    } // end of load SNAP



    template<typename Graph>
    bool load_edge_list(const std::string& filename,
                        graph_type& graph) {
      typedef typename Graph::vertex_id_type vertex_id_type;
      std::ifstream fin(filename.c_str());
      if(!fin.good()) return false;
      size_t self_edges = 0;
      // Loop through file reading each line
      while(fin.good() && !fin.eof()) {
        vertex_id_type source = 0, target = 0;
        fin >> source >> target;
        if(!fin.good()) break;
        // Ensure that the number of vertices is correct
        if(source >= graph.num_vertices() ||
           target >= graph.num_vertices())
          graph.resize(std::max(source, target) + 1);
        if(source != target) graph.add_edge(source, target);
        else if(self_edges++ == 0) 
          logstream(LOG_WARNING) 
            << "Self edge encountered but not supported!" << std::endl
            << "\t Further warnings will be surpressed." << std::endl
      }            
      fin.close();
      logstream(LOG_INFO) 
        << "Finished loading graph with: " << std::endl
        << "\t Vertices: " << graph.num_vertices() << std::endl
        << "\t Edges: " << graph.num_edges() << std::endl;        
      if(self_edges > 0) 
        logstream(LOG_INFO) << "\t Dropped self edges: " << self_edges 
                            << std::endl;  
      return true;
    } // end of load edge list

    

    template<typename Graph>
    bool save_metis(const std::string& filename,
                    const Graph& graph) { 
      typedef typename Graph::vertex_id_type vertex_id_type;
      std::ofstream fout(filename.c_str());
      if(!fout.good()) return false;
      size_t nverts = 0, nedges = 0;
      fout << graph.num_vertices() << graph.num_edges() << '\n';
      for(vertex_id_type source = 0; source < graph.num_vertices(); 
          ++source) {


      }
      for(vertex_id_type source = 0; source < nverts; ++source) {
        while(fin.peek() != '\n') {
          ASSERT_TRUE(fin.good());
          vertex_id_type target = 0;
          fin >> target; 
          ASSERT_GT(target, 0);
          // decrement the value since starting value is 1 not zero
          target--; 
          ASSERT_LT(target, graph.num_vertices());     
          if(source != target) graph.add_edge(source, target);
                  else if(self_edges++ == 0) 
          logstream(LOG_WARNING) 
            << "Self edge encountered but not supported!" << std::endl
            << "\t Further warnings will be surpressed." << std::endl;
          
        }
        skip_newline(fin);
      }
      fin.close();
      logstream(LOG_INFO) 
        << "Finished loading graph with: " << std::endl
        << "\t Vertices: " << graph.num_vertices() << std::endl
        << "\t Edges: " << graph.num_edges() << std::endl;      
      if(self_edges > 0) 
        logstream(LOG_INFO) << "\t Dropped self edges: " << self_edges 
                            << std::endl;
      return true;
    } // end of save metis



    template<typename Graph>
    bool save(const std::string& fname,
              const std::string& format,
              const Graph& graph) {
      if (format == "metis") return save_metis(fname, graph);
      else if (format == "snap") return save_snap(fname, graph);
      else if (format == "tsv") return save_edge_list(fname, graph);
      else {
        logstream(LOG_WARNING)
          << "Invalid format \"" << format << "\".  "
          << "Unable to save file \"" << fname << "\"!" << std::endl;     
      }
      return false;
    } // end of save


    template<typename Graph>
    bool save(const std::string& fname,
              const Graph& graph) {
      const size_t pos = fname.rfind('.');
      if(pos == std::string::npos || pos + 1 >= fname.size()) {
        logstream(LOG_WARNING) 
          << "Filename \"" << fname 
          << "\" does not have a suffix." << std::endl
          << "Unable to infer file format!" << std::endl;
        return false;
      }      
      const std::string format(fname.substr(pos+1, std::string::npos));
      logstream(LOG_INFO) << "File format: " << format << std::endl;
      return save(fname, format, graph);
    } // end of save



  }; // end of graph io
}; // end of namespace graphlab
#endif





