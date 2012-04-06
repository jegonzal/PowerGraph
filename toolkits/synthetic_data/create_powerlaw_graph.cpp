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


#include <iostream>
#include <algorithm>
#include <set>
#include <fstream>


#include <boost/unordered_set.hpp>


#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/plod_generator.hpp>
#include <boost/random/linear_congruential.hpp>

#include <graphlab/util/fast_multinomial.hpp>

#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>


std::string algorithm = "constant";
std::string fname = "graph.tsv";
double alpha = 2.1;
double beta = 10000;
size_t nverts = 10;
size_t fanout = 2;


void pdf2cdf(std::vector<double>& pdf) {
  double Z = 0;
  for(size_t i = 0; i < pdf.size(); ++i) Z += pdf[i];
  for(size_t i = 0; i < pdf.size(); ++i) 
    pdf[i] = pdf[i]/Z + ((i>0)? pdf[i-1] : 0);
} // end of pdf2cdf

size_t sample(const std::vector<double>& cdf) {
  return std::upper_bound(cdf.begin(), cdf.end(), 
                          graphlab::random::rand01()) - cdf.begin();  
} // end of sample

void constant_fanin(bool count_edges_only) {
  std::vector<double> prob(nverts-1, 0);
  std::cout << "constructing pdf" << std::endl;
  for(size_t i = 0; i < prob.size(); ++i) 
    prob[i] = std::pow(double(i+1), -alpha);
  std::cout << "constructing cdf" << std::endl;
  pdf2cdf(prob);
  std::cout << "sampling out_degrees" << std::endl;
  std::vector<size_t> out_degree(nverts, 0);
  size_t edge_count = 0;
  for(size_t i = 0; i < nverts; ++i) 
    edge_count += (out_degree[i] = sample(prob) + 1);
  std::cout << "Edges: " << edge_count << std::endl;
  if(count_edges_only) return;

  std::cout << "Saving out_degree distribution" << std::endl;
  const std::string out_degree_fname = fname + ".out_degree";
  std::ofstream out_degree_fout(out_degree_fname.c_str());
  for(size_t i = 0; i < out_degree.size(); ++i) 
    out_degree_fout << out_degree[i] << '\n';
  out_degree_fout.close();
  std::cout << "Sampling graph" << std::endl;
  std::ofstream fout(fname.c_str());  
  boost::unordered_set<size_t> targets;
  std::vector<size_t> in_degree(nverts);
  for(size_t source = 0; source < nverts; ++source) {
    targets.clear();
    // Determine if we are computing out_degree targets are nverts -
    // out_degree not targets
    const size_t tmp_out_degree = out_degree[source] < nverts/2? 
      out_degree[source] : nverts - out_degree[source];
    ASSERT_LE(tmp_out_degree, nverts);
    // uniformly sample targets (or anit-targets)
    while(targets.size() < tmp_out_degree) {
      const size_t target = graphlab::random::fast_uniform<size_t>(0, nverts-1);
      if(source != target) targets.insert(target);
    }
    // if tmp_out_degree == out_degree then we were sampling targets
    if(tmp_out_degree == out_degree[source]) {
      foreach(size_t target, targets) {
        fout << source << '\t' << target << '\n';
        in_degree[target]++;
      }
    } else {
      // Otherwise we sampled anti targets so everything that is not
      // in the targets set is a target.
      size_t count = 0;
      for(size_t target = 0; target < nverts; ++target) {
        if(targets.count(target) == 0) {
          fout << source << '\t' << target << '\n';
          count++;
          in_degree[target]++;
        }
      }
      ASSERT_EQ(count, out_degree[source]);
    }  
  }
  fout.close();

  std::cout << "Saving in_degree distribution" << std::endl;
  const std::string in_degree_fname = fname + ".in_degree";
  std::ofstream in_degree_fout(in_degree_fname.c_str());
  for(size_t i = 0; i < in_degree.size(); ++i) 
    in_degree_fout << in_degree[i] << '\n';
  in_degree_fout.close();

} // end of constant fanout powerlaw

void boost_powerlaw() {
  typedef boost::adjacency_list<> Graph;
  typedef boost::plod_iterator<boost::minstd_rand, Graph> SFGen;
  
  boost::minstd_rand gen;
  // Create graph with 100 nodes 
  Graph graph(SFGen(gen, nverts, alpha, fanout, false), 
              SFGen(), nverts);

  boost::graph_traits<Graph>::edge_iterator edge, edge_end;
  std::ofstream fout(fname.c_str());
  for( boost::tie(edge, edge_end) = boost::edges(graph); edge != edge_end; ++edge)
    fout << boost::source(*edge, graph) << '\t'
         << boost::target(*edge, graph) << '\n';
  fout.close();
} // end of boost powerlaw



void preferential_attachment() {
  graphlab::fast_multinomial multi(nverts, 1);
  for(size_t i = 0; i < nverts; ++i) multi.set(i, beta);
  boost::unordered_set<size_t> targets;
  std::ofstream fout(fname.c_str());
  for(size_t source = 0; source < nverts; ++source) {
    targets.clear();
    while(targets.size() < fanout) {
      size_t target(-1);
      multi.sample(target, 0);
      if(source != target) {
        multi.add(target, 1.0);
        targets.insert(target);
      }
    }
    multi.add(source, 1.0);
    foreach(size_t target, targets)
      fout << source << '\t' << target << '\n';
  }
  fout.close();
} // end of boost powerlaw



int main(int argc, char** argv) {
  
  graphlab::command_line_options 
    clopts("Generate synthetic graph data.", true);
  bool count_edges_only = false;
  clopts.attach_option("alg", &algorithm, algorithm, 
                       "The algorithm to use."); 
  clopts.attach_option("count_edges", &count_edges_only, 
                       count_edges_only, "just count edges and return");
  clopts.attach_option("graph", &fname, fname,
                       "The name of the graph file."); 
  clopts.attach_option("alpha", &alpha, alpha, 
                       "Power law constant:  beta * deg^{-alpha}.");
  clopts.attach_option("beta", &beta, beta, 
                       "Power law constant:  beta * deg^{-beta}.");

  clopts.attach_option("nverts", &nverts, nverts,
                       "Number of vertices.");
  clopts.attach_option("fanout", &fanout, fanout,
                       "The fanout of each page");

  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  if(algorithm == "constant") constant_fanin(count_edges_only);
  else if(algorithm == "boost") boost_powerlaw();
  else if(algorithm == "preferential") preferential_attachment();
  else {
    std::cout << "Invalid algorithm type \"" << algorithm  << "\" valid types are: " 
              << std::endl
              << "\t constant: draw out neigbors according too a powerlaw \n"
              << "\t           degree distribution.\n"
              << "\t boost: use the boost power law graph algorithm.\n"
              << "\t preferential: use the preferential attachment algorithm.\n"
              << std::endl;
  }
 
  return EXIT_SUCCESS;
} // end of main
