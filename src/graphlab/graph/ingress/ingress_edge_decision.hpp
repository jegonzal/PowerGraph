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
#ifndef GRAPHLAB_DISTRIBUTED_INGRESS_EDGE_DECISION_HPP
#define GRAPHLAB_DISTRIBUTED_INGRESS_EDGE_DECISION_HPP

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/graph/distributed_graph.hpp>

namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;
 
 template<typename VertexData, typename EdgeData>
 class ingress_edge_decision {

    public:
      typedef graphlab::vertex_id_type vertex_id_type;
      typedef distributed_graph<VertexData, EdgeData> graph_type;
      typedef fixed_dense_bitset<RPC_MAX_N_PROCS> bin_counts_type; 


    public:
      /** \brief A decision object for computing the edge assingment. */
      ingress_edge_decision(distributed_control& dc) {
      }

      /** Random assign (source, target) to a machine p in {0, ... numprocs-1} */
      procid_t edge_to_proc_random (const vertex_id_type source, 
          const vertex_id_type target,
          size_t numprocs) {
        typedef std::pair<vertex_id_type, vertex_id_type> edge_pair_type;
        boost::hash< edge_pair_type >  hash_function;
        const edge_pair_type edge_pair(std::min(source, target), 
            std::max(source, target));
        return hash_function(edge_pair) % (numprocs);
      };

      /** Random assign (source, target) to a machine p in a list of candidates */
      procid_t edge_to_proc_random (const vertex_id_type source, 
          const vertex_id_type target,
          const std::vector<procid_t> & candidates) {
        typedef std::pair<vertex_id_type, vertex_id_type> edge_pair_type;
        boost::hash< edge_pair_type >  hash_function;
        const edge_pair_type edge_pair(std::min(source, target), 
            std::max(source, target));
        return candidates[hash_function(edge_pair) % (candidates.size())];
      };


      /** Greedy assign (source, target) to a machine using: 
       *  bitset<MAX_MACHINE> src_degree : the degree presence of source over machines
       *  bitset<MAX_MACHINE> dst_degree : the degree presence of target over machines
       *  vector<size_t>      proc_num_edges : the edge counts over machines
       * */
      procid_t edge_to_proc_greedy (const vertex_id_type source, 
          const vertex_id_type target,
          bin_counts_type& src_degree,
          bin_counts_type& dst_degree,
          std::vector<size_t>& proc_num_edges,
          bool usehash = false,
          bool userecent = false
          ) {
        size_t numprocs = proc_num_edges.size();

        // Compute the score of each proc.
        procid_t best_proc = -1; 
        double maxscore = 0.0;
        double epsilon = 1.0; 
        std::vector<double> proc_score(numprocs); 
        size_t minedges = *std::min_element(proc_num_edges.begin(), proc_num_edges.end());
        size_t maxedges = *std::max_element(proc_num_edges.begin(), proc_num_edges.end());

        for (size_t i = 0; i < numprocs; ++i) {
          size_t sd = src_degree.get(i) + (usehash && (source % numprocs == i));
          size_t td = dst_degree.get(i) + (usehash && (target % numprocs == i));
          double bal = (maxedges - proc_num_edges[i])/(epsilon + maxedges - minedges);
          proc_score[i] = bal + ((sd > 0) + (td > 0));
        }
        maxscore = *std::max_element(proc_score.begin(), proc_score.end());

        std::vector<procid_t> top_procs; 
        for (size_t i = 0; i < numprocs; ++i)
          if (std::fabs(proc_score[i] - maxscore) < 1e-5)
            top_procs.push_back(i);

        // Hash the edge to one of the best procs.
        typedef std::pair<vertex_id_type, vertex_id_type> edge_pair_type;
        boost::hash< edge_pair_type >  hash_function;
        const edge_pair_type edge_pair(std::min(source, target), 
            std::max(source, target));
        best_proc = top_procs[hash_function(edge_pair) % top_procs.size()];

        ASSERT_LT(best_proc, numprocs);
        if (userecent) {
          src_degree.clear();
          dst_degree.clear();
        }
        src_degree.set_bit(best_proc);
        dst_degree.set_bit(best_proc);
        ++proc_num_edges[best_proc];
        return best_proc;
      };

      /** Greedy assign (source, target) to a machine using: 
       *  bitset<MAX_MACHINE> src_degree : the degree presence of source over machines
       *  bitset<MAX_MACHINE> dst_degree : the degree presence of target over machines
       *  vector<size_t>      proc_num_edges : the edge counts over machines
       * */
      procid_t edge_to_proc_greedy (const vertex_id_type source, 
          const vertex_id_type target,
          bin_counts_type& src_degree,
          bin_counts_type& dst_degree,
          std::vector<procid_t>& candidates,
          std::vector<size_t>& proc_num_edges,
          bool usehash = false,
          bool userecent = false
          ) {
        size_t numprocs = proc_num_edges.size();

        // Compute the score of each proc.
        procid_t best_proc = -1; 
        double maxscore = 0.0;
        double epsilon = 1.0; 
        std::vector<double> proc_score(candidates.size()); 
        size_t minedges = *std::min_element(proc_num_edges.begin(), proc_num_edges.end());
        size_t maxedges = *std::max_element(proc_num_edges.begin(), proc_num_edges.end());

        for (size_t j = 0; j < candidates.size(); ++j) {
          size_t i = candidates[j];
          size_t sd = src_degree.get(i) + (usehash && (source % numprocs == i));
          size_t td = dst_degree.get(i) + (usehash && (target % numprocs == i));
          double bal = (maxedges - proc_num_edges[i])/(epsilon + maxedges - minedges);
          proc_score[j] = bal + ((sd > 0) + (td > 0));
        }
        maxscore = *std::max_element(proc_score.begin(), proc_score.end());

        std::vector<procid_t> top_procs; 
        for (size_t j = 0; j < candidates.size(); ++j)
          if (std::fabs(proc_score[j] - maxscore) < 1e-5)
            top_procs.push_back(candidates[j]);

        // Hash the edge to one of the best procs.
        typedef std::pair<vertex_id_type, vertex_id_type> edge_pair_type;
        boost::hash< edge_pair_type >  hash_function;
        const edge_pair_type edge_pair(std::min(source, target), 
            std::max(source, target));
        best_proc = top_procs[hash_function(edge_pair) % top_procs.size()];

        ASSERT_LT(best_proc, numprocs);
        if (userecent) {
          src_degree.clear();
          dst_degree.clear();
        }
        src_degree.set_bit(best_proc);
        dst_degree.set_bit(best_proc);
        ++proc_num_edges[best_proc];
        return best_proc;
      };


      ////////////////////// Deprecated ////////////////////////////////////

      /** Greedy assign (source, target) to a machine using:
       *  vector<size_t> src_degree : the degree counts of source over machines
       *  vector<size_t> dst_degree : the degree counts of target over machines
       *  vector<size_t> proc_num_edges : the edge counts over machines
       * */
      procid_t edge_to_proc_greedy (const vertex_id_type source, 
          const vertex_id_type target,
          std::vector<size_t>& src_degree,
          std::vector<size_t>& dst_degree,
          std::vector<size_t>& proc_num_edges,
          bool usehash = false,
          bool userecent = false) {
        size_t numprocs = proc_num_edges.size();
        if (src_degree.size() == 0)
          src_degree.resize(numprocs, 0);
        if (dst_degree.size() == 0)
          dst_degree.resize(numprocs, 0);

        // Compute the score of each proc.
        procid_t best_proc = -1; 
        double maxscore = 0.0;
        double epsilon = 1.0; 
        std::vector<double> proc_score(numprocs); 
        size_t minedges = *std::min_element(proc_num_edges.begin(), proc_num_edges.end());
        size_t maxedges = *std::max_element(proc_num_edges.begin(), proc_num_edges.end());
        for (size_t i = 0; i < numprocs; ++i) {
          size_t sd = src_degree[i] + (usehash && (source % numprocs == i));
          size_t td = dst_degree[i] + (usehash && (target % numprocs == i));
          double bal = (maxedges - proc_num_edges[i])/(epsilon + maxedges - minedges);
          proc_score[i] = bal + ((sd > 0) + (td > 0));
        }
        maxscore = *std::max_element(proc_score.begin(), proc_score.end());

        std::vector<procid_t> top_procs; 
        for (size_t i = 0; i < numprocs; ++i)
          if (std::fabs(proc_score[i] - maxscore) < 1e-5)
            top_procs.push_back(i);

        if (top_procs.size() > 1) {
          if (maxscore >= 2) {
            //PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_BOTH_TIE, 1)
          } else if (maxscore >= 1) {
            //PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_ONE_TIE, 1)
          } else {
            //PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_NONE_TIE, 1); 
          }
        } else {
          if (maxscore >= 2) {
            //PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_BOTH_UNIQUE, 1);
          } else if (maxscore >= 1) {
            //PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_ONE_UNIQUE, 1)
          } else {
            //PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_NONE_UNIQUE, 1); 
          }
        }

        // Hash the edge to one of the best procs.
        typedef std::pair<vertex_id_type, vertex_id_type> edge_pair_type;
        boost::hash< edge_pair_type >  hash_function;
        const edge_pair_type edge_pair(std::min(source, target), 
            std::max(source, target));
        best_proc = top_procs[hash_function(edge_pair) % top_procs.size()];

        ASSERT_LT(best_proc, numprocs);

        // only counts the recent proc.
        if (userecent) {
          src_degree.clear();
          dst_degree.clear();
        }

        ++src_degree[best_proc];
        ++dst_degree[best_proc];
        ++proc_num_edges[best_proc];
        return best_proc;
      };

  };// end of ingress_edge_decision




}

#endif
