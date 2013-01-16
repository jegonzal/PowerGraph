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

#ifndef GRAPHLAB_DISTRIBUTED_SHARDING_CONSTRAINT_HPP
#define GRAPHLAB_DISTRIBUTED_SHARDING_CONSTRAINT_HPP

#include <graphlab/graph/graph_basic_types.hpp>
#include <algorithm>
#include <vector>

namespace graphlab {
  class sharding_constraint {
    size_t nshards;
    std::vector<std::vector<procid_t> > constraint_graph;
   public:
    sharding_constraint(size_t num_shards, std::string method) {
      nshards = num_shards;
      // ignore the method input for now, only construct grid graph. 
      // assuming nshards is perfect square
      make_grid_constraint();
      check();
    }

    bool get_neighbors (procid_t shard, std::vector<procid_t>& neighbors) {
      ASSERT_LT(shard, nshards);
      neighbors.clear();
      std::vector<procid_t>& ls = constraint_graph[shard];
      for (size_t i = 0; i < ls.size(); ++i)
        neighbors.push_back(ls[i]);
      return true;
    }

    bool get_joint_neighbors (procid_t shardi, procid_t shardj, std::vector<procid_t>& neighbors) {
      ASSERT_LT(shardi, nshards);
      ASSERT_LT(shardj, nshards);
      std::vector<procid_t>& ls1 = constraint_graph[shardi];
      std::vector<procid_t>& ls2 = constraint_graph[shardj];
      neighbors.clear();
      size_t i = 0;
      size_t j = 0;
      while (i < ls1.size() && j < ls2.size()) {
        if (ls1[i] == ls2[j]) {
          neighbors.push_back(ls1[i]);
          ++i; ++j;
        } else if (ls1[i] < ls2[j]) {
          ++i;
        } else {
          ++j;
        }
      }
      return true;
    }

   private:
    void make_grid_constraint() {
      size_t ncols, nrows;
      ncols = nrows = (size_t)sqrt(nshards);

      for (size_t i = 0; i < nshards; i++) {
        std::vector<procid_t> adjlist;
        // add self
        adjlist.push_back(i);

        // add the row of i
        size_t rowbegin = (i/ncols) * ncols;
        for (size_t j = rowbegin; j < rowbegin + ncols; ++j)
          if (i != j) adjlist.push_back(j); 

        // add the col of i
        for (size_t j = i % ncols; j < nshards; j+=ncols)
          if (i != j) adjlist.push_back(j); 

        std::sort(adjlist.begin(), adjlist.end());
        constraint_graph.push_back(adjlist);
      }
    }

    void check() {
      // debug 
      for (size_t i = 0; i < constraint_graph.size(); ++i) {
        std::vector<procid_t> adjlist = constraint_graph[i];
        std::cout << i << ": [";
        for (size_t j = 0; j < adjlist.size(); j++)
          std::cout << adjlist[j] << " ";
        std::cout << "]" << std::endl;
      }

      for (size_t i = 0; i < nshards; ++i) {
        for (size_t j = i+1; j < nshards; ++j) {
          std::vector<procid_t> ls;
          get_joint_neighbors(i, j, ls);
          ASSERT_GT(ls.size(), 0);
        }
      }
    }
    }; // end of sharding_constraint
}; // end of namespace graphlab
#endif
