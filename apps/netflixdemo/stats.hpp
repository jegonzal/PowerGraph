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
 * Class for collecting graph statistics
 */


#ifndef TK_STATS
#define TK_STATS

  
struct stats_info{
  size_t validation_edges;
  size_t training_edges;
  size_t max_user;
  size_t max_item;
  
  stats_info(){
    validation_edges = training_edges = max_user = max_item = 0;
  }

  stats_info & operator+=(const stats_info & other){
     validation_edges += other.validation_edges;
     training_edges += other.training_edges;
     max_user = std::max(max_user, other.max_user);
     max_item = std::max(max_item, other.max_item); 
     return *this;
  }

  /** \brief Save the values to a binary archive */
  void save(graphlab::oarchive& arc) const { arc << validation_edges << training_edges << max_user << max_item; }

  /** \brief Read the values from a binary archive */
  void load(graphlab::iarchive& arc) { arc >> validation_edges >> training_edges >> max_user >> max_item; }  

};

  
stats_info info;


#endif //TK_STATS
