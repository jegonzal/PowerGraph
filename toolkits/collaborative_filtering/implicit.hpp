#ifndef _IMPLICIT_HPP__
#define _IMPLICIT_HPP__
/**
 * @file
 * @author  Danny Bickson
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 * header file for handling the addition of implicit edges
 */

#include "eigen_wrapper.hpp"
#include "stats.hpp"

enum{
  IMPLICIT_RATING_DISABLED = 0,
  IMPLICIT_RATING_RANDOM = 1
};

double implicitratingweight;
double implicitratingvalue;
double implicitratingpercentage;
int    implicitratingtype;

template<typename als_edge_type>
uint add_implicit_edges4(int type, graph_type & graph, graphlab::distributed_control & dc){

  switch(type){
    case IMPLICIT_RATING_DISABLED: return 0;
    case IMPLICIT_RATING_RANDOM: break;
    default: assert(false);
  };

  uint added = 0;
  size_t M = info.max_user;
  size_t N = info.max_item;
  uint toadd  = implicitratingpercentage*N*M;
  dc.cout()<<"Going to add: " << toadd << " implicit edges. users: " << M << " items: " << N << std::endl;
  assert(toadd >= 1);
  for (uint j=0; j< toadd; j++){
    ivec item = ::randi(1,0,N-1);
    ivec user = ::randi(1,0,M-1);
    graph.add_edge(user[0], -(graphlab::vertex_id_type(item[0] + SAFE_NEG_OFFSET)), als_edge_type(implicitratingvalue, edge_data::TRAIN, implicitratingweight));
    added++;
  } 
  dc.cout()<<"Finished adding " << toadd << " implicit edges. " << std::endl;
  return added;
};

template<typename als_edge_type>
uint add_implicit_edges(int type, graph_type & graph, graphlab::distributed_control & dc){

  switch(type){
    case IMPLICIT_RATING_DISABLED: return 0;
    case IMPLICIT_RATING_RANDOM: break;
    default: assert(false);
  };

  uint added = 0;
  size_t M = info.max_user;
  size_t N = info.max_item;
  uint toadd  = implicitratingpercentage*N*M;
  dc.cout()<<"Going to add: " << toadd << " implicit edges. users: " << M << " items: " << N <<std::endl;
  assert(toadd >= 1);
  for (uint j=0; j< toadd; j++){
    ivec item = ::randi(1,0,N-1);
    ivec user = ::randi(1,0,M-1);
    graph.add_edge(user[0], -(graphlab::vertex_id_type(item[0] + SAFE_NEG_OFFSET)), als_edge_type(implicitratingvalue));
    added++;
  } 
  dc.cout()<<"Finished adding " << toadd << " implicit edges. " << std::endl;
  return added;
};

void parse_implicit_command_line(graphlab::command_line_options & clopts){
   clopts.attach_option("implicitratingweight", implicitratingweight,"implicit rating weight");
   clopts.attach_option("implicitratingvalue", implicitratingvalue, "implicit rating value");
   clopts.attach_option("implicitratingtype", implicitratingtype, "implicit rating type (-=disabled, 1=random)");
   if (implicitratingtype != IMPLICIT_RATING_RANDOM && implicitratingtype != IMPLICIT_RATING_DISABLED)
     logstream(LOG_FATAL)<<"Implicit rating type should be either 0 (IMPLICIT_RATING_DISABLED) or 1 (IMPLICIT_RATING_RANDOM)" << std::endl;
   clopts.attach_option("implicitratingpercentage", implicitratingpercentage, "implicit rating percentage (1e-8,0.8)");
   if (implicitratingpercentage < 1e-8 && implicitratingpercentage > 0.8)
     logstream(LOG_FATAL)<<"Implicit rating percentage should be (1e-8, 0.8)" << std::endl;
  clopts.attach_option("users", info.max_user, "max user id (for implicit ratings)");
  clopts.attach_option("items", info.max_item, "max item id (for implicit ratings)");

}
#endif //_IMPLICIT_HPP__
