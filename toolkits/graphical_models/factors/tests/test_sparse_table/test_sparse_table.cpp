/**  
 *  Software submitted by 
 *  Systems & Technology Research / Vision Systems Inc., 2013
 *
 *  Approved for public release; distribution is unlimited. [DISTAR Case #21428]
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
 * Test various functions of a sparse table
 *
 *  \author Scott Richardson 
 */

#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <iterator>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <factors/sparse_table.hpp>
#include <factors/dense_table.hpp>


const size_t MAX_DIM = 4;
typedef graphlab::dense_table<MAX_DIM>          dense_table_t;
typedef graphlab::sparse_table<MAX_DIM>         sparse_table_t;
typedef graphlab::discrete_assignment<MAX_DIM>  assignment_t;
typedef graphlab::discrete_domain<MAX_DIM>      domain_t;
typedef graphlab::discrete_variable             variable_t;


sparse_table_t create_rand_sparse_table(unsigned v0_id, unsigned v1_id, unsigned v2_id) 
{
  variable_t v0(v0_id, 4);
  variable_t v1(v1_id, 3);
  variable_t v2(v2_id, 2);

  std::vector<variable_t> args;
  args.push_back(v0);
  args.push_back(v1);
  args.push_back(v2);

  domain_t dom(args); 
  sparse_table_t st(dom);

  domain_t::const_iterator asg_it = dom.begin();
  domain_t::const_iterator end = dom.end();
  for( ; asg_it != end; ++asg_it) {
    if(rand() % 100 <= 20) {
      st.set_logP( *asg_it, -1 * (rand() % 100) );
    }
  }

  return st;
}

// [-Inf -Inf -3   -Inf]
// [-15  -4   -Inf -23 ]
// [-Inf -Inf -Inf -Inf]
//
// [-Inf -Inf -Inf -Inf]
// [-20  -12  -19  -78 ]
// [-Inf -Inf -Inf -32 ]
void create_data_vector(std::vector<std::pair<size_t, double> > &data, 
    const domain_t &dom) 
{
  {
  size_t sa[] = {2,1,1}; 
  assignment_t da(dom, std::vector<size_t>(sa, sa+3));
  data.push_back(std::make_pair(da.linear_index(), -19));
  } {  
  size_t sa[] = {3,2,1}; 
  assignment_t da(dom, std::vector<size_t>(sa, sa+3));
  data.push_back(std::make_pair(da.linear_index(), -32));
  } {
  size_t sa[] = {2,0,0}; 
  assignment_t da(dom, std::vector<size_t>(sa, sa+3));
  data.push_back(std::make_pair(da.linear_index(), -3));
  } {
  size_t sa[] = {3,1,0}; 
  assignment_t da(dom, std::vector<size_t>(sa, sa+3));
  data.push_back(std::make_pair(da.linear_index(), -23));
  } {
  size_t sa[] = {3,1,1}; 
  assignment_t da(dom, std::vector<size_t>(sa, sa+3));
  data.push_back(std::make_pair(da.linear_index(), -78));
  } {
  size_t sa[] = {1,1,1}; 
  assignment_t da(dom, std::vector<size_t>(sa, sa+3));
  data.push_back(std::make_pair(da.linear_index(), -12));
  } {
  size_t sa[] = {1,1,0}; 
  assignment_t da(dom, std::vector<size_t>(sa, sa+3));
  data.push_back(std::make_pair(da.linear_index(), -4));
  } {
  size_t sa[] = {0,1,0}; 
  assignment_t da(dom, std::vector<size_t>(sa, sa+3));
  data.push_back(std::make_pair(da.linear_index(), -15));
  } {  
  size_t sa[] = {0,1,1}; 
  assignment_t da(dom, std::vector<size_t>(sa, sa+3));
  data.push_back(std::make_pair(da.linear_index(), -20));
  }
}

void testDataReorder(unsigned v0_id, unsigned v1_id, unsigned v2_id) {
  assert(v0_id != v1_id && v1_id != v2_id && v0_id != v2_id);

  // setup the data as if its from belief prop
  std::vector<std::pair<size_t, double> > data;
  variable_t v0(0, 4);
  variable_t v1(1, 3);
  variable_t v2(2, 2);

  std::vector<variable_t> args;
  args.push_back(v0);
  args.push_back(v1);
  args.push_back(v2);
  domain_t dom(args); 

  create_data_vector(data, dom);

  // create a table with variables that are not resorted
  sparse_table_t st_gm(args, data);
  //std::cout << "st_gm " << st_gm << std::endl;


  variable_t r0(v0_id, 4);
  variable_t r1(v1_id, 3);
  variable_t r2(v2_id, 2);
  std::vector<variable_t> reordered_args;
  reordered_args.push_back(r0);
  reordered_args.push_back(r1);
  reordered_args.push_back(r2);
  sparse_table_t st(reordered_args, data);
  //std::cout << "st " << st << std::endl;

  assignment_t asg_gm(st_gm.domain());
  domain_t::const_iterator asg_it = dom.begin();
  domain_t::const_iterator end = dom.begin();
  for( ; asg_it != end; ++asg_it) {
    asg_gm.set_asg(v0, asg_it->asg(r0));
    asg_gm.set_asg(v1, asg_it->asg(r1));
    asg_gm.set_asg(v2, asg_it->asg(r2));
    
    //std::cout << "st_gm{" << asg_gm << "}=" << st_gm.logP(asg_gm) << "; "
    //          << "st{" << asg_it << "}=" << st.logP(asg_it) << std::endl;
    ASSERT_TRUE(st_gm.logP(asg_gm) == st.logP(*asg_it));
  }
}

sparse_table_t create_sparse_table(unsigned v0_id, unsigned v1_id, unsigned v2_id) {
  // create a table with dimensions ordered like the original table (e.g., from 
  // belief prop) to make computing the linear index easier
  std::vector<std::pair<size_t, double> > data;
  {
  variable_t v0(0, 4);
  variable_t v1(1, 3);
  variable_t v2(2, 2);

  std::vector<variable_t> args;
  args.push_back(v0);
  args.push_back(v1);
  args.push_back(v2);
  domain_t dom(args); 

  create_data_vector(data, dom);
  }

  variable_t v0(v0_id, 4);
  variable_t v1(v1_id, 3);
  variable_t v2(v2_id, 2);

  std::vector<variable_t> args;
  args.push_back(v0);
  args.push_back(v1);
  args.push_back(v2);

  sparse_table_t st(args, data);

  return st;
}


void multiplyTest(unsigned v0_id, unsigned v1_id, unsigned v2_id, dense_table_t& dt) 
{
  dense_table_t  dt_gm = dt;
  sparse_table_t st = create_sparse_table(v0_id, v1_id, v2_id);
  sparse_table_t st_gm = st;
  //std::cout << "st " << st << std::endl;

  dt *= dt;

  st *= st;
  //std::cout << "st *= st " << st << std::endl;

  {
  dense_table_t st_as_dt(st.domain());
  st.as_dense_table(st_as_dt); 
  double err = dt.l1_diff(st_as_dt);
  //std::cout << "err: " << err << std::endl;
  ASSERT_LT(err, 1e-4);
  }


  dt = dt_gm;
  st = st_gm;
  unsigned msg_length = st.var(st.domain().var_location(v0_id)).size();
  variable_t v0(v0_id, msg_length);
  dense_table_t msg(v0);

  domain_t::const_iterator asg = msg.domain().begin();
  domain_t::const_iterator end = msg.domain().end();
  for( ; asg != end; ++asg) {
    msg.set_logP( *asg, -1*(rand() % 100) );
  }

  dt *= msg;

  st *= msg;
  //std::cout << "st *= msg " << st << std::endl;

  {
  dense_table_t st_as_dt(st.domain());
  st.as_dense_table(st_as_dt);
  double err = dt.l1_diff(st_as_dt);
  ASSERT_LT(err, 1e-4);
  }
}

int main() {
  // create a table 
  sparse_table_t st_gm = create_sparse_table(2, 0, 1);
  //std::cout << "st_gm " << st_gm << std::endl;  
  dense_table_t dt(st_gm.domain());
  st_gm.as_dense_table(dt);

  // equals test
  sparse_table_t st = create_sparse_table(2, 0, 1);
  //std::cout << "st " << st << std::endl;
  ASSERT_TRUE(st == st_gm);

  // copy test
  sparse_table_t st_copy;
  st_copy = st;
  ASSERT_TRUE(st == st_copy);
  
  // test sparse table data reorder
  testDataReorder(2, 3, 4);
  testDataReorder(2, 4, 3);
  testDataReorder(3, 2, 4);
  testDataReorder(3, 4, 2);
  testDataReorder(4, 2, 3);
  testDataReorder(4, 3, 2);

  // multiply test - compare a dense table to the sparse table. 
  // (i've already compared a pre-ordered dense table to the 
  // re-ordered ones.)
  multiplyTest(2, 0, 1, dt);

  std::cout << "All tests passed" << std::endl;
}
