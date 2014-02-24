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
 * Test various functions of a dense table
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

#include <factors/dense_table.hpp>


const size_t MAX_DIM = 4;
typedef graphlab::dense_table<MAX_DIM>          dense_table_t;
typedef graphlab::discrete_domain<MAX_DIM>      domain_t;
typedef graphlab::discrete_assignment<MAX_DIM>  assignment_t;
typedef graphlab::discrete_variable             variable_t;


dense_table_t create_rand_dense_table(unsigned v0_id, unsigned v1_id, unsigned v2_id) 
{
  variable_t v0(v0_id, 4);
  variable_t v1(v1_id, 3);
  variable_t v2(v2_id, 2);

  std::vector<variable_t> vars;
  vars.push_back(v0);
  vars.push_back(v1);
  vars.push_back(v2);
  domain_t domain(vars); 

  dense_table_t dt(domain);
  assignment_t da(domain);

  for(int i=0; i < domain.var(0).size(); ++i) {
    da.set_asg(v0, i);
    for(int j=0; j < domain.var(1).size(); ++j) {
      da.set_asg(v1, j);
      for(int k=0; k < domain.var(2).size(); ++k) {
        da.set_asg(v2, k);

        if(rand() % 100 <= 20) {
          dt.set_logP(da, -1 * (rand() % 100));
        }
      }
    }
  }

  return dt;
}

std::vector<double> create_rand_data_vector(size_t d0, size_t d1, size_t d2) {
  std::vector<double> v(d0*d1*d2);
  for(size_t i=0; i < v.size(); ++i) {
    v[i] = -1 * (rand() % 100);
  }
  return v;
}

// [-Inf -Inf -3   -Inf]
// [-15  -4   -Inf -23 ]
// [-Inf -Inf -Inf -Inf]
//
// [-Inf -Inf -Inf -Inf]
// [-20  -12  -19  -78 ]
// [-Inf -Inf -Inf -32 ]
void create_data_vector(std::vector<double> &data) 
{
  data.resize(4*3*2);
  data[0] = -1000; data[1] = -1000; data[2]  = -3;    data[3]  = -1000;
  data[4] = -15;   data[5] = -4;    data[6]  = -1000; data[7]  = -23;
  data[8] = -1000; data[9] = -1000; data[10] = -1000; data[11] = -1000;

  data[12] = -1000; data[13] = -1000; data[14] = -1000; data[15] = -1000;
  data[16] = -20;   data[17] = -12;   data[18] = -19;   data[19] = -78;
  data[20] = -1000; data[21] = -1000; data[22] = -1000; data[23] = -32;
}

void testDataReorder(unsigned v0_id, unsigned v1_id, unsigned v2_id) {
  // setup the data as if its from belief prop
  std::vector<double> data = create_rand_data_vector(4, 4, 4);

  // create a table with variables that are not resorted
  variable_t v0(0, 4);
  variable_t v1(1, 4);
  variable_t v2(2, 4);
  std::vector<variable_t> vars;
  vars.push_back(v0);
  vars.push_back(v1);
  vars.push_back(v2);
  dense_table_t dt_gm(vars, data);
  //std::cout << "dt_gm " << dt_gm << std::endl;

  //std::cout << "v0_id " << v0_id << " v1_id " << v1_id << " v2_id " << v2_id << std::endl;
  variable_t r0(v0_id, 4);
  variable_t r1(v1_id, 4);
  variable_t r2(v2_id, 4);
  std::vector<variable_t> reordered_vars;
  reordered_vars.push_back(r0);
  reordered_vars.push_back(r1);
  reordered_vars.push_back(r2);
  dense_table_t dt(reordered_vars, data);
  //std::cout << "dt " << dt << std::endl;


  for(int i=0; i < dt.numel(); ++i) {
    assignment_t da_gm(domain_t(vars), i);
    // permute assignment
    std::vector<size_t> asgs(da_gm.begin(), da_gm.end());
    assignment_t da(reordered_vars, asgs);
    
    //std::cout << "dt_gm{" << da_gm << "}=" << dt_gm.logP(da_gm) << "; "
    //          << "dt{" << da << "}=" << dt.logP(da) << std::endl;
    ASSERT_TRUE(dt_gm.logP(da_gm) == dt.logP(da));
  }
}

dense_table_t create_dense_table(unsigned v0_id, unsigned v1_id, unsigned v2_id) {
  // create a table with dimensions ordered like the original table (e.g., from 
  // belief prop) to make computing the linear index easier
  std::vector<double> data;
  create_data_vector(data);

  variable_t v0(v0_id, 4);
  variable_t v1(v1_id, 3);
  variable_t v2(v2_id, 2);

  std::vector<variable_t> vars;
  vars.push_back(v0);
  vars.push_back(v1);
  vars.push_back(v2);

  dense_table_t dt(vars, data);

  return dt;
}


void multiplyTest(unsigned v0_id, unsigned v1_id, unsigned v2_id) 
{
  dense_table_t dt = create_dense_table(v0_id, v1_id, v2_id);
  dense_table_t dt_gm = dt;
  //std::cout << "dt " << dt << std::endl;

  dt *= dt;
  //std::cout << "dt *= dt " << dt << std::endl;
  {
  domain_t::const_iterator asg = dt.domain().begin();
  domain_t::const_iterator end = dt.domain().end();
  for( ; asg != end; ++asg)
    ASSERT_EQ(dt.logP(*asg), dt_gm.logP(*asg)+dt_gm.logP(*asg));
  }

  dt = dt_gm;
  unsigned msg_length = dt.var(dt.domain().var_location(v0_id)).size();
  variable_t v0(v0_id, msg_length);
  dense_table_t msg(v0);
  domain_t::const_iterator asg = msg.domain().begin();
  domain_t::const_iterator end = msg.domain().end();
  for( ; asg != end; ++asg) {
    msg.set_logP( *asg, -1*(rand() % 100) );
  }

  dt *= msg;
  //std::cout << "msg = " << msg << std::endl;
  //std::cout << "dt *= msg " << dt << std::endl;

  for(size_t i=0; i < dt.size(); ++i) {
    assignment_t dt_asg(dt.domain(), i);
    assignment_t msg_asg = dt_asg.restrict(msg.domain());
    ASSERT_EQ(dt.logP(dt_asg), dt_gm.logP(dt_asg)+msg.logP(msg_asg));
  }
}

int main() {
  // create a table 
  dense_table_t dt_gm = create_dense_table(2, 0, 1);
  //std::cout << "dt_gm " << dt_gm << std::endl;  

  double err;

  // equals test
  dense_table_t dt = create_dense_table(2, 0, 1);
  //std::cout << "dt " << dt << std::endl;
  err = dt.l1_diff(dt_gm);
  //std::cout << "err: " << err << std::endl;
  ASSERT_LT(err, 1e-4);

  // copy test
  dense_table_t dt_copy;
  dt_copy = dt;
  err = dt.l1_diff(dt_copy);
  //std::cout << "err: " << err << std::endl;
  ASSERT_LT(err, 1e-4);

  // test dense table data reorder
  testDataReorder(2, 3, 4);
  testDataReorder(2, 4, 3);
  testDataReorder(3, 2, 4);
  testDataReorder(3, 4, 2);
  testDataReorder(4, 2, 3);
  testDataReorder(4, 3, 2);

  // multiply test - compare a (pre-ordered) dense table to the dense table 
  multiplyTest(2, 3, 4);
  multiplyTest(2, 4, 3);
  multiplyTest(3, 2, 4);
  multiplyTest(3, 4, 2);
  multiplyTest(4, 2, 3);
  multiplyTest(4, 3, 2);

  std::cout << "All tests passed" << std::endl;
}
