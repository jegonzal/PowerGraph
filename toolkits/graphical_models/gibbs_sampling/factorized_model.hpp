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


#ifndef PGIBBS_FACTORIZED_MODEL_HPP
#define PGIBBS_FACTORIZED_MODEL_HPP

/**
 *
 * Represents a factorized model in alchemy format
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>


#include <iostream>
#include <iomanip>

#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cassert>



#include <graphlab.hpp>




// The maximum number of dimensions in a factor table
const size_t MAX_TREEWIDTH = 32;
const size_t MAX_DIM = MAX_TREEWIDTH + 1;

// Basic graphical model typedefs
typedef uint32_t                        factor_id_t;
typedef graphlab::discrete_variable     variable_t;
typedef graphlab::table_factor<MAX_DIM> factor_t;
typedef factor_t::domain_type           domain_t;
typedef factor_t::assignment_type       assignment_t;



// A class which represents a factorized distribution as a collection
// of factors.
class factorized_model {
private:
  std::set<variable_t> _variables;
  std::vector<factor_t> _factors;
  std::map< variable_t, std::set<factor_id_t> > _var_to_factor;
  std::vector<std::string> _var_name;

  /**
   * same as the stl string get line but this also increments a line
   * counter which is useful for debugging purposes
   */
  inline bool getline(std::ifstream& fin,
                      std::string& line,
                      size_t line_number) {
    return std::getline(fin, line).good();
  }

  /**
   * Removes trailing and leading white space from a string
   */
  inline std::string trim(const std::string& str) {
    std::string::size_type pos1 = str.find_first_not_of(" \t\r");
    std::string::size_type pos2 = str.find_last_not_of(" \t\r");
    return str.substr(pos1 == std::string::npos ? 0 : pos1,
                      pos2 == std::string::npos ? str.size()-1 : pos2-pos1+1);
  }

public:

  typedef std::vector<factor_t> factor_map_t;


  void reserve(const size_t num_factors);

  //! add a factor to the factorized model
  factor_t& add_factor(const factor_t& factor);

  //! add a factor to the factorized model
  factor_t& add_factor(const domain_t& dom);

  
  const factor_map_t& factors() const { return _factors; }
  const std::set<variable_t>& variables() const { return _variables; }

  const std::set<factor_id_t>& factor_ids(const variable_t& var) const;

  const std::string& var_name(size_t id) const {
    ASSERT_LT(id, _var_name.size());
    return _var_name[id];
  }

  
  void load_alchemy(const std::string& filename);

  //! Save the factor to a file
  void save(graphlab::oarchive &arc) const;

  //! Load the factor from a file
  void load(graphlab::iarchive &arc);

  //! save the alchemy file
  void save_alchemy(const std::string& filename) const;
 
}; //end of factorized model


#endif



