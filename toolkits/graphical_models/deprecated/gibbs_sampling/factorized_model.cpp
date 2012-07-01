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


#include "factorized_model.hpp"




#include <graphlab/macros_def.hpp>

void factorized_model::reserve(const size_t num_factors) {
  _factors.reserve(num_factors);
}

factor_t& factorized_model::add_factor(const factor_t& factor) {
  _factors.push_back(factor);
  // // normalize the factor
  // _factors.back().normalize();
  factor_id_t factor_id = _factors.size() - 1;
  for(size_t i = 0; i < factor.num_vars(); ++i) {
    variable_t var = factor.args().var(i); 
    _variables.insert(var);
    // add factor to reverse map
    _var_to_factor[var].insert(factor_id);
  }
  return _factors.back();
}

factor_t& factorized_model::add_factor(const domain_t& vars) {
  _factors.push_back(factor_t());
  factor_t& factor = _factors.back();
  factor.set_args(vars);
  // // normalize the factor
  // _factors.back().normalize();
  factor_id_t factor_id = _factors.size() - 1;
  for(size_t i = 0; i < factor.num_vars(); ++i) {
    variable_t var = factor.args().var(i); 
    _variables.insert(var);
    // add factor to reverse map
    _var_to_factor[var].insert(factor_id);
  }
  return _factors.back();
}


const std::set<factor_id_t>& 
factorized_model::factor_ids(const variable_t& var) const {
  typedef std::map<variable_t, std::set<factor_id_t> >::const_iterator iterator;
  iterator iter = _var_to_factor.find(var);
  ASSERT_TRUE(iter != _var_to_factor.end());
  return iter->second;
}




void factorized_model::load_alchemy(const std::string& filename) {
  // Open an input file stream
  std::ifstream fin(filename.c_str());
  ASSERT_TRUE(fin.good());
  std::string line;
  size_t line_number = 0;
  // Read the first line which should be "variable:"
  const bool success = getline(fin,line,line_number++);
  ASSERT_TRUE(success);
  line = trim(line);
  {
    const std::string variables_str("variables:");
    ASSERT_EQ(line, variables_str);
  }
  // Read all the variables and create a map from the variable name
  // (string) to the variable* prl variable pointer.
  typedef std::map<std::string, variable_t> var_map_type;
  typedef var_map_type::iterator var_map_iter_type;
  var_map_type var_map;
  size_t unique_var_id = 0;
  while(fin.good() &&
        getline(fin, line, line_number++) &&
        trim(line) != "factors:") {
    // Separate into name and size
    line = trim(line);
    ASSERT_GT(line.length(), 0);
    size_t namelen = line.find_last_of('\t');
    size_t varsize = 2;
    // if their is a '\t' character then the variable size follows it
    if(namelen != std::string::npos) {
      std::stringstream istrm(trim(line.substr(namelen)));
      istrm >> varsize;
    }
    // Get the variable name
    std::string var_name = trim(line.substr(0, namelen));
    ASSERT_GT(varsize, 0);
    // Create a new finite variable in the universe
    variable_t variable(unique_var_id++, varsize);
    // Store the variable in the local variable map
    var_map[var_name] = variable;
    _var_name.push_back(var_name);
    ASSERT_EQ(_var_name.size(), unique_var_id);
  }

  // Starting to read factors
  {
    const std::string factors_string("factors:");
    ASSERT_EQ(trim(line), factors_string);
  }

  while(fin.good() && getline(fin, line, line_number++)) {
    /// if the line is empty skip it
    if(trim(line).length() == 0) continue;
    //      std::cout << "Line: " << line << std::endl;
    // the factor being read may contain the same variable multiple
    // times to account for that, we first read a temporary factors,
    // making every variable unique artificially, and then convert
    // it to the factor we actually need

    // Process the arguments
    size_t end_of_variables = line.find("//")-1;
    std::vector<variable_t> args;
    std::set<variable_t> args_set;

    // Read in all the variables in the factor and store them
    for(size_t i = 0; i < end_of_variables;
        i = line.find_first_of('/', i) + 1) {
      // Read the next variable as a string
      std::string variable_name =
        trim(line.substr(i, line.find_first_of('/',i) - i));
      //        std::cout << "Variable Name: " << variable_name << std::endl;
      // Look up the variable in the variable map
      var_map_iter_type iter = var_map.find(variable_name);
      ASSERT_TRUE(iter != var_map.end());
      variable_t var = iter->second;
      // This argument must be unique
      if(args_set.count(var) > 0) {
        std::cout << "Line Number: " << line_number << std::endl;
        ASSERT_EQ(args_set.count(var), 0);
      }

      args_set.insert(var);
      // Save the arguments read from the file
      args.push_back(var);
    } // end of first pass through variable

      // Construct the arguments (which will remap the domain)
    domain_t domain(args);
    //      std::cout << "domain: " << domain << std::endl;
    // Build the factor
    factor_t factor(domain);
      
    // Now for the tricky part we need an assignment in the original
    // order
    domain_t orig_domain;
    for(size_t i = 0; i < args.size(); ++i) {
      orig_domain += variable_t(i, args[i].size());
    }


    // Advance to the correct location in the line
    std::istringstream tbl_values;
    size_t weightpos = line.find("///");
    if (weightpos == std::string::npos) {
      tbl_values.str(line.substr(line.find("//") + 2));
    } else {
      size_t startpos = line.find("//") + 2;
      tbl_values.str(line.substr(startpos, weightpos - startpos));
    }
      
    // Read in the weights
    for(assignment_t orig_asg = orig_domain.begin();
        orig_asg < orig_domain.end(); ++orig_asg) {
      assignment_t asg(domain);
      // Translate the original assignment into the sorted factor assignment
      for(size_t i = 0; i < domain.num_vars(); ++i) {
        size_t variable_id = args[i].id();
        asg.set_asg(variable_id, orig_asg.asg(i));
      }
      // Read then next value
      ASSERT_TRUE(tbl_values.good());
      double value = 0;
      tbl_values >> value;
      // Values are stored in log form      
      factor.logP(asg.linear_index()) = value;                
    }
    // Save the factor to the factor graph
    add_factor(factor);                
  } // End of outer while loop over factors should be end of file

  ASSERT_FALSE(fin.good());
  fin.close();
} // end of load alchemy



  //! Save the factor to a file
void factorized_model::save(graphlab::oarchive &arc) const {
  arc << _variables
      << _factors
      << _var_to_factor
      << _var_name;
}


//! Load the factor from a file
void factorized_model::load(graphlab::iarchive &arc) {
  arc >> _variables
      >> _factors
      >> _var_to_factor
      >> _var_name;
}  


//! save the alchemy file
void factorized_model::save_alchemy(const std::string& filename) const {
  std::ofstream fout(filename.c_str());
  ASSERT_TRUE(fout.good());
  fout << "variables:" << std::endl;
  foreach(variable_t var, _variables) {
    fout << var.id() << '\t' << var.size() << "\n";
  }
  fout << "factors:" << std::endl;
  foreach(const factor_t& factor, _factors) {
    domain_t domain = factor.args();
    for(size_t i = 0; i < domain.num_vars(); ++i) {
      fout << domain.var(i).id();
      if(i + 1 < domain.num_vars()) fout << " / ";
    }
    fout << " // ";
    for(size_t i = 0; i < factor.size(); ++i) {
      fout << factor.logP(i);
      if(i + 1 < factor.size()) fout << ' ';
    }
    fout << '\n';
  }
  fout.flush();
  fout.close();
}


#include <graphlab/macros_undef.hpp>
