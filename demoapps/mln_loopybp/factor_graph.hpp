#ifndef FACTOR_GRAPH_HPP
#define FACTOR_GRAPH_HPP

#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>


#include <graphlab/factors/table_factor.hpp>

class factor_graph {
public:
  static const size_t MAX_DIM = 16;
  typedef graphlab::table_factor<MAX_DIM>  factor_type;
  typedef factor_type::assignment_type     assignment_type;
  typedef factor_type::domain_type         domain_type;
  typedef factor_type::variable_type       variable_type;

  void add_factor(const factor_type& factor) {
    _factors.push_back(factor);
    for(size_t i = 0; i < factor.num_vars(); ++i) {
      _variables.insert(factor.args().var(i));
    }
  }
  
  const std::vector<factor_type>& factors() const { return _factors; }
  const std::set<variable_type>& variables() const { return _variables; }

  const std::string& var_name(size_t id) const {
    assert(id < _var_name.size());
    return _var_name[id];
  }

  
  void load_alchemy(const std::string& filename) {
    // Open an input file stream
    std::ifstream fin(filename.c_str());
    assert(fin.good());
    std::string line;
    size_t line_number = 0;
    // Read the first line which should be "variable:"
    assert(getline(fin,line,line_number++));
    line = trim(line);
    assert(line == "variables:");
    // Read all the variables and create a map from the variable name
    // (string) to the variable* prl variable pointer.
    typedef std::map<std::string, variable_type> var_map_type;
    typedef var_map_type::iterator var_map_iter_type;
    var_map_type var_map;
    size_t unique_var_id = 0;
    while(fin.good() &&
          getline(fin, line, line_number++) &&
          trim(line) != "factors:") {
      // Separate into name and size
      line = trim(line);
      assert(line.length() > 0);
      size_t namelen = line.find_last_of('\t');
      size_t varsize = 2;
      // if their is a '\t' character then the variable size follows it
      if(namelen != std::string::npos) {
        std::stringstream istrm(trim(line.substr(namelen)));
        istrm >> varsize;
      }
      // Get the variable name
      std::string var_name = trim(line.substr(0, namelen));
      // Create a new finite variable in the universe
      variable_type variable(unique_var_id++, varsize);
      // Store the variable in the local variable map
      var_map[var_name] = variable;
      _var_name.push_back(var_name);
      assert(_var_name.size() == unique_var_id);
    }

    // Starting to read factors
    assert(trim(line) == "factors:");

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
      std::vector<variable_type> args;
      std::set<variable_type> args_set;

      // Read in all the variables in the factor and store them
      for(size_t i = 0; i < end_of_variables;
          i = line.find_first_of('/', i) + 1) {
        // Read the next variable as a string
        std::string variable_name =
          trim(line.substr(i, line.find_first_of('/',i) - i));
        //        std::cout << "Variable Name: " << variable_name << std::endl;
        // Look up the variable in the variable map
        var_map_iter_type iter = var_map.find(variable_name);
        assert(iter != var_map.end());
        variable_type var = iter->second;
        // This argument must be unique
        if(args_set.count(var) > 0) {
          std::cout << "Line Number: " << line_number << std::endl;
          assert(args_set.count(var) == 0);
        }

        args_set.insert(var);
        // Save the arguments read from the file
        args.push_back(var);
      } // end of first pass through variable

      // Construct the arguments (which will remap the domain)
      domain_type domain(args);
      //      std::cout << "domain: " << domain << std::endl;
      // Build the factor
      factor_type factor(domain);
      
      // Now for the tricky part we need an assignment in the original
      // order
      domain_type orig_domain;
      for(size_t i = 0; i < args.size(); ++i) {
        orig_domain += variable_type(i, args[i].size());
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
      for(assignment_type orig_asg = orig_domain.begin();
          orig_asg < orig_domain.end(); ++orig_asg) {
        assignment_type asg(domain);
        // Translate the original assignment into the sorted factor assignment
        for(size_t i = 0; i < domain.num_vars(); ++i) {
          size_t variable_id = args[i].id();
          asg.set_asg(variable_id, orig_asg.asg(i));
        }
        // Read then next value
        assert(tbl_values.good());
        double value = 0;
        tbl_values >> value;
        // Values are stored in log form      
        factor.logP(asg.linear_index()) = value;                
      }

      // Save the factor to the factor graph
      add_factor(factor);                
    } // End of outer while loop over factors should be end of file

    assert(fin.good() == false);
    fin.close();
  } // end of load alchemy

  
  
private:
  std::set<variable_type> _variables;
  std::vector<factor_type> _factors;
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

}; //end of factor graph


#endif
