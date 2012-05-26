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

/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */


#ifndef GRAPHLAB_OPTIONS_MAP_HPP
#define GRAPHLAB_OPTIONS_MAP_HPP
#include <map>
#include <sstream>
#include <ostream>
#include <istream>
#include <boost/lexical_cast.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/generics/robust_cast.hpp>

namespace graphlab {

  /**
     options data structure.  Defines a collection of key->value pairs
     where the key is a string, and the value is an arbitrary data
     type.  The options_map class will invisibly cast between string,
     integer and double data types.
  */
  class options_map {
  public:
   
    options_map() {};

    explicit options_map(std::string &s) {
      parse_options(s);
    };

    /**
     * Add an option -> value pair where value is a string.
     * Don't use. add_option() prefered.
     */
    inline void add_option_str(const std::string &opt, 
                               const std::string &val) {
      options[opt].strval = val;
      try {
        options[opt].intval = boost::lexical_cast<int>(val);
      } catch(boost::bad_lexical_cast& error) {options[opt].intval = 0; }
      try {
        options[opt].dblval = boost::lexical_cast<double>(val);
      } catch(boost::bad_lexical_cast& error) { options[opt].dblval = 0.0; }

      if (val == "true" || val == "TRUE" || 
          val == "yes" || val == "YES" || val == "1") options[opt].boolval = true;
    }

    template <typename T>
    void add_option(const std::string& opt, const T& val) {
      if (boost::is_convertible<T, std::string>::value) {
        add_option_str(opt, robust_cast<std::string>(val));
      } else {
        options[opt].strval  = robust_cast<std::string>(val);
        options[opt].intval  = robust_cast<int>(val);
        options[opt].dblval  = robust_cast<double>(val);
        options[opt].boolval = robust_cast<bool>(val);
      }
    }

    /**
     * Test if the option has been created
     */
    inline bool is_set(const std::string& opt) const { 
      return options.find(opt) != options.end();
    }


    /**
     * Reads a string option
     */
    inline bool get_option(const std::string& opt, std::string& val) const {
      std::map<std::string, option_values>::const_iterator i = options.find(opt);
      if (i == options.end()) return false;
      val = i->second.strval;
      return true;
    }

   /**
     * Reads a string option
     */
    inline bool get_option(const std::string& opt, bool& val) const {
      std::map<std::string, option_values>::const_iterator i = options.find(opt);
      if (i == options.end()) return false;
      val = i->second.boolval;
      return true;
    }


    /**
     * Reads a integer option
     */
    template <typename IntType>
    inline bool get_option(const std::string& opt, IntType& val) const {
      std::map<std::string, option_values>::const_iterator i = options.find(opt);
      if (i == options.end()) return false;
      val = i->second.intval;
      return true;
    }

    /**
     * Reads a float option
     */
    inline bool get_option(const std::string& opt, float& val) const {
      std::map<std::string, option_values>::const_iterator i = options.find(opt);
      if (i == options.end()) return false;
      val = i->second.dblval;
      return true;
    }

    /**
     * Reads a double option
     */
    inline bool get_option(const std::string& opt, double& val) const {
      std::map<std::string, option_values>::const_iterator i = options.find(opt);
      if (i == options.end()) return false;
      val = i->second.dblval;
      return true;
    }


    /**
     * Erases an option
     */
    template <typename T>
    inline void erase_option(const std::string &opt) { options.erase(opt); }


    /**
     * Parses an option stream  of the form "a=b c=d ..."
     */
    inline void parse_options(const std::string& s) {
      std::stringstream strm(s);
      parse_options(strm);
    }

  
    /**
     * Parses an option stream  of the form "a=b c=d ..."
     */
    inline void parse_options(std::istream& s) {
      std::string opt, value;
      // read till the equal
      while(s.good()) {
        getline(s, opt, '=');
        if (s.bad() || s.eof()) break;
        getline(s, value, ' ');
        if (s.bad()) break;
        add_option_str(trim(opt), trim(value));
      }
    }

    /// The internal storage of the options
    struct option_values{
      std::string strval;
      int intval;
      double dblval;
      bool boolval;
      option_values () : intval(0), dblval(0), boolval(false) { }
    };


    
    /**
     * Parse a comma delimited series of key1=value1,key2=value2 
     */
    void parse_string(std::string arguments);

    std::map<std::string, option_values> options;

  };


  std::ostream& operator<<(std::ostream& out,
                           const graphlab::options_map& opts);


} // end of graphlab namespace




#endif

