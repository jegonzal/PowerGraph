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
#include <graphlab/util/generics/any.hpp>
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
      try{
        options[opt].intval = boost::lexical_cast<size_t>(val);
      }
      catch(boost::bad_lexical_cast &) {
        options[opt].intval = 0;
      }
      try{
        options[opt].dblval = boost::lexical_cast<double>(val);
      }
      catch(boost::bad_lexical_cast &) {
        options[opt].dblval = 0.0;
      }
      options[opt].anyval = val;
    }

    /**
     * Add an option -> value pair where value is a string.
     * There are two version of this function implemented since
     * The any cannot store functions, but only function pointers,
     * we need to be able to differentiate between them
     */
    template <typename T>
    typename boost::disable_if_c<boost::is_function<T>::value, void>::type
    add_option(const std::string &opt, const T &val) {
      if (boost::is_convertible<T, std::string>::value) {
        add_option_str(opt, robust_cast<std::string>(val));
      }
      else {
        options[opt].strval = robust_cast<std::string>(val);
        options[opt].intval = robust_cast<size_t>(val);
        options[opt].dblval = robust_cast<double>(val);
        options[opt].anyval = val;
      }
    }

    template <typename T>
    typename boost::enable_if_c<boost::is_function<T>::value, void>::type
    add_option(const std::string &opt, const T &val) {
      if (boost::is_convertible<T, std::string>::value) {
        add_option_str(opt, robust_cast<std::string>(val));
      }
      else {
        options[opt].strval = robust_cast<std::string>(val);
        options[opt].intval = robust_cast<size_t>(val);
        options[opt].dblval = robust_cast<double>(val);
        options[opt].anyval = &val;
      }
    }






    /**
     * Reads a string option
     */
    inline bool get_string_option(const std::string& opt, 
                                  std::string& val) const {
      std::map<std::string, option_values>::const_iterator i = 
        options.find(opt);
      if (i == options.end()) return false;
      val = i->second.strval;
      return true;
    }

    /**
     * Reads a integer option
     */
    template <typename IntType>
    inline bool get_int_option(const std::string& opt, 
                               IntType& val) const {
      std::map<std::string, option_values>::const_iterator i = 
        options.find(opt);
      if (i == options.end()) return false;
      val = i->second.intval;
      return true;
    }

    /**
     * Reads a float option
     */
    template <typename FloatType>
    inline bool get_float_option(const std::string& opt, 
                                 FloatType& val) const {
      std::map<std::string, option_values>::const_iterator i = 
        options.find(opt);
      if (i == options.end()) return false;
      val = i->second.dblval;
      return true;
    }

    /**
     * Reads an any option
     */
    inline bool get_any_option(const std::string& opt, any &val) const {
      std::map<std::string, option_values>::const_iterator i = 
        options.find(opt);
      if (i == options.end()) return false;
      val = i->second.anyval;
      return true;
    }

    /**
     * Erases an option
     */
    template <typename T>
    inline void erase_option(const std::string &opt) {
      options.erase(opt);
    }

    /**
     * Combines two options. Warns if options intersects
     */
    inline void apply_options(const options_map& other) {
      std::set<std::string> commonopts = 
        set_intersect(keys(options),  keys(other.options));
      if (commonopts.size() > 0) {
        std::set<std::string>::const_iterator i = commonopts.begin();
        while(i != commonopts.end()) {
          std::string myval, otherval;
          this->get_string_option(*i, myval);
          other.get_string_option(*i, otherval);
          if (myval != otherval) {
            logger(LOG_WARNING,
                   "Common options detected between options set\n"
                   "programmatically and options set on the command line.\n");
            logstream(LOG_WARNING) << "\t" << *i << " = " << myval << " || " 
                                   << otherval << "\n";
          }
          ++i;
        }
      }
      std::map<std::string, option_values> res;
      res = map_union(other.options, options);
      options = res;
    }


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
      size_t intval;
      double dblval;
      any anyval;
    };


    
    /**
     * Parse the scheduler string returning the scheduler and storing
     * all the options.
     */
    std::string parse_string(std::string options_raw);

    std::map<std::string, option_values> options;

  };


  std::ostream& operator<<(std::ostream& out,
                           const graphlab::options_map& opts);


} // end of graphlab namespace




#endif

