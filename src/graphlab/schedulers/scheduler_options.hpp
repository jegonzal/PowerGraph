#ifndef SCHEDULER_OPTIONS
#define SCHEDULER_OPTIONS
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
     Scheduler options data structure.
     Defines a collection of key->value pairs where the key 
     is a string, and the value is an arbitrary data type.
     The scheduler_options class will invisibly cast between string, integer
     and double data types. For instance, if a scheduler expects a 
  */
  class scheduler_options {
  public:
   
    scheduler_options() {};

    explicit scheduler_options(std::string &s) {
      parse_options(s);
    };

    /**
     * Add an option -> value pair where value is a string.
     * Don't use. add_option() prefered.
     */
    inline void add_option_str(const std::string &opt, 
                               const std::string &val) {
      options[opt].strval = val;
      options[opt].intval = boost::lexical_cast<size_t>(val);
      options[opt].dblval = boost::lexical_cast<double>(val);
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
    inline bool get_string_option(const std::string &opt, 
                                  std::string &val) const {
      std::map<std::string, scheduler_option_values>::const_iterator i = 
        options.find(opt);
      if (i == options.end()) return false;
      val = i->second.strval;
      return true;
    }

    /**
     * Reads a integer option
     */
    template <typename IntType>
    inline bool get_int_option(const std::string &opt, 
                               IntType &val) const {
      std::map<std::string, scheduler_option_values>::const_iterator i = 
        options.find(opt);
      if (i == options.end()) return false;
      val = i->second.intval;
      return true;
    }

    /**
     * Reads a float option
     */
    template <typename FloatType>
    inline bool get_float_option(const std::string &opt, 
                                 FloatType &val) const {
      std::map<std::string, scheduler_option_values>::const_iterator i = 
        options.find(opt);
      if (i == options.end()) return false;
      val = i->second.dblval;
      return true;
    }

    /**
     * Reads an any option
     */
    inline bool get_any_option(const std::string &opt, any &val) const {
      std::map<std::string, scheduler_option_values>::const_iterator i = 
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
     * Combines two scheduler options. Warns if options intersects
     */
    inline void merge_options(const scheduler_options &other) {
      std::set<std::string> commonopts = set_intersect(keys(options),
                                                       keys(other.options));
      if (commonopts.size() > 0) {
        std::set<std::string>::const_iterator i = commonopts.begin();
        while(i != commonopts.end()) {
          std::string myval, otherval;
          this->get_string_option(*i, myval);
          other.get_string_option(*i, otherval);
          if (myval != otherval) {
            logger(LOG_WARNING,
                   "Common scheduler options detected between options set\n"
                   "programmatically and options set on the command line.\n");
            logstream(LOG_WARNING) << "\t" << *i << " = " << myval << " || " 
                                   << otherval << "\n";
          }
          ++i;
        }
      }
      std::map<std::string, scheduler_option_values> res;
      res = map_union(options, other.options);
      options = res;
    }


    /**
     * Parses an option stream  of the form "a=b c=d ..."
     */
    inline void parse_options(const std::string &s) {
      std::stringstream strm(s);
      parse_options(strm);
    }

  
    /**
     * Parses an option stream  of the form "a=b c=d ..."
     */
    inline void parse_options(std::istream &s) {
      std::string opt, value;
      // read till the equal
      while(s.good()) {
        getline(s, opt, '=');
        if (s.bad() || s.eof()) break;
        getline(s, value);
        if (s.bad()) break;
        add_option_str(trim(opt), trim(value));
      }
    }

    /// The internal storage of the scheduler options
    struct scheduler_option_values{
      std::string strval;
      size_t intval;
      double dblval;
      any anyval;
    };


    std::map<std::string, scheduler_option_values> options;

  };



  std::ostream& operator<<(std::ostream& out,
                           const graphlab::scheduler_options& opts);

  std::pair<std::string, scheduler_options> 
  parse_scheduler_string(std::string scheduler_raw);
}
#endif
