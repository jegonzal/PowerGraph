#ifndef GRAPHLAB_COMMAND_LINE_OPTIONS
#define GRAPHLAB_COMMAND_LINE_OPTIONS

#include <string>


#include <boost/program_options.hpp>


#include <graphlab/engine/engine_options.hpp>


namespace boost {
  
  template<typename T>
  std::string graphlab_vec_to_string(const std::vector<T>& vec) {
    std::stringstream strm;
    strm << "{" ;
    for(size_t i = 0; i < vec.size(); ++i) {
      strm << vec[i];
      if(i < vec.size() - 1) strm << ", ";
    }
    strm << "}";
    return strm.str();
  }
  
  template<>
  std::string lexical_cast<std::string>(const std::vector<size_t>& vec);
  template<>
  std::string lexical_cast< std::string>(const std::vector<int>& vec);
  template<>
  std::string lexical_cast< std::string >(const std::vector<double>& vec);
  template<>
  std::string lexical_cast< std::string>(const std::vector<float>& vec);
  template<>
  std::string lexical_cast< std::string>(const std::vector<std::string>& vec);
};


namespace graphlab {
  
  class command_line_options : public engine_options {

    boost::program_options::options_description desc;
    boost::program_options::positional_options_description 
    pos_opts;
    boost::program_options::variables_map vm;
    
  public:
    command_line_options(const std::string& desc_str = "GraphLab program.",
                         size_t default_ncpus = 2,
                         const std::string default_engine = "async",
                         const std::string default_scope = "edge",
                         const std::string default_scheduler = "fifo") : desc(desc_str) {
      ncpus = default_ncpus;
      engine_type = default_engine;
      scope_type = default_scope;
      scheduler_type = default_scheduler;
      // Add documentation for help
      namespace boost_po = boost::program_options;
      
      desc.add_options()("help", "Print this help message.");
      
      desc.add_options()("schedhelp",
                         boost_po::value<std::string>()->implicit_value(""),
                        "Display help for a particular scheduler.");

    } // End constructor


    
    //! Print the description
    inline void print_description() const { std::cout << desc << std::endl; }


    //! Parse the command line arguments
    bool parse(int argc, char** argv);

    //! test whether the program option was set
    bool is_set(const std::string& option);


    //! Attach a user define option to parse
    template<typename T>
    void attach_option(const std::string& option,
                       T* ret_cont,
                       const std::string& description) {
      namespace boost_po = boost::program_options;
      assert(ret_cont != NULL);
      desc.add_options()
        (option.c_str(), 
         boost_po::value<T>(ret_cont), 
         description.c_str());
    }


    //! attach a user define option to parse
    template<typename T>
    void attach_option(const std::string& option,
                       T* ret_cont,
                       const T& default_value, 
                       const std::string& description) {
      namespace boost_po = boost::program_options;
      assert(ret_cont != NULL);
      desc.add_options()
        (option.c_str(),
         boost_po::value<T>(ret_cont)->default_value(default_value),
         description.c_str());
    }
    
    //! add a positional argument
    void add_positional(const std::string& str);

    
  }; // end class command line options


  




}; // end namespace graphlab


#endif
