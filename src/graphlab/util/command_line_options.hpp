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

}; // end of namespace boost


namespace graphlab {
  /**
  Because of the many engine options GraphLab has relatively
  sophisticated command line parsing tools. However we have found that
  many of our ML applications had poorly written command line parsing
  support so we tried to generalize the GraphLab command line tools to
  enable user applications to benefit from sophisticated and still
  easy to use command line parsing.

  The command_line_options data-structure extends (wrapping) the
  boost::program_options library. We have tried to retain much of the
  functionality of the boost::program_options library while hiding
  some of the less "friendly" template meta-programming "features".
  
  Here is an example of how the library is used:
  
  \code
  int main(int argc, char** argv) {
  
    std::string filename;
    size_t dimensions = 20;
    double bound = 1E-5;
    bool use_x = false;
    std::vector<size_t> nsamples(1,10000);
  
    // Parse command line options
    graphlab::command_line_options clopts("Welcome to a the HelloWorld");
    clopts.attach_option("file", &filename, 
                         "The input filename (required)");
    clopts.add_positional("file");
    clopts.attach_option("dim",
                         &dimensions, dimensions,
                         "the dimension of the grid");
    clopts.attach_option("bound",
                         &bound, bound,
                         "The termination bound");
    clopts.attach_option("usex",
                         &use_x, use_x,
                         "Use algorithm x");
    clopts.attach_option("nsamples",
                         &nsamples, nsamples,
                         "A vector of the number of samples"); 
    clopts.set_scheduler_type("fifo");
    clopts.set_scope_type("edge");
  
    if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  
    if(!clopts.is_set("file")) {
      std::cout << "Input file not provided" << std::endl;
  
      clopts.print_description();
      return EXIT_FAILURE;
    }
  }
  \endcode
  
  
  */
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


    
    /// Print the same message that is printed when the --help command line argument is provided.
    inline void print_description() const { std::cout << desc << std::endl; }


    /**
    \brief This function should be called AFTER all the options have been seen
    (including positionals). The parse function reads the standard command line 
    arguments and fills in the attached variables. If there is an error in the 
    syntax or parsing fails the parse routine will print the error and return false. 
    */
    bool parse(int argc, char** argv);

    /** The is set function is used to test if the user provided the option. 
    The option string should match one of the attached options. 
    */
    bool is_set(const std::string& option);


    /**
    \brief attach a user defined option to the command line options
    parser.
    
    The attach option command is used to attach a user defined option
    to the command line options parser.
    
    \param option The name of the command line flag for that option.

    \param ret_cont A pointer to an "arbitrary" type which can be any
                    of the basic types (char, int, size_t, float,
                    double, bool, string...) or an std::vector of
                    basic types. It is important that the ret_cont
                    point to a memory block that will exist when parse
                    is invoked.
                    
    \param desc Used to describe the option when --help 
          is called or when print_description is invoked.
    */
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


    /**
    \brief attach a user defined option to the command line options
    parser.
    
    The attach option command is used to attach a user defined option
    to the command line options parser.
    
    \param option The name of the command line flag for that option.

    \param ret_cont A pointer to an "arbitrary" type which can be any
                    of the basic types (char, int, size_t, float,
                    double, bool, string...) or an std::vector of
                    basic types. It is important that the ret_cont
                    point to a memory block that will exist when parse
                    is invoked.

    \param default_value The default value of the parameter if the
                         user does not provide this parameter on the
                         command line.

    \param desc Used to describe the option when --help 
          is called or when print_description is invoked.
    */
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
    
    /** This function adds the option as a positional argument.  A
    positional argument does not require --option and instead is read
    based on its location. Each add_positional call adds to the next
    position. */
    void add_positional(const std::string& str);

    
  }; // end class command line options


  




}; // end namespace graphlab


#endif
