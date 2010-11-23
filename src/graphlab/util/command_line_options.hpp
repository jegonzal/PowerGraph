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
      desc.add_options()("help", "Print this help message.");
    } // End constructor


    
    //! Print the description
    void print_description() const { std::cout << desc << std::endl; }


    //! Parse the command line arguments
    bool parse(int argc, char** argv) {
      namespace boost_po = boost::program_options;
      
      // Set the program options
      desc.add_options()
        ("ncpus",
         boost_po::value<size_t>(&(ncpus))->
         default_value(ncpus),
         "Number of cpus to use.")
        ("engine",
         boost_po::value<std::string>(&(engine_type))->
         default_value(engine_type),
         "Options are {async, async_sim, synchronous}")
        ("affinities",
         boost_po::value<bool>(&(enable_cpu_affinities))->
         default_value(enable_cpu_affinities),
         "Enable forced assignment of threads to cpus")
        ("schedyield",
         boost_po::value<bool>(&(enable_sched_yield))->
         default_value(enable_sched_yield),
         "Enable yeilding when threads conflict in the scheduler.")
        ("scope",
         boost_po::value<std::string>(&(scope_type))->
         default_value(scope_type),
         "Options are {none, vertex, edge, full}")
        ("scheduler",
         boost_po::value<std::string>(&(scheduler_type))->
         default_value(scheduler_type),
         "There are several scheduler supported by the graphlab framework:"
         "{fifo, sweep({permute,linear}), multiqueue_fifo, priority, "
         "sampling,  multiqueue_priority, "
         "clustered_priority(one of {metis,bfs,random}, vertices perpartition), "
         "round_robin(iterations), colored(iterations)}");
      // Parse the arguments
      try{
        boost_po::store(boost_po::command_line_parser(argc, argv).
                        options(desc).positional(pos_opts).run(), vm);
        boost_po::notify(vm);
      } catch( boost_po::error error) {
        std::cout << "Invalid syntax:\n" 
                  << "\t" << error.what()
                  << "\n\n" << std::endl 
                  << "Description:"
                  << std::endl;        
        print_description();
        return false;
      }
      if(vm.count("help")) {
        print_description();
        return false;
      }
      return true;
    } // end of parse


    //! test whether the program option was set
    bool is_set(const std::string& option) {
      return vm.count(option);
    }


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
    void add_positional(const std::string& str) {
      pos_opts.add(str.c_str(), 1);
    }

    
  }; // end class command line options


  




}; // end namespace graphlab


#endif
