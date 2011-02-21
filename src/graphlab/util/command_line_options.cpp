#include <graphlab/util/command_line_options.hpp>
#include <graphlab/schedulers/scheduler_list.hpp>

namespace boost {  
  template<>
  std::string lexical_cast<std::string>(const std::vector<size_t>& vec) {
    return graphlab_vec_to_string(vec);
  }

  template<>
  std::string lexical_cast< std::string>(const std::vector<int>& vec) {
    return graphlab_vec_to_string(vec);
  }

  template<>
  std::string lexical_cast< std::string >(const std::vector<double>& vec) {
    return graphlab_vec_to_string(vec);
  }

  template<>
  std::string lexical_cast< std::string>(const std::vector<float>& vec) {
    return graphlab_vec_to_string(vec);
  }

  template<>
  std::string lexical_cast< std::string>(const std::vector<std::string>& vec) {
    return graphlab_vec_to_string(vec);
  }
};


namespace graphlab {
  
  bool command_line_options::parse(int argc, char** argv) {
    namespace boost_po = boost::program_options;
    
    size_t ncpus(get_ncpus());
    std::string enginetype(get_engine_type());
    bool cpuaffin(get_cpu_affinities());
    bool schedyield(get_sched_yield());
    std::string scopetype(get_scope_type());
    std::string schedulertype(get_scheduler_type());
    std::string metricstype(get_metrics_type());

    if(!surpress_graphlab_options) {
      // Set the program options
      desc.add_options()
        ("ncpus",
         boost_po::value<size_t>(&(ncpus))->
         default_value(ncpus),
         "Number of cpus to use.")
        ("engine",
         boost_po::value<std::string>(&(enginetype))->
         default_value(enginetype),
         "Options are {async, async_sim, synchronous}")
        ("affinities",
         boost_po::value<bool>(&(cpuaffin))->
         default_value(cpuaffin),
         "Enable forced assignment of threads to cpus")
        ("schedyield",
         boost_po::value<bool>(&(schedyield))->
         default_value(schedyield),
         "Enable yielding when threads conflict in the scheduler.")
        ("scope",
         boost_po::value<std::string>(&(scopetype))->
         default_value(scopetype),
         "Options are {none, vertex, edge, full}")
        ("metrics",
         boost_po::value<std::string>(&(metricstype))->
         default_value(metricstype),
         "Options are {none, basic, file, html}")
        ("schedhelp",
         boost_po::value<std::string>()->implicit_value(""),
         "Display help for a particular scheduler.")
        ("scheduler",
         boost_po::value<std::string>(&(schedulertype))->
         default_value(schedulertype),
         (std::string("Supported schedulers are: ")
          + get_scheduler_names_str() +
          ". Too see options for each scheduler, run the program with the option"
          " ---schedhelp=[scheduler_name]").c_str());
    }
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
    if (vm.count("schedhelp")) {
      std::string schedname = vm["schedhelp"].as<std::string>();
      if (schedname != "") {
        print_scheduler_info(schedname, std::cout);
      }
      else {
        std::vector<std::string> schednames = get_scheduler_names();
        for(size_t i = 0;i < schednames.size(); ++i) {
          print_scheduler_info(schednames[i], std::cout);
        }
      }
      return false;
    }    

    set_ncpus(ncpus);

    if(!set_engine_type(enginetype)) {
      std::cout << "Invalid engine type! : " << enginetype 
                << std::endl;
      return false;
    }

    set_cpu_affinities(cpuaffin);
    set_sched_yield(schedyield);

    if(!set_scope_type(scopetype)) {
      std::cout << "Invalid scope type! : " << scopetype
                << std::endl;
      return false;
    }

    if(!set_scheduler_type(schedulertype)) {
      std::cout << "Invalid scheduler type! : " << schedulertype 
                << std::endl;
      return false;
    }
    
    if(!set_metrics_type(metricstype)) {
      std::cout << "Invalid metrics type! : " << metricstype
                << std::endl;
      return false;
    }

    return true;
  } // end of parse


  bool command_line_options::is_set(const std::string& option) {
    return vm.count(option);
  }

  void command_line_options::add_positional(const std::string& str) {
    pos_opts.add(str.c_str(), 1);
  }
}

