

/**
  * Main method for starting up an application in graphlab
  */

#include "aapo_lassoapp.hpp"
#include <boost/program_options.hpp>
#include <graphlab/distributed/distributed_control.hpp>


int main(int argc,  char *argv[]) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "aapolasso starting\n");
  
  std::string mode;
  std::string disable_objcalc;
  double obj_termination_threshold;
  int K;
  
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
  desc("Shotgun parameters");
  // Set the program options
  desc.add_options()
  ("mode",  boost_po::value<std::string>(&(mode))->default_value("accurate"),
   "Mode = (accurate,cs)")
  ("obj_termination_threshold",  boost_po::value<double>(&(obj_termination_threshold))->default_value(1e-4),
   "obj_termination_threshold = (default = 1e-4)")
   ("K",  boost_po::value<int>(&(K))->default_value(0),
   "K (number of lambda steps) = (default nx/2000)")
  ("disable_objcalc",  boost_po::value<std::string>(&(disable_objcalc))->default_value("true"),
   "Mode = (true,false)");

  boost_po::variables_map vm;
  boost_po::store(boost_po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
  boost_po::notify(vm);
  
  
  std::cout << "=== Initializing distributed Shotgun ===" << std::endl;
  graphlab::distributed_control dc(&argc, &argv);
  std::cout << "=== Entering init_message_processing" << std::endl;

  dc.init_message_processing(1);
  dc.barrier();
  
  aapolassoapp * app = new aapolassoapp(std::string(argv[1]), 0, &dc, mode, disable_objcalc, 
  	obj_termination_threshold, K);  
  app->parse_args(argc, argv);
  app->start();
 } 
