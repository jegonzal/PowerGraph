/*
 * \author akyrola
 *  Starter for CoEM
 */

#include <cstdio>
#include <string>

#include <graphlab/logger/logger.hpp>

#include "distmulticoemapp.hpp"


int main(int argc,  char *argv[]) {
  global_logger().set_log_level(0);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "Distributed MultiCoEM / ReadTheWeb starting\n");
  
  std::string ROOT = "/mnt/bigbrofs/usr5/graphlab/testdata/coem/justin";
  
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
  desc(" ");
  // Set the program options
  desc.add_options()
  ("root",  boost_po::value<std::string>(&(ROOT))->default_value(ROOT),
   "Input directory ");

  boost_po::variables_map vm;
  boost_po::store(boost_po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
  boost_po::notify(vm);

  
  
  std::cout << "root = " << ROOT << std::endl;
   // Initialize distributed stuff
  distributed_control dc(&argc, &argv);
  dc.init_message_processing(1);
  
  distmulticoemapp * app = new distmulticoemapp(&dc, 
    ROOT + "/cat_nps.txt",
  	ROOT + "/cat_contexts.txt", 
  	ROOT + "/cat_pairs_cont-idx.txt", 
  	ROOT + "/seeds/", 
  	ROOT + "/seeds-neg/");
  app->parse_args(argc, (char**) argv);
  app->start();
}
