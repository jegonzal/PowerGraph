/*
 *  pagerank.cpp
 *  GraphLab_saved
 *
 */

#include <string>
#include <boost/program_options.hpp>

#include "distpagerankapp.hpp"

int main(int argc,  char *argv[]) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  
  
   // Initialize distributed stuff
  distributed_control dc(&argc, &argv);
  dc.init_message_processing(1);
 
  logger(LOG_INFO, "PageRank starting\n");
  
  if (argc == 1) {
    logger(LOG_ERROR, "Usage: eigenvector --infile=[filename-csv]\n"
                      "                   [--binoutfile=output binary graph file]\n"
                      "                   [other options]\n");
    return 1;
  }

  std::string filename;
  std::string binoutfile;
  bool optimize = false;
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
  desc("PageRank input file. Can be a Jure text file, or a graphlab serialized binary file");
  // Set the program options
  desc.add_options()
  ("infile",  boost_po::value<std::string>(&(filename))->default_value(""),
   "Input filename ")
   ("binoutfile",  boost_po::value<std::string>(&(binoutfile))->default_value(""),
   "Optional Output binary graph")
   ("optimize", boost_po::value<bool>(&optimize)->default_value(false), "Optimize Graph Layout");

  boost_po::variables_map vm;
  boost_po::store(boost_po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
  boost_po::notify(vm);

  if (filename == "") {
    std::cout << "Input Graph needed\n";
    return 0;
  }

  pagerankapp app(&dc, filename, binoutfile, optimize);
  
  app.parse_args(argc, argv);
  app.start();
}
