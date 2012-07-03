#define OUTPUTLEVEL   LOG_INFO

#include <graphlab.hpp>
#include <v8.h>

#include "cvv8.hpp"
#include "pilot.hpp"

using namespace v8;
using namespace graphlab;
namespace cv = cvv8;

/**
 * Initializes v8 shell and setups bindings.
 */
static int v8_main(const std::string& script){

  cv::Shell shell;
  shell.SetupDefaultBindings();
  HandleScope scope;
  pilot::setup_bindings(shell.Global());
  
  std::cout << "pre-setup bindings" << std::endl;  
  if (script.empty()){
    // if no script provided, read from STDIN
    shell.ExecuteStream(std::cin, "standard in");
  }else {
    logstream(LOG_EMPH) << "Flying script [" << script << "]" << std::endl;
    shell.ExecuteFile(script.c_str());
  }

  return EXIT_SUCCESS;

}

int main(int argc, char **argv){

  global_logger().set_log_level(LOG_DEBUG);
  graphlab::command_line_options clopts("GraphLab Javascript Shell");
  std::string script;
  clopts.attach_option("script", script,
                       "The javascript file.  If none is provided, then I will read "
                       "the script from STDIN");
 
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  
  // graphlab - build communication layers
  mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  pilot_dc = &dc;
  pilot::set_clopts(clopts);

  int rc = EXIT_FAILURE;
  try {
    logstream(LOG_INFO) << "Initializing pilot ..." << std::endl;
    rc = v8_main(script);
    logstream(LOG_INFO) << "Launch completed." << std::endl;
  } catch (std::exception const & e){
    logstream(LOG_ERROR) << "Exception : " << e.what() << std::endl;
    logstream(LOG_ERROR) << "Launch failed." << std::endl;
  }
  
  // graphlab - teardown communication layers
  mpi_tools::finalize();
  logstream(LOG_INFO) << "Bye." << std::endl;
  
  return rc;

}
