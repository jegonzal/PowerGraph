#include <cassert>
#include <iostream>

#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/distributed2/graph/atom_index_file.hpp>
#include <graphlab.hpp>


int main(int argc, char** argv) {

  // set the global logger
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Initialize the mpi tools
  graphlab::mpi_tools::init(argc, argv);




  // Parse the command lines
  std::string path;
  std::string aindex("atom_index.txt");
  graphlab::command_line_options 
    clopts("Construct an atom index file from a collection of "
           "potentially distributed files. ", true);
  clopts.attach_option("path", &path, path, 
                       "Location of the atom files.  "
                       "This could be a local path.");
  clopts.attach_option("aindex", &aindex, aindex,
                       "[output] atom index file.");

 

  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }

  graphlab::build_atom_index_file(path, aindex);

  graphlab::mpi_tools::finalize();  
  return EXIT_SUCCESS;
}
