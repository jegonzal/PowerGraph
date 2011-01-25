#include <iostream>
#include <vector>



#include <graphlab.hpp>

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

