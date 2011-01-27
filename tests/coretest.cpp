#include <string>
#include <iostream>

#include <graphlab.hpp>

typedef graphlab::core<int, int> core_type;



typedef core_type::types::graph           graph;
typedef core_type::types::update_function update_function;

int main(int argc, char** argv) {
  
  graphlab::command_line_options clopts;

  std::cout << "Flags: " << clopts.get_compile_flags() << std::endl;
  
  int opt1 = 2;
  std::string opt2 = "nope";
  double opt3 = 3.0;
  std::vector<size_t> numb3rs;
  numb3rs.push_back(1);
  numb3rs.push_back(2);

  clopts.attach_option("test1", &opt1, 123,  "some integer option");
  clopts.attach_option("test2", &opt2, std::string("hello world"), "some string option");
  clopts.attach_option("test3", &opt3, 1.0E-5, "some float option");
  clopts.attach_option("numb3rs", &numb3rs, numb3rs, "some numbers");

  bool success = clopts.parse(argc, argv);
  if(!success) {
    clopts.print_description();
    return EXIT_FAILURE;
  }


  
  core_type glcore;
  
  glcore.set_engine_options(clopts);

  std::cout << "Engine       " << clopts.get_engine_type() << "\n"
            << "Option 1:    " << opt1 << "\n"
            << "Option 2:    " << opt2 << "\n"
            << "Option 3:    " << opt3 << "\n"
            << "Option 4:    " << boost::lexical_cast<std::string>(numb3rs) << std::endl;
  
}
