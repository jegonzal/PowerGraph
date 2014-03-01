#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/options/command_line_options.hpp>
#include <graphlab/rpc/dc.hpp>
int __real_main(int argc, char** argv);
graphlab::command_line_options __glopts("");


int actual_main(int argc, char** argv) {
  graphlab::mpi_tools::init(argc, argv);
  if (!__glopts.parse(argc, argv, true)) return false;
  // rebuild argc argv with unrecognized options
  std::vector<std::string> vs = __glopts.unrecognized();
  char** newargv = new char*[vs.size() + 1];
  newargv[0] = argv[0];
  for (size_t i = 0;i < vs.size(); ++ i) {
    newargv[i + 1] = (char*)(vs[i].c_str());
  }
  int ret = __real_main(vs.size() + 1, newargv);
  if (graphlab::dc_impl::get_last_dc()) delete graphlab::dc_impl::get_last_dc();
  graphlab::mpi_tools::finalize();
  return ret;
}

#if 1
// don't seem to be able to get -wrap main working correctly
int main(int argc, char** argv) {
  return actual_main(argc, argv);
}
#else
int __wrap_main(int argc, char** argv) {
  return actual_main(argc, argv);
}

#endif
 
