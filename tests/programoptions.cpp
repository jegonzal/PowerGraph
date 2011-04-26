#include <iostream>
#include <boost/program_options.hpp>

// This tests if it is possible to make pass argc and argv through 2 sets of
// program options

struct options{
  size_t a;
  size_t b;
  size_t c;
};

void parse1(int argc, char** argv, options& opts) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
  desc("Denoise a randomly generated image using Gibbs Sampling.");
  // Set the program options
  desc.add_options()
  ("a",  boost_po::value<size_t>(&(opts.a))->default_value(1),
   "a")
  ("b",  boost_po::value<size_t>(&(opts.b))->default_value(2),
   "b");

  boost_po::variables_map vm;
  boost_po::store(boost_po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
  boost_po::notify(vm);
}

void parse2(int argc, char** argv, options& opts) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
  desc("Denoise a randomly generated image using Gibbs Sampling.");
  // Set the program options
  desc.add_options()
  ("c",  boost_po::value<size_t>(&(opts.c))->default_value(3),
   "c");

  boost_po::variables_map vm;
  boost_po::store(boost_po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
  boost_po::notify(vm);
}

int main(int argc, char** argv) {
  options opt;
  parse1(argc, argv, opt);
  parse2(argc, argv, opt);
  std::cout << "a=" << opt.a << "\n";
  std::cout << "b=" << opt.b << "\n";
  std::cout << "c=" << opt.c << "\n";
}

