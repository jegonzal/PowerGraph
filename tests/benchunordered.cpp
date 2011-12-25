
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstring>

#include <cxxtest/TestSuite.h>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>


#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/timer.hpp>

using namespace graphlab;
std::string randstr() {
  std::string ret;
  for (size_t i = 0;i < 8; ++i) {
    ret += char(rand());
  }
  return ret;
}

int main(int argc, char** argv){
  {
    boost::unordered_map<std::string, size_t> m;
    for (size_t i  =0;i < 10000000; ++i)  m[randstr()] = rand();
    timer ti; ti.start();
    std::ofstream f;
    f.open("test.bin",std::fstream::binary);
    oarchive a(f);
    a << m;
    f.close();
    std::cout << "write in: " << ti.current_time() << "\n";
    
    m.clear();
  }
  {
    boost::unordered_map<std::string, size_t> m;
    timer ti; ti.start();
    std::ifstream f;
    f.open("test.bin",std::fstream::binary);
    iarchive a(f);
    a >> m;
    f.close();
    std::cout << "read in: " << ti.current_time() << "\n";
  }
  return 0;
}

