#include <graphlab/util/cuckoo_map.hpp>
#include <graphlab/util/cuckoo_map_pow2.hpp>
#include <graphlab/util/timer.hpp>
#include <boost/unordered_map.hpp>

void sanity_checks() {
  boost::unordered_map<size_t, size_t> um;
  graphlab::cuckoo_map<size_t, size_t, (size_t)(-1)> cm;
  for (size_t i = 0;i < 10000; ++i) {
    cm[17 * i] = i;
    um[17 * i] = i;
  }

  for (size_t i = 0;i < 10000; ++i) {
    assert(cm[17 * i] == i);
    assert(um[17 * i] == i);
  }
  assert(cm.size() == 10000);
  assert(um.size() == 10000);

  for (size_t i = 0;i < 10000; i+=2) {
    cm.erase(17*i);
    um.erase(17*i);
  }
  for (size_t i = 0;i < 10000; i+=2) {
    assert(cm.count(17*i) == i % 2);
    assert(um.count(17*i) == i % 2);
    if (cm.count(17*i)) {
      assert(cm.find(17*i)->second == i);
    }
  }

  assert(cm.size() == 5000);
  assert(um.size() == 5000);
}

void benchmark() {
  boost::unordered_map<size_t, size_t> um;
  graphlab::cuckoo_map_pow2<size_t, size_t, (size_t)(-1), 4, uint32_t> cm(32);
  //graphlab::cuckoo_map<size_t, size_t, (size_t)(-1), 2, uint32_t> cm;

  graphlab::timer ti;
  ti.start();
  for (size_t i = 0;i < 10000000; ++i) {
    um[17 * i] = i;
  }
  std::cout << "10M unordered map inserts in " << ti.current_time() << " (Load factor = " << um.load_factor() << ")" << std::endl;

  ti.start();
  for (size_t i = 0;i < 10000000; ++i) {
    size_t t = um[17 * i];
    assert(t == i);
  }
  std::cout << "10M unordered map successful probes in " << ti.current_time() << std::endl;


  ti.start();
  for (size_t i = 10000000;i < 20000000; ++i) {
    assert(um.count(17 * i) == 0);
  }
  std::cout << "10M unordered map failed probes in " << ti.current_time() << std::endl;

  um.clear();


  //cm.reserve(12000000);
  ti.start();
  for (size_t i = 0;i < 10000000; ++i) {
    cm[17 * i] = i;
  }
  std::cout << "10M cuckoo map inserts in " << ti.current_time() << " (Load factor = " << cm.load_factor() << ")" << std::endl;

  ti.start();
  for (size_t i = 0;i < 10000000; ++i) {
    size_t t = cm[17 * i];
    assert(t == i);
  }
  std::cout << "10M cuckoo map successful probes in " << ti.current_time() << std::endl;

  ti.start();
  for (size_t i = 10000000;i < 20000000; ++i) {
    assert(cm.count(17 * i) == 0);
  }
  std::cout << "10M cuckoo map failed probes in " << ti.current_time() << std::endl;
  
}

int main(int argc, char** argv) {
  std::cout << "Basic Sanity Checks... ";
  std::cout.flush();
  sanity_checks();
  std::cout << "Done" << std::endl;


  std::cout << "Running Benchmarks" << std::endl;
  benchmark();

  
  
}