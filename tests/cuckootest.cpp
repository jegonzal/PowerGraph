#include <graphlab/util/cuckoo_map.hpp>
#include <graphlab/util/cuckoo_map_pow2.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/memory_info.hpp>
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
  graphlab::timer ti;

  size_t NUM_ELS = 10000000;
  
  std::vector<uint32_t> v;
  uint32_t u = 0;
  for (size_t i = 0;i < NUM_ELS; ++i) {
    v.push_back(u);
    u += 1 + rand() % 8;
  }
  std::random_shuffle(v.begin(), v.end());
  graphlab::memory_info::print_usage();

  {
    boost::unordered_map<uint32_t, uint32_t> um;
    ti.start();
    for (size_t i = 0;i < NUM_ELS; ++i) {
      um[v[i]] = i;
    }
    std::cout <<  NUM_ELS / 1000000 << "M unordered map inserts in " << ti.current_time() << " (Load factor = " << um.load_factor() << ")" << std::endl;

    graphlab::memory_info::print_usage();
    
    ti.start();
    for (size_t i = 0;i < 10000000; ++i) {
      size_t t = um[v[i]];
      // assert(t == i);
    }
    std::cout << "10M unordered map successful probes in " << ti.current_time() << std::endl;
    um.clear();
  }

  {
    graphlab::cuckoo_map<uint32_t, uint32_t, (uint32_t)(-1), 2, uint32_t> cm(1024);

    //cm.reserve(102400);
    ti.start();
    for (size_t i = 0;i < NUM_ELS; ++i) {
      cm[v[i]] = i;
      if (i % 1000000 == 0) std::cout << cm.load_factor() << std::endl;

    }
    std::cout <<  NUM_ELS / 1000000 << "M cuckoo map inserts in " << ti.current_time() << " (Load factor = " << cm.load_factor() << ")" << std::endl;

    graphlab::memory_info::print_usage();

    ti.start();
    for (size_t i = 0;i < 10000000; ++i) {
      size_t t = cm[v[i]];
      // assert(t == i);
    }
    std::cout << "10M cuckoo map successful probes in " << ti.current_time() << std::endl;

  }
  
  {
    graphlab::cuckoo_map_pow2<uint32_t, uint32_t, (size_t)(-1), 3, uint32_t> cm(128);
    
    //cm.reserve(102400);
    ti.start();
    for (size_t i = 0;i < NUM_ELS; ++i) {
      cm[v[i]] = i;
      if (i % 1000000 == 0) std::cout << cm.load_factor() << std::endl;

    }
    std::cout << NUM_ELS / 1000000 << "M cuckoo map pow2 inserts in " << ti.current_time() << " (Load factor = " << cm.load_factor() << ")" << std::endl;

    graphlab::memory_info::print_usage();

    ti.start();
    for (size_t i = 0;i < 10000000; ++i) {
      size_t t = cm[v[i]];
      // assert(t == i);
    }
    std::cout << "10M cuckoo map pow2 successful probes in " << ti.current_time() << std::endl;

  }

}

int main(int argc, char** argv) {
  std::cout << "Basic Sanity Checks... ";
  std::cout.flush();
  sanity_checks();
  std::cout << "Done" << std::endl;


  std::cout << "Running Benchmarks" << std::endl;
  benchmark();



}