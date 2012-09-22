#include <graphlab/util/hopscotch_table.hpp>
#include <graphlab/util/hopscotch_map.hpp>
#include <graphlab/util/cuckoo_map_pow2.hpp>
#include <boost/unordered_set.hpp>
#include <boost/bind.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/parallel/pthread_tools.hpp>

#include <graphlab/util/memory_info.hpp>
#include <graphlab/macros_def.hpp>

const size_t NINS = 1500000;
boost::unordered_set<uint32_t> um;
graphlab::hopscotch_table<uint32_t> cm(1.2 * NINS);
 
void parallel_inserter(size_t start, size_t end) {
  for (size_t i = start; i < end; ++i) {
    assert(cm.put_sync(17 * i));
  }
}

void parallel_erase(size_t start, size_t end) {
  for (size_t i = start; i < end; i+=2) {
    assert(cm.erase_sync(17 * i));
    assert(cm.put_sync(17 * i));
    assert(cm.get_sync(17 * i).first);
    assert(cm.erase_sync(17 * i));
    assert(cm.get_sync(17 * i).first == false);
    assert(cm.put_sync(17 * i));
    assert(cm.get_sync(17 * i).first);
    assert(cm.erase_sync(17 * i));
    assert(cm.get_sync(17 * i).first == false);
    assert(cm.put_sync(17 * i));
    assert(cm.get_sync(17 * i).first);
    assert(cm.erase_sync(17 * i));
    assert(cm.get_sync(17 * i).first == false);
  }
}



void sanity_checks() {
 ASSERT_TRUE(cm.begin() == cm.end());

  for (size_t i = 0;i < NINS; ++i) {
    um.insert(17 * i);
  }
  
  graphlab::thread_group thrgroup;
  for (size_t i = 0; i < 10; ++i) {
    thrgroup.launch(boost::bind(parallel_inserter, 
                                i * NINS/10, 
                                (i + 1) * NINS / 10));

  }

  thrgroup.join();
   
  std::cout << "Size: " << cm.size() << std::endl;
  std::cout << "Capacity: " << cm.capacity() << std::endl;
  std::cout << "Load Factor: " << cm.load_factor() << std::endl;

  for (size_t i = 0;i < NINS; ++i) {
    if (cm.get_sync(17 * i).first == false) {
      std::cout << cm.count(17 * i) << "\n";
      std::cout << "Failure on: " << 17 * i << std::endl;
      assert(cm.get_sync(17 * i).first == true);
    }
    assert(um.count(17 * i) == 1);
  }
  assert(cm.size() == NINS);
  assert(um.size() == NINS);

  for (size_t i = 0;i < NINS; i+=2) {
    um.erase(17*i);
  }

  for (size_t i = 0; i < 10; ++i) {
    thrgroup.launch(boost::bind(parallel_erase, 
                                i * NINS/10, 
                                (i + 1) * NINS / 10));
  }

  thrgroup.join();
   

  for (size_t i = 0;i < NINS; i+=2) {
    assert(cm.get_sync(17*i).first == (bool)(i % 2));
    assert(um.count(17*i) == i % 2);
  }

  assert(cm.size() == NINS / 2);
  assert(um.size() == NINS / 2);
} 












boost::unordered_map<uint32_t, uint32_t> um2;
graphlab::hopscotch_map<uint32_t, uint32_t> cm2;
 
void hopscotch_map_sanity_checks() {
  ASSERT_TRUE(cm2.begin() == cm2.end());
  for (size_t i = 0;i < NINS; ++i) {
    cm2[17 * i] = i;
    um2[17 * i] = i;
  }

  for (size_t i = 0;i < NINS; ++i) {
    assert(cm2[17 * i] == i);
    assert(um2[17 * i] == i);
  }
  assert(cm2.size() == NINS);
  assert(um2.size() == NINS);

  for (size_t i = 0;i < NINS; i+=2) {
    cm2.erase(17*i);
    um2.erase(17*i);
  }
  for (size_t i = 0;i < NINS; i+=2) {
    assert(cm2.count(17*i) == i % 2);
    assert(um2.count(17*i) == i % 2);
    if (cm2.count(17*i)) {
      assert(cm2.find(17*i)->second == i);
    }
  }

  assert(cm2.size() == NINS / 2);
  assert(um2.size() == NINS / 2);

  typedef graphlab::hopscotch_map<uint32_t, uint32_t>::value_type vpair;
  {
    size_t cnt = 0;
    foreach(vpair &v, cm2) {
      ASSERT_EQ(v.second, um2[v.first]);
      ++cnt;
    }
    ASSERT_EQ(cnt, NINS / 2);
  }
  {
    size_t cnt = 0;
    foreach(const vpair &v, cm2) {
      ASSERT_EQ(v.second, um2[v.first]);
      ++cnt;
    }
    ASSERT_EQ(cnt, NINS / 2);
  }

  std::stringstream strm;
  graphlab::oarchive oarc(strm);
  oarc << cm2;
  strm.flush();

  cm2.clear();
  ASSERT_EQ(cm2.size(), 0);
  graphlab::iarchive iarc(strm);
  iarc >> cm2;
  ASSERT_EQ(cm2.size(), NINS / 2);

}






void parallel_map_inserter(graphlab::hopscotch_map<uint32_t, uint32_t>* cm,
                        std::vector<uint32_t>* v,
                       size_t start, size_t end) {
  for (size_t i = start; i < end; ++i) {
    cm->put_sync((*v)[i], i);
  }
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
      assert(t == i);
    }
    std::cout << "10M unordered map successful probes in " << ti.current_time() << std::endl;
    um.clear();
  }

  {
    graphlab::cuckoo_map_pow2<uint32_t, uint32_t, 3, uint32_t> cm(-1, 128);
    
    //cm.reserve(102400);
    ti.start();
    for (size_t i = 0;i < NUM_ELS; ++i) {
      cm[v[i]] = i;
  //    if (i % 1000000 == 0) std::cout << cm.load_factor() << std::endl;

    }
    std::cout << NUM_ELS / 1000000 << "M cuckoo map pow2 inserts in " << ti.current_time() << " (Load factor = " << cm.load_factor() << ")" << std::endl;

    graphlab::memory_info::print_usage();

    ti.start();
    for (size_t i = 0;i < 10000000; ++i) {
      size_t t = cm[v[i]];
      assert(t == i);
    }
    std::cout << "10M cuckoo map pow2 successful probes in " << ti.current_time() << std::endl;
  }

{
    graphlab::hopscotch_map<uint32_t, uint32_t> cm;
    
    ti.start();
    for (size_t i = 0;i < NUM_ELS; ++i) {
      cm[v[i]] = i;
//      if (i % 1000000 == 0) std::cout << cm.load_factor() << std::endl;

    }
    std::cout << NUM_ELS / 1000000 << "M hopscotch inserts in " << ti.current_time() << " (Load factor = " << cm.load_factor() << ")" << std::endl;

    graphlab::memory_info::print_usage();

    ti.start();
    for (size_t i = 0;i < 10000000; ++i) {
      size_t t = cm[v[i]];
      assert(t == i);
    }
    std::cout << "10M hopscotch successful probes in " << ti.current_time() << std::endl;

  }


{
    graphlab::hopscotch_map<uint32_t, uint32_t> cm;
    
    graphlab::thread_group thrgroup;
    ti.start();
    for (size_t i = 0;i < 2; ++i) {
      thrgroup.launch(boost::bind(parallel_map_inserter, &cm, &v, 
                                  i * NUM_ELS / 2, (i + 1) * NUM_ELS / 2));
    }
    thrgroup.join();
    std::cout << NUM_ELS / 1000000 << "M hopscotch parallel inserts in " << ti.current_time() << " (Load factor = " << cm.load_factor() << ")" << std::endl;

    graphlab::memory_info::print_usage();

    ti.start();
    for (size_t i = 0;i < 10000000; ++i) {
      std::pair<bool, uint32_t> res = cm.get_sync(v[i]);
      assert(res.first && res.second == i);
    }
    std::cout << "10M hopscotch successful probes in " << ti.current_time() << std::endl;

  }
}



int main(int argc, char** argv) {
  std::cout << "Hopscotch Table Parallel Access Sanity Checks... \n";
  //sanity_checks();

  std::cout << "Hopscotch Map Sequential Access Sanity Checks... \n";
  //hopscotch_map_sanity_checks();

  std::cout << "Map Benchmarks... \n";
  benchmark();
  std::cout << "Done" << std::endl;
}
