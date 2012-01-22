#include <graphlab/util/lock_free_pool.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/assertions.hpp>

using namespace graphlab;
lock_free_pool<size_t> pool;

void exec() {
  size_t *s = NULL;
  for (size_t i = 0;i < 1000000; ++i) {
    while(1) {
      s = pool.alloc();
      if (s == NULL) continue;
      else {
       for (size_t j = 0;j < 10; ++j) (*s)++;
       pool.free(s);
       break;
      }
    } 
  }
}


class LockFreePoolTestSuite: public CxxTest::TestSuite {
 public:  
  void test_lock_free_pool() {
    size_t nthreads = 8;
    
    pool.reset_pool(32);
    thread_group g;
    for (size_t i = 0; i < nthreads; ++i) {
      g.launch(exec);
    }
    while(1) {
      try {
        g.join();
        break;
      }
      catch(const char* c ) {
        std::cout << c << "\n";
      }
    }
    
    std::vector<size_t> alldata = pool.unsafe_get_pool_ref();
    size_t total = 0;
    for (size_t i = 0;i < alldata.size(); ++i) {
      total += alldata[i];
    }
    TS_ASSERT_EQUALS(total, 10000000 * nthreads);
  }
};