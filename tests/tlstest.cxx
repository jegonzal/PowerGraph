#include <vector>
#include <algorithm>
#include <boost/bind.hpp>


#include <cxxtest/TestSuite.h>


#include <graphlab.hpp>


using namespace graphlab;


class worker {
public:
  std::vector<size_t>* bits;
  void run() {
    sleep(1);
    std::cout << "thread ID is " << thread::thread_id() << "\n";
    bits->at(thread::thread_id()) = 1;
    std::cout << graphlab::random::rand01() << std::endl;
    std::cout << graphlab::random::rand01() << std::endl;

  }
};



class TLSTestSuite : public CxxTest::TestSuite {
  public:
    void test_tls(void) {
 
      size_t num_workers(100);
      std::vector<worker> workers(num_workers);
      std::vector<size_t> bits(num_workers, 0 );
      thread_group g;
      for (size_t i = 0; i < num_workers; ++i) {
        std::cout <<      graphlab::random::rand01() << std::endl;
        workers[i].bits = &bits;
        g.launch(boost::bind(&worker::run, &workers[i]));
      }
      g.join();
      size_t total = 0;
      for(size_t i = 0; i < bits.size(); ++i) total += bits[i];
      TS_ASSERT_EQUALS(total, num_workers);
    }
};
