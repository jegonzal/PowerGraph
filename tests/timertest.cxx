#include <vector>
#include <algorithm>


#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>


using namespace graphlab;


class TimerTestSuite : public CxxTest::TestSuite {
  public:
    void test_timer(void) {
      float t1 = lowres_time_seconds();
      usleep(150000);
      float t2 = lowres_time_seconds();
      usleep(150000);
      float t3 = lowres_time_seconds();
      std::cout << t1 << std::endl;
      std::cout << t2 << std::endl;
      std::cout << t3 << std::endl;
      TS_ASSERT_DELTA(t1 + 0.1, t2, 0.00001);
      TS_ASSERT_DELTA(t2 + 0.1, t3, 0.00001);
    }
};

