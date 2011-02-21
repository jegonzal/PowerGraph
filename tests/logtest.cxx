#include <string>

#include <cxxtest/TestSuite.h>


#include <graphlab/logger/logger.hpp>

void test_basic_assertions() {
  int i = 1;
  int j = 2;
  ASSERT_LT(i, j);
  ASSERT_LE(i, j);
  
  //ASSERT_MSG(i>j, "%d not greater than %d", i, j);
  std::string a = "abc";
  std::string b = "cde";
  ASSERT_EQ(a, a);
  
}

class LogTestSuite: public CxxTest::TestSuite {
 public:
  void test_log() {

    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_file("logtest.logger");
    global_logger().set_log_to_console(false);
    logger(LOG_INFO, "this should only be in the file");

    global_logger().set_log_to_console(true);
    logger(LOG_WARNING, "you should see this both the console and file");
    logstream(LOG_INFO) << "log info again! but with the stream" << std::endl;
    
    global_logger().set_log_file("");
    logger(LOG_ERROR, "this is only in the console");
    logger(LOG_INFO, "console only too");

    logger(LOG_FATAL, "test format strings: %d %s %f", 1, "123", 99.5);
    test_basic_assertions();
 }
};

