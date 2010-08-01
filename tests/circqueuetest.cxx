#include <cxxtest/TestSuite.h>
#include <graphlab/util/synchronized_circular_queue.hpp>

class CircQueueTestSuite : public CxxTest::TestSuite {
  public:
    void test_circ_queue_basic(void) {
      graphlab::synchronized_circular_queue<size_t> cq(128);
      // push 128 elements
      for (size_t i = 0;i < 128;++i) cq.push(i);
      
      // pop 128 elements
      size_t val = 0;      
      for (size_t i = 0;i < 128;++i) {
          TS_ASSERT(cq.safepop(&val));
          TS_ASSERT_EQUALS(val, i);
      }
      
      // queue should be empty now
      TS_ASSERT(! cq.safepop(&val));
    }
    
    
    void test_circ_queue_extension(void) {
      // test self-resizing
      graphlab::synchronized_circular_queue<size_t> cq(16);
      // push 128 elements
      for (size_t i = 0;i < 128;++i) cq.push(i);
      
      // pop 128 elements
      size_t val = 0;      
      for (size_t i = 0;i < 128;++i) {
          TS_ASSERT(cq.safepop(&val));
          TS_ASSERT_EQUALS(val, i);
      }
      
      // queue should be empty now
      TS_ASSERT(! cq.safepop(&val));
    }
    
    void test_circ_queue_loop(void) {
      // test self-resizing
      graphlab::synchronized_circular_queue<size_t> cq(34);

      for (size_t i = 0;i < 24;++i) cq.push(i);
      
      size_t val = 0;
      for (size_t i = 0;i < 24;++i) {
          TS_ASSERT(cq.safepop(&val));
      }
      
      for (size_t i = 0;i < 24;++i) cq.push(i);
      for (size_t i = 0;i < 24;++i) TS_ASSERT(cq.safepop(&val));

      for (size_t i = 0;i < 87;++i) cq.push(i);
      for (size_t i = 0;i < 87;++i) TS_ASSERT(cq.safepop(&val));

      // queue should be empty now
      TS_ASSERT(! cq.safepop(&val));
    }
};
