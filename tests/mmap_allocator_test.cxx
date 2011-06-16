#include <cxxtest/TestSuite.h>
#include <graphlab/util/mmap_allocator.hpp>

using namespace graphlab;

struct teststruct{
  int a;
  double b;
};

class MmapAllocatorTestSuite: public CxxTest::TestSuite {
  public:
    MmapAllocatorTestSuite() { 
      unlink("mmaptest.bin");
    }
    
    void test_basic_vector() {
      mmap_allocator allocator("mmaptest.bin");
      mmap_allocator_vector rootvec = allocator.get_vector(0, 8);
      rootvec.resize(128);
      // basic set
      for (size_t i = 0;i < 128; ++i) {
        uint64_t t = (uint64_t)i;
        rootvec.set_entry(i, &t);
      }
      
      // push_back
      for (size_t i = 0;i < 1024; ++i) {
        uint64_t t = i + 128;
        rootvec.push_back(&t);
      }
      
      // read back
      for (size_t i = 0;i < 128 + 1024; ++i) {
        uint64_t t;
        rootvec.get_entry(i, &t);
        TS_ASSERT_EQUALS(t, i);
      }
      
      // get_all
      uint64_t* all = new uint64_t[128 + 1024];
      rootvec.get_all(all, 128 + 1024);
      
      // verify all
      for (size_t i = 0;i < 128 + 1024; ++i) {
        TS_ASSERT_EQUALS(all[i], i);
        all[i] = 129;
      }
      
      // set all to 129
      rootvec.set_all(all, 128 + 1024);
      // verify entries
      for (size_t i = 0;i < 128 + 1024; ++i) {
        uint64_t t;
        rootvec.get_entry(i, &t);
        TS_ASSERT_EQUALS(t, (uint64_t)129);
      }
      rootvec.release();
      allocator.close();
    }
    
    void test_basic_persistance() {
      mmap_allocator allocator("mmaptest.bin");
      mmap_allocator_vector rootvec = allocator.get_vector(0, 8);
      for (size_t i = 0;i < 128 + 1024; ++i) {
        uint64_t t;
        rootvec.get_entry(i, &t);
        TS_ASSERT_EQUALS(t, (uint64_t)129);
      }
        // get_all
      uint64_t* all = new uint64_t[128 + 1024];
      rootvec.get_all(all, 128 + 1024);
      
      // verify all
      for (size_t i = 0;i < 128 + 1024; ++i) {
        TS_ASSERT_EQUALS(all[i], (uint64_t)129);
      }
      rootvec.release();
      allocator.close();
    }
    
    void test_vector_creation() {
      mmap_allocator allocator("mmaptest.bin");
      mmap_allocator_offset_t chararray_offset = allocator.create_vector(1, 4);
      mmap_allocator_vector chararray = allocator.get_vector(chararray_offset, 1);
      
      mmap_allocator_offset_t ts_offset = allocator.create_vector(sizeof(teststruct), 4);
      mmap_allocator_vector ts = allocator.get_vector(ts_offset, sizeof(teststruct));
      
      mmap_allocator_vector rootvec = allocator.get_vector(0, 8);
      // store the offsets
      uint64_t coffset = chararray_offset;
      uint64_t toffset = ts_offset;
      rootvec.set_entry(0, &coffset);
      rootvec.set_entry(1, &toffset);
      
      chararray.set_all("hello world", strlen("hello world") + 1);
      for (size_t i = 0;i < 128; ++i) {
        teststruct t;
        t.a = i;
        t.b = i;
        ts.push_back(&t);
      }
      chararray.set_all("hello world hello world", strlen("hello world hello world") + 1);
      chararray.release();
      ts.release();
      rootvec.release();
      allocator.close();
    }
    
    void test_vector_persistance2() {
      mmap_allocator allocator("mmaptest.bin");
      mmap_allocator_vector rootvec = allocator.get_vector(0, 8);
      
      uint64_t coffset;
      uint64_t toffset;
      rootvec.get_entry(0, &coffset);
      rootvec.get_entry(1, &toffset);
      
      mmap_allocator_vector chararray = allocator.get_vector(coffset, 1);
      mmap_allocator_vector ts = allocator.get_vector(toffset, sizeof(teststruct));
      
      for (size_t i = 3;i < 128 + 1024; ++i) {
        uint64_t t;
        rootvec.get_entry(i, &t);
        TS_ASSERT_EQUALS(t, (uint64_t)129);
      }
      char c[128] = {0};
      chararray.get_all(c, strlen("hello world hello world") + 1);
      TS_ASSERT(strcmp(c, "hello world hello world") == 0);
      for (size_t i = 0;i < 128; ++i) {
        teststruct t;
        ts.get_entry(i, &t);
        TS_ASSERT_EQUALS(t.a, (int)i);
        TS_ASSERT_EQUALS(t.b, i);
      }
      chararray.release();
      ts.release();
      rootvec.release();
      allocator.close();
    }
};
  