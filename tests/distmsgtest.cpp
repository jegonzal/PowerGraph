#include <graphlab/distributed/distributed_control.hpp>
#include <logger/assertions.hpp>
#include <logger/logger.hpp>
#include <graphlab/util/timer.hpp>
using namespace graphlab;

/**********************    Handler Tests    ******************************/
// h0 to h6. A tests all 6 possible handler lengths
// h6_with_data. Tests the complete set of 6 integer arguments + data in the ptr
void h0(distributed_control& dc, size_t source, void* ptr, size_t len) {
  logstream(LOG_INFO) << "h0 received from " << source << "\n";
}

void h1(distributed_control& dc, size_t source, void* ptr, size_t len,
        handlerarg_t a) {
  ASSERT_EQ(a, 1);
  logstream(LOG_INFO) << "h1 received from " << source << "\n";
}

void h2(distributed_control& dc, size_t source, void* ptr, size_t len,
        handlerarg_t a, handlerarg_t b) {
  ASSERT_EQ(a, 1); ASSERT_EQ(b, 2);
  logstream(LOG_INFO) << "h2 received from " << source << "\n";
}

void h3(distributed_control& dc, size_t source, void* ptr, size_t len,
        handlerarg_t a, handlerarg_t b, handlerarg_t c) {
  ASSERT_EQ(a, 1); ASSERT_EQ(b, 2); ASSERT_EQ(c, 3);
  logstream(LOG_INFO) << "h3 received from " << source << "\n";
}

void h4(distributed_control& dc, size_t source, void* ptr, size_t len,
        handlerarg_t a, handlerarg_t b, handlerarg_t c,
        handlerarg_t d) {
  ASSERT_EQ(a, 1); ASSERT_EQ(b, 2); ASSERT_EQ(c, 3);
  ASSERT_EQ(d, 4);
  logstream(LOG_INFO) << "h4 received from " << source  << "\n";
}

void h5(distributed_control& dc, size_t source, void* ptr, size_t len,
        handlerarg_t a, handlerarg_t b, handlerarg_t c,
        handlerarg_t d, handlerarg_t e) {
  ASSERT_EQ(a, 1); ASSERT_EQ(b, 2); ASSERT_EQ(c, 3);
  ASSERT_EQ(d, 4); ASSERT_EQ(e, 5);
  logstream(LOG_INFO) << "h5 received from " << source << "\n";
}

void h6(distributed_control& dc, size_t source, void* ptr, size_t len,
        handlerarg_t a, handlerarg_t b, handlerarg_t c,
        handlerarg_t d, handlerarg_t e, handlerarg_t f) {
  ASSERT_EQ(a, 1); ASSERT_EQ(b, 2); ASSERT_EQ(c, 3);
  ASSERT_EQ(d, 4); ASSERT_EQ(e, 5); ASSERT_EQ(f, 6);
  logstream(LOG_INFO) << "h6 received from " << source  << "\n";
}


void h6_with_data(distributed_control& dc, size_t source, void* ptr, size_t len,
        handlerarg_t a, handlerarg_t b, handlerarg_t c,
        handlerarg_t d, handlerarg_t e, handlerarg_t f) {
  ASSERT_EQ(a, 1); ASSERT_EQ(b, 2); ASSERT_EQ(c, 3);
  ASSERT_EQ(d, 4); ASSERT_EQ(e, 5); ASSERT_EQ(f, 6);
  char* strptr = (char*)ptr;
  ASSERT_STREQ(strptr, "helloworld");
  logstream(LOG_INFO) << "h8_with_data received " << a<< "\n";
}




// Tests the extended handler, passing data in the ptr and 2 floating point args
// of different lengths.
void exth1(distributed_control& dc, size_t source, void* ptr, size_t len,
        float a, double b) {
  char* strptr = (char*)ptr;
  ASSERT_STREQ(strptr, "helloworld");
  ASSERT_EQ(a, 2);
  ASSERT_GT(b, 3.141592); ASSERT_LE(b, 3.141593);
  logstream(LOG_INFO) << "exth1 received from " << source << "\n";
}

// Tests the extended handler, passing data in a struct
struct teststruct {
  double a;
  float b;
  char c[255];
  size_t d;
};
void exth2(distributed_control& dc, size_t source, void* ptr, size_t len,
           teststruct ts) {
  ASSERT_EQ(ts.a, 1); ASSERT_EQ(ts.b, 1); 
  ASSERT_STREQ(ts.c, "helloworld");
  ASSERT_EQ(ts.d, 1);  
  logstream(LOG_INFO) << "exth2 received from " << source << "\n";
}




// Tests the extended handler, passing data in the ptr and 2 floating point args
// of different lengths.
void exts1(distributed_control& dc, size_t source, void* ptr, size_t len,
        float a, double b, std::vector<size_t> &test) {
  char* strptr = (char*)ptr;
  ASSERT_STREQ(strptr, "helloworld");
  ASSERT_EQ(a, 2);
  ASSERT_GT(b, 3.141592); ASSERT_LE(b, 3.141593);
  ASSERT_EQ(test.size() , 10);
  for (size_t i = 0;i < test.size(); ++i) {
    ASSERT_EQ(test[i], i);
  }
  logstream(LOG_INFO) << "exts1 received from " << source << "\n";
}


// test function
void test_handlers(distributed_control &dc) {
  if (dc.procid() == 0) {
    logger(LOG_INFO, "Testing standard remote_call");
    dc.remote_call(1, h0, NULL, 0);
    dc.remote_call(1, h1, NULL, 0, 1);
    dc.remote_call(1, h2, NULL, 0, 1, 2);
    dc.remote_call(1, h3, NULL, 0, 1, 2, 3);
    dc.remote_call(1, h4, NULL, 0, 1, 2, 3, 4);
    dc.remote_call(1, h5, NULL, 0, 1, 2, 3, 4, 5);
    dc.remote_call(1, h6, NULL, 0, 1, 2, 3, 4, 5, 6);

    char helloworld[] = "helloworld";
    dc.remote_call(1, h6_with_data, (void*)("helloworld"), strlen(helloworld)+1, 1, 2, 3, 4, 5, 6);
  
    logger(LOG_INFO, "Testing extended remote_call on standard handlers");
    dc.remote_callx(1, h0, NULL, 0);
    dc.remote_callx(1, h1, NULL, 0, 1);
    dc.remote_callx(1, h2, NULL, 0, 1, 2);
    dc.remote_callx(1, h3, NULL, 0, 1, 2, 3);
    dc.remote_callx(1, h4, NULL, 0, 1, 2, 3, 4);
    dc.remote_callx(1, h5, NULL, 0, 1, 2, 3, 4, 5);
    dc.remote_callx(1, h6, NULL, 0, 1, 2, 3, 4, 5, 6);
    dc.remote_callx(1, h6_with_data, (void*)("helloworld"), strlen(helloworld)+1, 1, 2, 3, 4, 5, 6);

    logger(LOG_INFO, "Testing extended remote_call on extended handlers");
    dc.remote_callx(1, exth1, (void*)("helloworld"), strlen(helloworld)+1, 2, 3.1415926);
    
    teststruct ts = {0};
    ts.a = 1;
    ts.b = 1;
    strcpy(ts.c, helloworld);
    ts.d= 1;
    dc.remote_callx(1, exth2, NULL, 0, ts);

    
    std::vector<size_t> v;
    for (size_t i = 0; i < 10; ++i) {
      v.push_back(i);
    }
    dc.remote_callxs(1, exts1, (void*)("helloworld"), strlen(helloworld)+1, 2, 3.1415926, v);
/*    tmpstruct ts;
    ts.a = 10;
    ts.b = 100.5;
    dc.remote_callx(1, exth2, NULL, 0, ts);*/
  }
  else {
    // instead of sleep 1 here I can implement a more complicated
    // way to wait till done but you know... this is so not necessary
    sleep(1);
  }
  dc.barrier();
}



/************************ Latency test *************************/
double sumlatency;
volatile size_t recvpongs;
timer ti;

void pong_handler(distributed_control& dc, size_t source, void* ptr, size_t len) {
  sumlatency+=ti.current_time();
  recvpongs++;
}

void ping_handler(distributed_control& dc, size_t source, void* ptr, size_t len) {
  dc.remote_call(source, pong_handler, NULL, 0);
}

void test_latency(distributed_control &dc) {
  if (dc.procid() == 0) {
    recvpongs = 0;
    sumlatency = 0;
    for (size_t i = 0;i < 10000; ++i) {
      ti.start();
      dc.remote_call(1, ping_handler, NULL, 0);
      if (i % 100 == 0)  std::cout << ".";
      while(recvpongs != (i+1)) {
        sched_yield();
      }
    }
    std::cout << "\n";
    logstream(LOG_INFO) << "RTT: " << 1000.0 * (sumlatency / recvpongs) << "ms\n";
    dc.print_stats(1);
  }  
  dc.barrier();
}

void test_flood(distributed_control &dc) {
  if (dc.procid() == 0) {
    recvpongs = 0;
    sumlatency = 0;
    ti.start();
    size_t msgsize = 10*1024*1024;
    //size_t msgsize = 1;
    size_t numpings = 20;
    char* c = new char[msgsize];
    for (size_t i = 0;i < numpings ; ++i) {
      dc.remote_call(1, ping_handler, c, msgsize);
    }
    while(recvpongs != (numpings)) {
      sched_yield();
    }
    sumlatency = ti.current_time();
    std::cout << "\n";
    size_t bytestransmitted = msgsize * numpings ;
    double MBps = double(bytestransmitted) / sumlatency / (1024*1024);
    logstream(LOG_INFO) << numpings << "*"<< msgsize << " packet flood in " << sumlatency << " seconds\n";
    logstream(LOG_INFO) << "Flood Rate : " << MBps << " MBps\n";
    dc.print_stats(1);
  }  
  dc.barrier();
}




int main(int argc, char **argv) {
  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(&argc, &argv);
  dc.init_message_processing(2);
  dc.barrier();
  test_handlers(dc);
  
  test_latency(dc);
  test_latency(dc);
  test_latency(dc);
  test_flood(dc);
  //test_latency(dc);
  //test_latency(dc);
  dc.barrier();
}
