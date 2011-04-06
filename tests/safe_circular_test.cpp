#include <cassert>
#include <iostream>
#include <string>
#include <boost/bind.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/safe_circular_char_buffer.hpp>
#include <graphlab/logger/logger.hpp>


using namespace graphlab;

const size_t BUFFER_SIZE(50000);

struct producer  {
  size_t id;
  safe_circular_char_buffer* buffer;
  size_t bytes_written;
  void run() {
    ASSERT_NE(buffer, NULL);
    bytes_written = 0;
    std::string string(10,'a');
    for(size_t i = 0; i < 200000; ++i) {
      //    std::stringstream strm;
      //    strm.good();
      const char* body = string.c_str();
      size_t size = string.size();
      while( ! (size_t(buffer->write(body, size)) == size) );
      bytes_written += size;
      ASSERT_LT(buffer->size(), BUFFER_SIZE);
    }
  }
};

struct consumer_thread  {
  safe_circular_char_buffer* buffer;
  bool running;
  size_t bytes_read;
  void run() {
    ASSERT_NE(buffer, NULL);
    running = true;
    bytes_read = 0;
    size_t count = 0;
    while( running || !buffer->empty() ) {
      ASSERT_LT(buffer->size(), BUFFER_SIZE);
      char* cstr(NULL);
      std::streamsize ammount_read =
        buffer->blocking_introspective_read(cstr, 100);
      bytes_read += ammount_read;
      std::string str(cstr, ammount_read);
      std::sort(str.begin(), str.end());
      
      buffer->advance_head(ammount_read);
      count++;
      if(count % 1000 == 0) 
        std::cout << "Buffer size: " << buffer->size() << std::endl;
    }
  }
};



int main(int argc, char** argv) {
  safe_circular_char_buffer cbuf(BUFFER_SIZE);

  // Initialize producers
  std::vector<producer> producers(8);
  for(size_t i = 0; i < producers.size(); ++i) {
    producers[i].id = i;
    producers[i].buffer = &cbuf;
  }
  // Launch all the producers
  thread_group threads;
  for(size_t i = 0; i < producers.size(); ++i) {
    threads.launch(boost::bind(&producer::run, &producers[i]));
  }
  
  // Launch the consumer;
  consumer_thread consumer;
  consumer.buffer = &cbuf;
  consumer.running = true;
  thread consumerthr = launch_in_new_thread(boost::bind(&consumer_thread::run, &consumer));
  
  threads.join();
  cbuf.stop_reader();
  consumer.running = false;
  consumerthr.join();

  std::cout << "Finished. Comparing:" << std::endl;

  size_t bytes_written(0);
  for(size_t i = 0; i < producers.size(); ++i) 
    bytes_written += producers[i].bytes_written;
  
  std::cout << "Bytes Written: " << bytes_written << std::endl;
  std::cout << "Bytes Read:    " << consumer.bytes_read << std::endl;

  ASSERT_EQ(bytes_written, consumer.bytes_read);
  
  // // one more time with an introspective write, expand, read cycle
  // size_t wrotebytes = 0;
  // size_t readbytes = 0;
  // for (size_t i = 0;i < 1024; ++i) {
  //   char v[128];
  //   char* tmp;
  //   size_t s = cbuf.write(v, 100);
  //   wrotebytes += 100;
  //   while (1) {
  //     size_t r = cbuf.introspective_read(tmp, 30);
  //     cbuf.advance_head(r);
  //     readbytes += r;
  //     if (r == 0) break;
  //   }
  // }
  // ASSERT_EQ(readbytes, wrotebytes);
 
}

