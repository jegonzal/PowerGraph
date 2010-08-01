
#include <iostream>
#include <graphlab.hpp>

struct worker : public graphlab::runnable {
  size_t id;
  void run() {
    for(size_t i = 0; i < 5; ++i) {
      std::cout << id << "\t"
                << graphlab::random::rand01()
                << "\t"
                << graphlab::random::rand_int(6)
                << std::endl;
    }
  }
};


int main(int argc, char** argv) {
  graphlab::mutable_queue<size_t, double> queue;

  std::cout << "Testing queue ordering in mutable queue." << std::endl;

  for(size_t i = 0; i < 10; ++i) {
    queue.push(i, 1.0);
  }

  for(size_t i = 0; i < 10; ++i) {
    queue.update(i, 1.0);
  }


  for(size_t i = 0; i < 10; ++i) {
    size_t j = queue.pop().first;
    std::cout << j << std::endl;
  }


  // for(size_t i = 0; i < 10; ++i) {
  //   size_t x = 0;
  //   std::cout << graphlab::random::rand01()
  //             << "\t"
  //             << graphlab::random::rand_int(6)
  //             << std::endl;
  // }

  graphlab::thread_group threads;
  std::vector<worker> workers(3);
  for(size_t i = 0; i < workers.size(); ++i) {
    workers[i].id = i;
    threads.launch(&workers[i]);
  }
  threads.join();

  
  
  
  return 0; 
}
