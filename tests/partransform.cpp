#include <vector>
#include <algorithm>
#include <graphlab.hpp>
#include <graphlab/util/par_transform.hpp>

#include <graphlab/macros_def.hpp>

double vec_sum(const std::vector<double>& data) {
  double sum = 0;
  foreach(double d, data) sum += d;
  return sum;
}

void single_op(double& val) {
  val += (cos(val) + 1) * (val + 1.5) * (sqrt(abs(val)) + 1) + sin(val);
}

class iterative_functor {
 public:
  size_t local_ops;
  iterative_functor(size_t local_ops) :
    local_ops(local_ops){};
  void operator()(double& val) {
    for(size_t i = 0; i < local_ops; ++i) {
      single_op(val);
    }
  }
};



int main(int argc, char** argv) {
  graphlab::command_line_options clopts;

  size_t num_elems = 100000;
  size_t local_ops = 10;
  clopts.attach_option("elems", &num_elems, num_elems, "number of elements");
  clopts.attach_option("ops", &local_ops, local_ops, "number of local operations");

  assert(clopts.parse(argc, argv));
  
  graphlab::timer timer;  
  std::vector<double> data(num_elems, 0);
  iterative_functor iterative_fun(local_ops);


  // std::cout << "Timing standard for loop with single op: " << std::endl;
  // timer.start();
  // foreach(double& val, data) {
  //   single_op(val);
  // }
  // std::cout << "Runtime: " << timer.current_time() << std::endl;
  // std::cout << "Sum: " << vec_sum(data) << std::endl;  
  // std::fill(data.begin(), data.end(), 0);
  
  // std::cout << "Timing standard for loop with multiple op: " << std::endl; 
  // timer.start();
  // foreach(double& val, data) {
  //   iterative_fun(val);
  // }
  // std::cout << "Runtime: " << timer.current_time() << std::endl;
  // std::cout << "Sum: " << vec_sum(data) << std::endl;  
  // std::fill(data.begin(), data.end(), 0);



  // std::cout << "Timing parallel for loop with single op: " << std::endl;
  // timer.start();
  // graphlab::par_transform(data.begin(), data.end(), single_op, clopts);
  // std::cout << "Runtime: " << timer.current_time() << std::endl;
  // std::cout << "Sum: " << vec_sum(data) << std::endl;  
  // std::fill(data.begin(), data.end(), 0);


  
  std::cout << "Timing parallel for loop with multiple op:" << std::endl; 
  timer.start();
  graphlab::par_transform(data.begin(), data.end(), iterative_fun, clopts);
  std::cout << "Runtime: " << timer.current_time() << std::endl;
  std::cout << "Sum: " << vec_sum(data) << std::endl;  
  std::fill(data.begin(), data.end(), 0);



  std::cout << "Timing standard for loop with multiple op: " << std::endl; 
  timer.start();
  foreach(double& val, data) {
    iterative_fun(val);
  }
  std::cout << "Runtime: " << timer.current_time() << std::endl;
  std::cout << "Sum: " << vec_sum(data) << std::endl;  
  std::fill(data.begin(), data.end(), 0);


  
}
