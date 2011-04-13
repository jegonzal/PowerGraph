#include <vector>
#include <algorithm>
#include <iostream>



struct discrete_variable {
  typedef uint32_t data_type;
  discrete_variable(data_type i, data_type j) : a(i), b(j) { }
  data_type a, b;
  bool operator==(const discrete_variable& other) const { return a == other.a; }
  bool operator!=(const discrete_variable& other) const { return a != other.a; } 
};

std::ostream& operator<<(std::ostream& out, const discrete_variable& df) {
  return out << df.a;
}


#include <graphlab/logger/assertions.hpp>

//#include <graphlab/parallel/pthread_tools.hpp>

//#include <graphlab/factors/discrete_variable.hpp>

// #include <graphlab/factors/unary_factor.hpp>
// #include <graphlab/factors/binary_factor.hpp>
// #include <graphlab/factors/table_factor.hpp>

//#include <graphlab.hpp>



// using namespace graphlab;
                         

void test_variables() {
  std::cout << "Test Variables" << std::endl;
  discrete_variable v1(1, 3);
  std::cout << v1 << std::endl;

  discrete_variable v2(2, 4);
  std::cout << v2 << std::endl;
  
  discrete_variable v3(3, 2);
  std::cout << v3 << std::endl;

  std::cout << v1 << v2 << v3 << std::endl;
  std::cout << (v1 == v1) << (v1 != v2)  << (v1 != v3) << std::endl;
  
  
  int x = 0;
  x++;


  ASSERT_NE(1, 2);  
  ASSERT_NE(v1, v2);
  ASSERT_NE(v1, v3);


  // assert(v1 == v1);
  // assert(v1 != v2);
  // assert(v1 != v3);


  // ASSERT_NE(v1, v2);
  // ASSERT_NE(v1, v3);

  // assert(v1 == v1);
  // assert(v1 != v2);
  // assert(v1 != v3);



  // ASSERT_EQ(v1, v1);
  // ASSERT_LT(v1, v2);
  // ASSERT_LT(v2, v3);       
}



int main(int argc, char** argv) {
  // sets the logging level of graphlab
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  test_variables(); 
 
  // test_domain();
  // test_assignment();

  return 0;
}

