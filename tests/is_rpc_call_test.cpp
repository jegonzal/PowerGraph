#include <string>
#include <distributed/is_rpc_call.hpp>
#include <distributed/dc.hpp>
#include <distributed/dc_types.hpp>

using namespace graphlab;

int rpc0(distributed_control &dc, procid_t p) { 
  return 0;
}

void rpc1(distributed_control &dc, procid_t p, size_t i) {
}


void basic0(distributed_control &dc, std::string s) {
}


void basic1(std::string s, procid_t p) {
}


void basic2(std::string s) {
}

void basic3() {
}

void basic4(std::string s, ...) {
}


void print( const boost::mpl::bool_<true> &t) {
  std::cout << "yes\n";
}

void print( const boost::mpl::bool_<false> &t) {
  std::cout << "no\n";
}
int main(int argc, char** argv) {
   
  print(dc_impl::is_rpc_call<typeof(rpc0)>::type());
  print(dc_impl::is_rpc_call<typeof(rpc1)>::type());
  print(dc_impl::is_rpc_call<typeof(basic0)>::type());
  print(dc_impl::is_rpc_call<typeof(basic1)>::type());
  print(dc_impl::is_rpc_call<typeof(basic2)>::type());
  print(dc_impl::is_rpc_call<typeof(basic3)>::type());
  print(dc_impl::is_rpc_call<typeof(basic4)>::type());


}
