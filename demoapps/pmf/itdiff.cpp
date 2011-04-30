#include <itpp/itbase.h>
#include <cassert>
#include <iostream>
using namespace itpp;
int main(int argc, char** argv) {
  assert(argc == 3);
  
  mat U1, U2, V1, V2;
  it_file input(argv[1]);
  input >> Name("User") >> U1;
  input >> Name("Movie") >> V1;
  it_file input2(argv[2]);
  input2 >> Name("User") >> U2;
  input2 >> Name("Movie") >> V2;
  
  assert(U1.cols() == U2.cols());
  assert(U1.rows() == U2.rows());
  assert(V1.cols() == V2.cols());
  assert(V1.rows() == V2.rows());
  
  assert(sumsum(U1 - U2) < U1.rows() * U1.cols() * 1E-10);
  assert(sumsum(V1 - V2) < V1.rows() * V1.cols() * 1E-10);
}
