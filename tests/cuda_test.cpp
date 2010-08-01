
#include <cassert>
#include <iostream>

#include <cuda_runtime.h>

#include <cublas.h>
//#include <cutil.h>

void check(cublasStatus status, const std::string& msg) {
  if (status != CUBLAS_STATUS_SUCCESS) {
    std::cerr << msg << "\n"
              << "CUBLAS status: " << status << std::endl;
    assert(false);
  }
}

int main(int argc, char** argv) {

  size_t n(1000);
  assert(n > 0);

  // Init CUBLAS, and set device vector to {0,...,n-1}.
  check(cublasInit(), "CUBLAS initialization error!");
  float* d_vec = NULL;
  check(cublasAlloc(n, sizeof(float), (void**)&d_vec),
        "CUBLAS device memory allocation error!\n");
  float* h_vec = new float[n];
  for (size_t i(0); i < n; ++i)
    h_vec[i] = i;
  check(cublasSetVector(n, sizeof(float), h_vec, 1, d_vec, 1),
        "CUBLAS failed when copying host vector to device.\n");
  delete [] h_vec;
  h_vec = NULL;

  // Compute sum(device vector).
  float total(cublasSasum(n, d_vec, 1));
  check(cublasGetError(), "CUBLAS error while summing device vector!\n");
  std::cout << "CUBLAS computed sum(0,...," << (n-1) << ") = " << total
            << std::endl;

  check(cublasFree(d_vec),
        "CUBLAS error while freeing vec on device!\n");
  d_vec = NULL;

  return 0;

}
