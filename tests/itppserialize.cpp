#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include <graphlab.hpp>

namespace graphlab {

template<>
oarchive& operator<< <itpp::Vec<double> > (oarchive& arc, const itpp::Vec<double> &vec) {
  arc << vec.length();
  serialize(arc, vec._data(), sizeof(double)*vec.length());
  return arc;
}
template<>
oarchive& operator<< <itpp::Mat<double> > (oarchive& arc, const itpp::Mat<double> &mat) {
  arc << mat.rows() << mat.cols();
  serialize(arc, mat._data(), sizeof(double)*mat._datasize());  
  return arc;
}

template<>
iarchive& operator>> <itpp::Vec<double> > (iarchive& arc, itpp::Vec<double> &vec) {
  size_t vlength;
  arc >> vlength;
  vec.set_size(vlength);
  deserialize(arc, vec._data(), sizeof(double)*vec.length());
  return arc;
}

template<>
iarchive& operator>> <itpp::Mat<double> > (iarchive& arc, itpp::Mat<double> &mat) {
  size_t rows, cols;
  arc >> rows >> cols;
  mat.set_size(rows,cols);
  deserialize(arc, mat._data(), sizeof(double)*mat._datasize());   
  return arc;
}

}

int main(int argc, char** argv) {
  srand(time(NULL));
  
  // generate a random matrix and vector
  itpp::Vec<double> test;
  itpp::Mat<double> testmat;
  test.set_size(100);
  testmat.set_size(150,150);
  for (size_t i = 0;i < 100; ++i) {
    test(i) = rand();
  }
  for (size_t i = 0;i < 150; ++i) {
    for (size_t j = 0;j < 150; ++j) {
      testmat(i,j) = rand();
    }
  }
 // write to disk
  std::ofstream fout("test.bin");
  graphlab::oarchive oarc(fout);
  oarc << test << testmat;
  fout.close();
  
  // read it back in
  std::ifstream fin("test.bin");
  graphlab::iarchive iarc(fin);
  itpp::Vec<double> test2;
  itpp::Mat<double> testmat2;
  iarc >> test2 >> testmat2;
  fin.close(); 
  
  // check all the values
  for (size_t i = 0;i < 100; ++i) {
    ASSERT_EQ(test(i), test2(i));
  }
  for (size_t i = 0;i < 150; ++i) {
    for (size_t j = 0;j < 150; ++j) {
      ASSERT_EQ(testmat(i,j), testmat2(i,j));
    }
  }
}