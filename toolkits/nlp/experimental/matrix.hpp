#ifndef MATRIX
#define MATRIX


template<typename T>
class matrix {
private:
  size_t rows, cols;
  std::vector<T> data;

  const size_t linear_index(const size_t& i, const size_t& j) const {
    assert(i < rows && j < cols);
    return i + j * rows;
  }

public:
  matrix(const size_t& rows, const size_t& cols, const T& zero = T(0)) :
    rows(rows), cols(cols), data(rows*cols, zero) { };
  const T& operator()(const size_t& i, const size_t& j) const {
    return data[linear_index(i,j)];
  }
  T& operator()(const size_t& i, const size_t& j) {
    return data[linear_index(i,j)];
  }
  const T& operator()(const size_t& i) const {
    assert(i < data.size());
    return data[i];
  }
  T& operator()(const size_t& i) {
    assert(i < data.size());
    return data[i];
  }
  void zeros() {
    std::fill(data.begin(), data.end(), T(0));
  }

  T sum() const {
    T z(0);
    for(size_t i = 0; i < data.size(); ++i) z += data[i];
    return z;
  }
};


#endif
