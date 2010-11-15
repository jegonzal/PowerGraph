#ifndef IMAGE_HPP
#define IMAGE_HPP

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

#include <boost/random.hpp>


#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>

/** A simple struct represent a gray scale image */
class image {
  size_t _rows, _cols;
  std::vector<double> data;
public:

  /** Create an empty image */
  image() : _rows(0), _cols(0), data(0,0) { }
  
  /** Create an image of a fixed size */
  image(size_t rows, size_t cols) : 
    _rows(rows), _cols(cols), data(rows * cols, 0) { }

  void resize(size_t rows, size_t cols);

  /** Get the number of rows */
  size_t rows() const { return _rows; }

  /** Get the number of columns */
  size_t cols() const { return _cols; }

  /** get the number of pixels */
  size_t pixels() const { return _rows * _cols; }

  /** A function to read a pixel */
  double& pixel(size_t i, size_t j) { return data[vertid(i,j)]; }
  double pixel(size_t i, size_t j) const { return data[vertid(i,j)]; }
  
  /** Linear indexing */
  double& pixel(size_t i) { return data.at(i); }
  double pixel(size_t i) const { return data.at(i); }

  /** Get the vertex id of a pixel */
  size_t vertid(size_t i, size_t j) const { 
    return vertid(_rows, _cols, i, j); 
  }

  /** get the vertex id from the pixel location */
  static size_t vertid(size_t rows, size_t cols, size_t i, size_t j) {
    assert(i < rows);
    assert(j < cols);    
    return i * cols + j; 
  }
  
  
  /** Get the pixel address from the vertex id */
  std::pair<size_t, size_t> loc(size_t vertex) const;

  
  /** A function to save the image to a file in pgm format */
  void save(const char* filename) const;

  void save_vec(const char* filename) const {
    std::ofstream os(filename);
    assert(os.good());
    for(size_t i = 0; i < pixels(); ++i) {
      os << pixel(i) << "\n";
    }
    os.flush();
    os.close();
  }

  
  /** paint a beautiful sunset */
  void paint_sunset(size_t states);
  
  void paint_checkerboard(size_t states, size_t blocks);

  void paint_uniform(double value);

  
  /** Add random noise to the image */
  void gaussian_corrupt(double sigma);

  void flip_corrupt(size_t states, double flip_prob);

  double min() {
    return *std::min_element(data.begin(), data.end());
  }

  double max() {
    return *std::max_element(data.begin(), data.end());
  }

  void save(graphlab::oarchive &oarc) const {
    oarc << _rows;
    oarc << _cols;
    oarc << data;
  }
  
  void load(graphlab::iarchive &iarc) {
    iarc >> _rows;
    iarc >> _cols;
    iarc >> data;
  }


};



/** Generate a normally distributed random number N(mu, sigma^2) */
// std::pair<double, double> randn(double mu = 0, double sigma = 1 ); 


// IMPLEMENTATION =============================================================>


void image::resize(size_t rows, size_t cols) {
  _rows = rows;
  _cols = cols;
  data.resize(rows * cols, 0);
}
  



// static size_t image::vertid(size_t rows, size_t cols, size_t i, size_t j)  {
//   assert(i < rows);
//   assert(j < cols);    
//   return i * cols + j; 
// }


/** Get the vertex id of a pixel */
std::pair<size_t, size_t> image::loc(size_t vertexid) const {
  assert(vertexid < _rows * _cols);
  return std::make_pair( vertexid / _cols, vertexid % _cols);
}


void image::save(const char* filename) const {
  assert(_rows > 0 && _cols > 0);
  std::ofstream os(filename);
  os << "P2" << std::endl
     << _cols << " " << _rows << std::endl
     << 255 << std::endl;
  // Compute min and max pixel intensities
  double min = data[0]; double max = data[0];
  for(size_t i = 0; i < _rows * _cols; ++i) {
    min = std::min(min, data[i]);
    max = std::max(max, data[i]);
  }

  // Save the image (rescaled)
  for(size_t r = 0; r < _rows; ++r) {
    for(size_t c = 0; c < _cols; c++) {
      if(min != max) {
        int color = 
          static_cast<int>(255.0 * (pixel(r,c) - min)/(max-min));
        os << color;
      } else { os << min; }
      if(c != _cols-1) os << "\t";
    }
    os << std::endl;
  } 
  os.flush();
  os.close();
} // end of save



void image::paint_sunset(size_t states) {
  const double center_r = rows() / 2.0;
  const double center_c = cols() / 2.0;
  const double max_radius = std::min(rows(), cols()) / 2.0;
  // Fill out the image
  for(size_t r = 0; r < rows(); ++r) {
    for(size_t c = 0; c < cols(); ++c) {
      double distance = sqrt((r-center_r)*(r-center_r) + 
                             (c-center_c)*(c-center_c));
      // If on top of image
      if(r < rows() / 2) {
        // Compute ring of sunset
        size_t ring = 
          static_cast<size_t>(std::floor(std::min(1.0, distance/max_radius)
                                         * (states - 1) ) );
        pixel(r,c) = ring;
      } else {
        size_t blockx = r / 20;
        size_t blocky = (c + 20 * sin(10.0*r/rows())) / 20;
        size_t index = blockx + 2*blocky;
        pixel(r,c) = index % states;

      }
    }
  }
} // end of paint_beatiful_sunset


void image::paint_checkerboard(size_t states, size_t blocks = 10) {
  size_t block_size = std::min(rows(), cols() ) / blocks;
  // Fill out the image
  for(size_t r = 0; r < rows(); ++r) {
    for(size_t c = 0; c < cols(); ++c) {
      size_t blockx = r / block_size;
      size_t blocky = c / block_size;
      size_t index = blockx + blocky * block_size;
      pixel(r,c) = index % states;
    }
  }
} // end of paint_beatiful_sunset


void image::paint_uniform(double value) {
  for(size_t i = 0; i < rows()*cols(); ++i) {
    pixel(i) = value;    
  }
} // end of paint_uniform





/** corrupt the image with gaussian noise */
void image::gaussian_corrupt(double sigma) {
  //  boost::mt19937 rng;
  boost::lagged_fibonacci607 rng;
  boost::normal_distribution<double> noise_model(0, sigma);
  for(size_t i = 0; i < rows() * cols();  ) {
    // Corrupt two pixels at a time.
    pixel(i++) += noise_model(rng);
  }
} // end of corrupt_image

/** flip_corrupt */
void image::flip_corrupt(size_t states, double prob_flip) {
  boost::mt19937 rng;
  boost::uniform_real<double> dist01;

  for(size_t i = 0; i < rows() * cols();  ++i) {
    double p = dist01(rng);
    if(p < prob_flip) pixel(i) = rand() % states;
  }
} // end of corrupt_image



// /** generate a normally distributed iid pair */
// std::pair<double, double> randn(double mu , double sigma ) {
//   // Generate a N(0,1) from a Unif(0,1) using Box-Muller generator:
//   double u1 = static_cast<double>(rand()) / RAND_MAX;
//   double u2 = static_cast<double>(rand()) / RAND_MAX;
//   double coeff = std::sqrt(-2.0 * std::log(u1));
//   double n1 = coeff * std::cos(2.0 * M_PI * u2) ;
//   double n2 = coeff * std::sin(2.0 * M_PI * u2) ;
//   // Adjust for mean and variance
//   n1 = sigma * n1 + mu;
//   n2 = sigma * n2 + mu; 
//   return std::make_pair(n1, n2);
// } // end of randn

#include <graphlab/macros_undef.hpp>

#endif
