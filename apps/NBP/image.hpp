/*
 * image.hpp
 *
 *  Created on: Oct 10, 2010
 *      Author: lesong
 */

#ifndef IMAGE_HPP_
#define IMAGE_HPP_

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

#include <boost/random.hpp>


#include <graphlab/macros_def.hpp>

/** A simple struct represent a gray scale image */
class image {
  size_t _rows, _cols;
  std::vector<double> data;
public:
  /** Create an image of a fixed size */
  image(size_t rows, size_t cols);

  /** Get the number of rows */
  size_t rows() const;

  /** Get the number of columns */
  size_t cols() const;

  /** get the number of pixels */
  size_t pixels() const;

  /** A function to read a pixel */
  double& pixel(size_t i, size_t j);
  double pixel(size_t i, size_t j) const;

  /** Linear indexing */
  double& pixel(size_t i);
  double pixel(size_t i) const;

  /** Get the vertex id of a pixel */
  size_t vertid(size_t i, size_t j) const;

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
  void paint_sunset(size_t num_rings);

  /** Add random noise to the image */
  void corrupt(double sigma);


};



/** Generate a normally distributed random number N(mu, sigma^2) */
std::pair<double, double> randn(double mu = 0, double sigma = 1 );


// IMPLEMENTATION =============================================================>

/** Create an image of a fixed size */
image::image(size_t rows, size_t cols) :
    _rows(rows), _cols(cols), data(rows * cols) {
  for(size_t i = 0; i < data.size(); ++i)  data[i] = 0;
}

/** The number of rows */
size_t image::rows() const { return _rows; }

/** The number of columns */
size_t image::cols() const { return _cols; }

size_t image::pixels() const { return _rows * _cols; }

/** A function to read a pixel */
double& image::pixel(size_t i, size_t j) {
  return data[vertid(i,j)];
}

double image::pixel(size_t i, size_t j) const {
  return data[vertid(i,j)];
}

/** A function to read a pixel */
double& image::pixel(size_t i) {
  assert(i < data.size());
  return data[i];
}

double image::pixel(size_t i) const {
  assert(i < data.size());
  return data[i];
}


/** Get the vertex id of a pixel */
size_t image::vertid(size_t i, size_t j) const {
  assert(i < _rows);
  assert(j < _cols);
  return i * _cols + j;
}

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



void image::paint_sunset(size_t num_rings) {
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
                                         * (num_rings - 1) ) );
        pixel(r,c) = ring;
      } else {
        pixel(r,c) = 0;
      }
    }
  }
} // end of paint_beatiful_sunset


/** corrupt the image with gaussian noise */
void image::corrupt(double sigma) {
  //  boost::mt19937 rng;
  boost::lagged_fibonacci607 rng;
  boost::normal_distribution<double> noise_model(0, sigma);
  for(size_t i = 0; i < rows() * cols();  ) {
    // Corrupt two pixels at a time.
    pixel(i++) += noise_model(rng);
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

#endif /* IMAGE_HPP_ */

