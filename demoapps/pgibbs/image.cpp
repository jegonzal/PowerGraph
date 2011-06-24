/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */




#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

#include <boost/random.hpp>


#include <graphlab.hpp>


#include "image.hpp"


#include <graphlab/macros_def.hpp>



void image::save_vec(const char* filename) const {
  std::ofstream os(filename);
  ASSERT_TRUE(os.good());
  for(size_t i = 0; i < pixels(); ++i) {
    os << pixel(i) << "\n";
  }
  os.flush();
  os.close();
}


double image::min() const {
  return *std::min_element(data.begin(), data.end());
}


double image::max() const {
  return *std::max_element(data.begin(), data.end());
}


void image::save(graphlab::oarchive &oarc) const {
  oarc << _rows;
  oarc << _cols;
  oarc << data;
}


void image::load(graphlab::iarchive &iarc) {
  iarc >> _rows;
  iarc >> _cols;
  iarc >> data;
}


/** Generate a normally distributed random number N(mu, sigma^2) */
// std::pair<double, double> randn(double mu = 0, double sigma = 1 ); 


// IMPLEMENTATION =============================================================>


void image::resize(size_t rows, size_t cols) {
  _rows = rows;
  _cols = cols;
  data.resize(rows * cols, 0);
}
  


/** Get the vertex id of a pixel */
size_t image::vertid(size_t i, size_t j) const {
  ASSERT_LT(i, _rows);
  ASSERT_LT(j, _cols);    
  return i * _cols + j; 
}

// static size_t image::vertid(size_t rows, size_t cols, size_t i, size_t j)  {
//   assert(i < rows);
//   assert(j < cols);    
//   return i * cols + j; 
// }


/** Get the vertex id of a pixel */
std::pair<size_t, size_t> image::loc(size_t vertexid) const {
  ASSERT_LT(vertexid, _rows * _cols);
  return std::make_pair( vertexid / _cols, vertexid % _cols);
}


void image::save(const char* filename) const {
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


void image::paint_checkerboard(size_t states, size_t blocks) {
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
