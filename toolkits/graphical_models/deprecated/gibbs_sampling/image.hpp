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


#ifndef PGIBBS_IMAGE_HPP
#define PGIBBS_IMAGE_HPP

#include <cassert>
#include <iostream>
#include <vector>

#include <graphlab.hpp>


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
  size_t vertid(size_t i, size_t j) const;

  static size_t vertid(size_t rows, size_t cols, size_t i, size_t j) {
    ASSERT_LT(i, rows);
    ASSERT_LT(j, cols);    
    return i * cols + j; 
  }
  
  
  /** Get the pixel address from the vertex id */
  std::pair<size_t, size_t> loc(size_t vertex) const;

  
  /** A function to save the image to a file in pgm format */
  void save(const char* filename) const;

  void save_vec(const char* filename) const;
  
  /** paint a beautiful sunset */
  void paint_sunset(size_t states);
  
  void paint_checkerboard(size_t states, size_t blocks = 10);
  
  /** Add random noise to the image */
  void gaussian_corrupt(double sigma);

  void flip_corrupt(size_t states, double flip_prob);

  double min() const;

  double max() const;

  void save(graphlab::oarchive &oarc) const;
  
  void load(graphlab::iarchive &iarc);

};


#endif

