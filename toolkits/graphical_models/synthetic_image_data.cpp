
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


/**
 * This application creates a synthetic create and process synthetic
 * image data.
 * 
 * \todo Finish documenting
 *
 *  \author Joseph Gonzalez
 */

#include <Magick++.h> 
#undef restrict

#include <graphlab.hpp>



graphlab::vertex_id_type sub2ind(size_t rows, size_t cols,
                                 size_t r, size_t c) {
  return r * cols + c;
}; // end of sub2ind

std::pair<int,int> ind2sub(size_t rows, size_t cols,
                           size_t ind) {
  return std::make_pair(ind / cols, ind % cols);
}; // end of sub2ind


void make_data(const size_t rows, const size_t cols,
               const size_t ncolors, const double error_rate,
               const std::string& vdata_fn,
               const std::string& edata_fn,
               const std::string& orig_img_fn,
               const std::string& noisy_img_fn) {
  using namespace Magick;

  const double center_r = rows / 2.0;
  const double center_c = cols / 2.0;
  const double max_radius = std::min(rows, cols) / 2.0;

  Image orig_img(Magick::Geometry(rows, cols), "white");
  Image noisy_img(Magick::Geometry(rows, cols), "white");
  
  std::ofstream vdata_fout(vdata_fn.c_str());
  std::ofstream edata_fout(edata_fn.c_str());

  for(size_t r = 0; r < rows; ++r) {
    for(size_t c = 0; c < cols; ++c) {
      // determine the true pixel id
      const graphlab::vertex_id_type vid = sub2ind(rows,cols,r,c);
      // Compute the true pixel value
      const double distance = sqrt((r-center_r)*(r-center_r) + 
                                   (c-center_c)*(c-center_c));
      // Compute ring of sunset
      const uint16_t ring_color =  
        std::floor(std::min(1.0, distance/max_radius) * (ncolors - 1) );
      // Compute the true pixel color by masking with the horizon
      const uint16_t true_color = r < rows/2 ? ring_color : 0;
      // compute the predicted color
      const uint16_t obs_color = graphlab::random::rand01() < error_rate?
        graphlab::random::fast_uniform<uint16_t>(0, ncolors) : true_color;

      const double c1p = double(true_color)/(ncolors-1);
      const Color c1(MaxRGB * c1p, MaxRGB * c1p, MaxRGB * c1p, 0);
      orig_img.pixelColor(c,r,c1);

      const double c2p = double(obs_color)/(ncolors-1);
      const Color c2(MaxRGB * c2p, MaxRGB * c2p, MaxRGB * c2p, 0);
      noisy_img.pixelColor(c,r,c2);

      // Save the prior
      vdata_fout << vid << '\t';
      for(size_t pred = 0; pred < ncolors; ++pred) {
        const double prior = obs_color == pred? error_rate : (error_rate) /(ncolors - 1);
        vdata_fout << prior << (pred+1 < ncolors? '\t' : '\n');
      }

      // Add the edges
      if(r + 1 < rows) 
        edata_fout << vid << '\t'
                   << sub2ind(rows,cols,r+1,c) << '\n';
      if(c + 1 < cols) 
        edata_fout << vid << '\t'
                   << sub2ind(rows,cols,r,c+1) << '\n';
    } // end of loop over cols
  } // end of loop over rows

  vdata_fout.close();
  edata_fout.close();
  orig_img.write(orig_img_fn);
  noisy_img.write(noisy_img_fn);

} // end of make data



int main(int argc, char** argv) {
  std::cout << "Create a synthetic noisy image." << std::endl;

  // Set initial values for members ------------------------------------------->
  size_t ncolors = 5;
  double error_rate = 0.5;
  size_t nrows = 200;
  size_t ncols = 200;
 
  std::string vdata_fn = "synth_vdata.tsv";
  std::string edata_fn = "synth_edata.tsv";

  std::string orig_img_fn = "orig_img.jpeg";
  std::string noisy_img_fn = "noisy_img.jpeg";
  std::string pred_fn = "pred_img.jpeg";

 



  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Create synthetic prediction", false);
 
  clopts.attach_option("vdata", &vdata_fn, vdata_fn,
                       "Vertex prior filename");
  clopts.attach_option("edata", &edata_fn, edata_fn,
                       "Adjacency information");

  clopts.attach_option("ncolors",
                       &ncolors, ncolors,
                       "The number of colors in the noisy image");
  clopts.attach_option("error_rate",
                       &error_rate, error_rate,
                       "Standard deviation of noise.");
  clopts.attach_option("nrows",
                       &nrows, nrows,
                       "The number of rows in the noisy image");
  clopts.attach_option("ncols",
                       &ncols, ncols,
                       "The number of columns in the noisy image");
  
  clopts.attach_option("orig",
                       &orig_img_fn, orig_img_fn,
                       "Original image file name.");
  clopts.attach_option("noisy",
                       &noisy_img_fn, noisy_img_fn,
                       "Noisy image file name.");
  clopts.attach_option("pred",
                       &pred_fn, pred_fn,
                       "Predicted image file name.");
    
  ///! Initialize control plain using mpi
  const bool success = clopts.parse(argc, argv);
  if(!success) {
    return EXIT_FAILURE;
  }
  
  make_data(nrows, ncols, ncolors, error_rate,
            vdata_fn, edata_fn,
            orig_img_fn, noisy_img_fn);

  return EXIT_SUCCESS;
} // End of main





// void save_image(const size_t rows, const size_t cols,
//                 const std::vector<pred_pair_type>& values,
//                 const std::string& fname) {
//   using namespace Magick;
//   std::cout << "NPixels: " << values.size() << std::endl;
//   // determine the max and min colors
//   float max_color = -std::numeric_limits<float>::max();
//   float min_color =  std::numeric_limits<float>::max();
//   foreach(pred_pair_type pair, values) {
//     max_color = std::max(max_color, pair.second);
//     min_color = std::min(min_color, pair.second);
//   }
//   Image img(Magick::Geometry(rows, cols), "white");
//   // img.modifyImage();
//   // Pixels img_cache(img);
//   // PixelPackets* pixels = img_cache.
//   foreach(pred_pair_type pair, values) {
//     std::pair<int,int> coords = ind2sub(rows,cols, pair.first);
//     float value = (pair.second - min_color) / (max_color - min_color);
//     Color color(MaxRGB * value, MaxRGB * value, MaxRGB * value, 0);
//     img.pixelColor(coords.second, coords.first, color);
//   }
//   img.write(fname);
// }
