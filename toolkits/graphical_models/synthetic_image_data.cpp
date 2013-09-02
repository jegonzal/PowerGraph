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
 * This application creates a synthetic data to test and demonstrate
 * the structred prediction applications.  The synthetic task is to
 * remove noise from a synthetic noisy image.
 *
 * In addition this application can be used to take the output of
 * the structured prediction tools and rendering the predicted noise
 * free image.
 *
 *
 * 
 *  \author Joseph Gonzalez
 */

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>


#include <cv.h>
#include <highgui.h>  

#include <graphlab.hpp>

using namespace cv;


/**
 * The pixel struct encodes a pixel location and value.
 *
 * Because each pixel corresponds to vertex in the graph we need a
 * mapping between integers and coordinats.  This is accomplished by
 * using the first two lowest order bytes to encode the column and the
 * two highest order bytes to encode the row.
 */
struct pixel {
  uint16_t row, col;
  double value;
  pixel(uint32_t ind = 0, double value = 0) :
    row(ind >> 16), col( ind & ((1 << 16)-1)), value(value) {  }
}; // end of sub2ind



graphlab::vertex_id_type sub2ind(uint16_t r, uint16_t c) {
  ASSERT_LT(r, ((1 << 16)-1));
  ASSERT_LT(c, ((1 << 16)-1));
  return (r << 16) | c;
}; // end of sub2ind




void make_data(const uint16_t rows, const uint16_t cols,
               const size_t ncolors, const double error_rate,
               const std::string& vdata_fn,
               const std::string& edata_fn,
               const std::string& orig_img_fn,
               const std::string& noisy_img_fn) {
  const double center_r = rows / 2.0;
  const double center_c = cols / 2.0;
  const double max_radius = std::min(rows, cols) / 2.0;

  Mat orig_img(cols, rows, CV_8UC1);
  Mat noisy_img(cols, rows, CV_8UC1);
  std::ofstream vdata_fout(vdata_fn.c_str());
  std::ofstream edata_fout(edata_fn.c_str());

  for(size_t r = 0; r < rows; ++r) {
    for(size_t c = 0; c < cols; ++c) {
      // determine the true pixel id
      const graphlab::vertex_id_type vid = sub2ind(r,c);
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
        graphlab::random::fast_uniform<uint16_t>(0, ncolors-1) : true_color;

      const double c1p = double(true_color)/(ncolors-1);
      unsigned char c1 = (unsigned char)(255 * c1p > 255 ? 255 : 255 * c1p);
      orig_img.at<unsigned char>(r, c) = c1;
      const double c2p = double(obs_color)/(ncolors-1);
      unsigned char c2 = (unsigned char)(255 * c2p > 255 ? 255 : 255 * c2p);
      noisy_img.at<unsigned char>(r, c) = c2;

      // Save the prior
      vdata_fout << vid << '\t';
      for(size_t pred = 0; pred < ncolors; ++pred) {
        const double prior = obs_color == pred? error_rate : (error_rate) /(ncolors - 1);
        vdata_fout << prior << (pred+1 < ncolors? '\t' : '\n');
      }

      // Add the edges
      if(r + 1 < rows) 
        edata_fout << vid << '\t' << sub2ind(r+1,c) << '\n';
      if(c + 1 < cols) 
        edata_fout << vid << '\t' << sub2ind(r,c+1) << '\n';
    } // end of loop over cols
  } // end of loop over rows

  vdata_fout.close();
  edata_fout.close();
  imwrite(orig_img_fn, orig_img);
  imwrite(noisy_img_fn, noisy_img);
} // end of make data




void read_data(const std::string& pred_img_fn) {
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  std::vector<pixel> pixels;

  std::string line;
  size_t line_counter = 0;
  uint16_t nrows = 0, ncols = 0;
  size_t min_pixel(-1);
  size_t max_pixel(0);
  while (std::getline(std::cin, line)) {
    graphlab::vertex_id_type vid(-1);
    std::vector<double> values;
    const bool success = qi::phrase_parse
      (line.begin(), line.end(),       
       //  Begin grammar
       (
        qi::ulong_[phoenix::ref(vid) = qi::_1] >> -qi::char_(",") >>
        (qi::double_[phoenix::push_back(phoenix::ref(values), qi::_1)] % -qi::char_(",") )
        )
       ,
       //  End grammar
       ascii::space); 
    if(!success) {
      logstream(LOG_ERROR) << "Error parsing line: " << line_counter << std::endl
                           << "\t\"" << line << "\"";
    }
    ASSERT_GT(values.size(), 0);
    const size_t pred = 
      std::max_element(values.begin(), values.end()) - values.begin();
    min_pixel = std::min(min_pixel, pred);
    max_pixel = std::max(max_pixel, pred);
    const pixel pix(vid, double(pred) / (values.size() - 1) );
    pixels.push_back(pix);  
    nrows = std::max(nrows, pix.row);
    ncols = std::max(ncols, pix.col);
  }
  nrows++; ncols++;
  std::cout << "nrows: " << nrows << std::endl
            << "ncols: " << ncols << std::endl
            << "minp:  " << min_pixel << std::endl
            << "maxp:  " << max_pixel << std::endl;
  Mat pred_img(ncols, nrows, CV_8UC1);
  for(size_t i = 0; i < pixels.size(); ++i) {
    int s = 255 * pixels[i].value;
    pred_img.at<unsigned char>(pixels[i].row, pixels[i].col) = 
                                     (unsigned char)(s >= 255 ? 255 : s);
  }
  imwrite(pred_img_fn, pred_img);
} // end of make data



int main(int argc, char** argv) {
  std::cout << "Create a synthetic noisy image." << std::endl;

  // Set initial values for members ------------------------------------------->
  size_t ncolors = 5;
  double error_rate = 0.5;
  uint16_t nrows = 200;
  uint16_t ncols = 200;
 
  std::string vdata_fn = "synth_vdata.tsv";
  std::string edata_fn = "synth_edata.tsv";

  std::string orig_img_fn = "orig_img.jpeg";
  std::string noisy_img_fn = "noisy_img.jpeg";
  std::string pred_img_fn;

 



  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Create synthetic prediction", false);
 
  clopts.attach_option("vdata", vdata_fn,
                       "Vertex prior filename");
  clopts.attach_option("edata", edata_fn,
                       "Adjacency information");
  clopts.attach_option("ncolors", ncolors,
                       "The number of colors in the noisy image");
  clopts.attach_option("error_rate", error_rate,
                       "Standard deviation of noise.");
  clopts.attach_option("nrows", nrows,
                       "The number of rows in the noisy image");
  clopts.attach_option("ncols", ncols,
                       "The number of columns in the noisy image");
  clopts.attach_option("orig", orig_img_fn,
                       "Original image file name.");
  clopts.attach_option("noisy", noisy_img_fn,
                       "Noisy image file name.");
  clopts.attach_option("pred", pred_img_fn,
                       "Predicted image file name.");
    
  ///! Initialize control plain using mpi
  const bool success = clopts.parse(argc, argv);
  if(!success) {
    return EXIT_FAILURE;
  }

  if(!pred_img_fn.empty()) {
    std::cout << "Reading in predictions" << std::endl;
    read_data(pred_img_fn);
  } else {
    std::cout << "Generating synthetic data" << std::endl;
    make_data(nrows, ncols, ncolors, error_rate,
              vdata_fn, edata_fn,
              orig_img_fn, noisy_img_fn);
  }
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
