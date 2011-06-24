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
 * Run parallel junction tree gibbs sampling on a factorized model
 */

#include <cstdlib>
#include <iostream>


#include <graphlab.hpp>

#include "image.hpp"
#include "factorized_model.hpp"

#include <graphlab/macros_def.hpp>



/** Construct denoising ising model based on the image */
void construct_denoise_graph(image& img,
                             size_t num_rings,
                             double sigma,
                             const std::string& corruption,
                             factor_t edge_factor,
                             factorized_model& model) {
 
} // End of construct graph




int main(int argc, char** argv) {
  std::cout << "make the image denoising alchemy problem"
            << std::endl;

  std::string model_filename = "image";
  std::string drawing = "sunset";
  std::string corruption = "gaussian";
  std::string smoothing = "square";
  double lambda = 3;
  double sigma = 1;
  size_t rows = 200;
  size_t rings = 7;
  


  

  // Command line parsing
  graphlab::command_line_options clopts("Make the alchemy image", true);
  clopts.attach_option("model", 
                       &model_filename, model_filename,
                       "Alchemy formatted model file");
  clopts.attach_option("drawing", 
                       &drawing, drawing,
                       "drawing type");
  clopts.attach_option("corruption", 
                       &corruption, corruption,
                       "corruption type");
  clopts.attach_option("smoothing", 
                       &smoothing, smoothing,
                       "smoothing type");
  clopts.attach_option("lambda", 
                       &lambda, lambda,
                       "edge parameter");
  clopts.attach_option("sigma", 
                       &sigma, sigma,
                       "noise parameter");
  clopts.attach_option("rows", 
                       &rows, rows,
                       "number of rows and cols");
  clopts.attach_option("rings", 
                       &rings, rings,
                       "number of rings");

  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }


  
  std::cout << "Creating a synethic image." << std::endl;
  image original(rows, rows);
  if(drawing == "sunset") 
    original.paint_sunset(rings);
  else if(drawing == "checkerboard")
    original.paint_checkerboard(rings);
  else {
    std::cout << "Invalid drawing type!" << std::endl;
    exit(1);
  }
  std::cout << "Saving original image. " << std::endl;
  original.save("original.pgm");    

    
  std::cout << "Corrupting Image. " << std::endl;
  image noisy = original;
  if(corruption == "gaussian") 
    noisy.gaussian_corrupt(sigma);
  else if(corruption == "flip")
    noisy.flip_corrupt(rings, 0.75);
  else if(corruption == "ising") 
    noisy = image(rows, rows);
  else {
    std::cout << "Invalid corruption type!" << std::endl;
    exit(1);
  }
  std::cout << "Saving corrupted image. " << std::endl;
  noisy.save("corrupted.pgm");
  

  // dummy variables 0 and 1 and num_rings by num_rings
  std::cout << "Creating edge factor" << std::endl;
  factor_t edge_factor(domain_t(variable_t(0, rings), variable_t(1, rings)));
  // Set the smoothing type
  if(smoothing == "square") {
    edge_factor.set_as_agreement(lambda);
  } else if (smoothing == "laplace") {
    edge_factor.set_as_laplace(lambda);
  } else  {
    std::cout << "Invalid smoothing stype!" << std::endl;
    assert(false);
  }
  std::cout << edge_factor << std::endl;
  
  std::cout << "Constructing factor graph." << std::endl;
  factorized_model model;
  // Add all the node factors
  double sigmaSq = sigma*sigma;
  for(size_t i = 0; i < noisy.rows(); ++i) {
    for(size_t j = 0; j < noisy.cols(); ++j) {
      // initialize the potential and belief
      uint32_t pixel_id = noisy.vertid(i, j);
      variable_t var(pixel_id, rings);
      factor_t factor(var);
      // Set the node potential
      double obs = noisy.pixel(i, j);
      if(corruption == "gaussian") {
        for(size_t pred = 0; pred < rings; ++pred) {
          factor.logP(pred) = 
            -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
        }
      } else if(corruption == "flip") {
        for(size_t pred = 0; pred < rings; ++pred) {
          factor.logP(pred) = obs == pred? 0 : -sigma;
        }
      } else if(corruption == "ising") {
        // Do nothing since we want a uniform node potential
        factor.uniform();
      } else {
        std::cout << "Invalid corruption!" << std::endl;
        exit(1);
      }
      factor.normalize();
      model.add_factor(factor);
    } // end of for j in cols
  } // end of for i in rows

  // Construct edge_factors  
  for(size_t i = 0; i < noisy.rows(); ++i) {
    for(size_t j = 0; j < noisy.cols(); ++j) {
      size_t source = noisy.vertid(i,j);
      variable_t source_var(source, rings);
      if(i+1 < noisy.rows()) {
        vertex_id_t target = noisy.vertid(i+1, j);
        variable_t target_var(target, rings);
        domain_t dom(source_var, target_var);
        edge_factor.set_args(dom);
        model.add_factor(edge_factor);
      }
      if(j+1 < noisy.cols()) {
        vertex_id_t target = noisy.vertid(i, j+1);
        variable_t target_var(target, rings);
        domain_t dom(source_var, target_var);
        edge_factor.set_args(dom);
        model.add_factor(edge_factor);
      }
    } // end of for j in cols
  } // end of for i in rows

  std::cout << "Saving model in alchemy format" << std::endl;
  model.save_alchemy(model_filename + ".alchemy");


  return EXIT_SUCCESS;
} // end of main













#include <graphlab/macros_undef.hpp>


