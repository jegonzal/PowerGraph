/**
 * Run parallel junction tree gibbs sampling on a factorized model
 */

#include <cstdlib>
#include <iostream>


#include <graphlab.hpp>



#include "data_structures.hpp"

#include "image_denoise.hpp"

#include <graphlab/macros_def.hpp>





int main(int argc, char** argv) {
  std::cout << "make the image denoising alchemy problem"
            << std::endl;

  std::string model_filename = "image";
  std::string drawing = "sunset";
  std::string corruption = "gaussian";
  std::string smoothing = "laplace";
  double lambda = 2;
  double sigma = 1;
  size_t rows = 100;
  size_t rings = 5;
  


  

  // Command line parsing
  graphlab::command_line_options clopts("Make the alchemy image");
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


  denoise_problem problem(rings, rows, rows,
                          sigma, lambda,
                          smoothing,
                          drawing,
                          corruption);

  problem.save(model_filename + ".bin");
  problem.save_all(model_filename + ".alchemy");




  return EXIT_SUCCESS;
} // end of main













#include <graphlab/macros_undef.hpp>


