#ifndef IMAGE_DENOISE_HPP
#define IMAGE_DENOISE_HPP




#include <iostream>
#include <fstream>


#include <graphlab.hpp>



#include "image.hpp"
#include "data_structures.hpp"




#include <graphlab/macros_def.hpp>


/** Construct denoising ising model based on the image */
void construct_denoise_graph(image& img,
                             size_t num_rings,
                             double sigma,
                             const std::string& corruption,
                             factor_t edge_factor,
                             factorized_model& model,
                             mrf::graph_type& graph) {
  // Add all the node factors
  double sigmaSq = sigma*sigma;
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      // initialize the potential and belief
      uint32_t pixel_id = img.vertid(i, j);
      variable_t var(pixel_id, num_rings);
      factor_t factor(var);
      // Set the node potential
      double obs = img.pixel(i, j);
      if(corruption == "gaussian") {
        for(size_t pred = 0; pred < num_rings; ++pred) {
          factor.logP(pred) = 
            -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
        }
      } else if(corruption == "flip") {
        for(size_t pred = 0; pred < num_rings; ++pred) {
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
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      size_t source = img.vertid(i,j);
      variable_t source_var(source, num_rings);
      if(i+1 < img.rows()) {
        vertex_id_t target = img.vertid(i+1, j);
        variable_t target_var(target, num_rings);
        domain_t dom(source_var, target_var);
        edge_factor.set_args(dom);
        model.add_factor(edge_factor);
      }
      if(j+1 < img.cols()) {
        vertex_id_t target = img.vertid(i, j+1);
        variable_t target_var(target, num_rings);
        domain_t dom(source_var, target_var);
        edge_factor.set_args(dom);
        model.add_factor(edge_factor);
      }
    } // end of for j in cols
  } // end of for i in rows

  // Construct the clique graph
  std::cout << "Constructing the clique graph" << std::endl;
  construct_mrf(model, graph);

  // Finalize the graph
  std::cout << "Finalizing graph" << std::endl;
  graph.finalize();

} // End of construct graph


/**
 * This struct represents a denoise problem
 */
struct denoise_problem {
  size_t num_rings;
  size_t rows;
  size_t cols;
  double sigma;
  double lambda;
  std::string smoothing;
  std::string drawing;
  std::string corruption;
  image original;
  image noisy;
  factor_t edge_factor;
  factorized_model model;
  mrf::graph_type graph;


  denoise_problem() : num_rings(0), rows(0), cols(0),
                      sigma(0), lambda(0) { }
  
  denoise_problem(size_t num_rings, size_t rows, size_t cols,
                  double sigma, double lambda, 
                  const std::string& smoothing,
                  const std::string& drawing,
                  const std::string& corruption) :
    num_rings(num_rings), rows(rows), cols(cols),
    sigma(sigma), lambda(lambda),
    smoothing(smoothing),
    drawing(drawing),
    corruption(corruption),
    original(rows, cols),
    noisy(rows, cols),
    edge_factor(domain_t(variable_t(0, num_rings), variable_t(1, num_rings))) {

    std::cout << "Creating a synethic image." << std::endl;
    if(drawing == "sunset") 
      original.paint_sunset(num_rings);
    else if(drawing == "checkerboard")
      original.paint_checkerboard(num_rings);
    else {
      std::cout << "Invalid drawing type!" << std::endl;
      exit(1);
    }
    
    std::cout << "Corrupting Image. " << std::endl;
    noisy = original;
    if(corruption == "gaussian") 
      noisy.gaussian_corrupt(sigma);
    else if(corruption == "flip")
      noisy.flip_corrupt(num_rings, 0.5);
    else if(corruption == "ising") 
      noisy = image(rows, cols);
    else {
      std::cout << "Invalid corruption type!" << std::endl;
      exit(1);
    }

    // dummy variables 0 and 1 and num_rings by num_rings
    std::cout << "Creating edge factor" << std::endl;
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
    double edge_weight = 1; // edge_factor.mk_derivative() ;
    std::cout << "Edge Factor weight: " << edge_weight << std::endl;

    std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
    construct_denoise_graph(noisy,
                            num_rings,
                            sigma,
                            corruption,
                            edge_factor,
                            model,
                            graph);
  }


  void save_all(const std::string& alchemy_fn);




  void save(graphlab::oarchive &arc) const {
    arc << num_rings;
    arc << rows;
    arc << cols;
    arc << sigma;
    arc << lambda;
    arc << smoothing;
    arc << drawing;
    arc << corruption;
    arc << original;
    arc << noisy;
    arc << edge_factor;
    arc << graph;
    arc << model;
  }
  
  void load(graphlab::iarchive &arc) {
    arc >> num_rings;
    arc >> rows;
    arc >> cols;
    arc >> sigma;
    arc >> lambda;
    arc >> smoothing;
    arc >> drawing;
    arc >> corruption;
    arc >> original;
    arc >> noisy;
    arc >> edge_factor;
    arc >> graph;
    arc >> model;
  }


  void load(const std::string& filename) {
    std::ifstream fin(filename.c_str());
    graphlab::iarchive iarc(fin);
    iarc >> *this;
    fin.close();
  } // end of load
  
  
  /**
   * \brief save the graph to the file given by the filename
   * 
   */    
  void save(const std::string& filename) const {
    std::ofstream fout(filename.c_str());
    graphlab::oarchive oarc(fout);
    oarc << *this;
    fout.close();
  } // end of save
  
};

std::ostream&
operator<<(std::ostream& out, const denoise_problem& problem) {
  out << "Rings:          " << problem.num_rings << std::endl
      << "Rows:           " << problem.rows << std::endl
      << "Cols:           " << problem.cols << std::endl
      << "Sigma:          " << problem.sigma << std::endl
      << "Lambda:         " << problem.lambda << std::endl
      << "Smoothing:      " << problem.smoothing << std::endl
      << "Edge Factor:  \n" << problem.edge_factor << std::endl;
  return out;
}


void denoise_problem::save_all(const std::string& alchemy_fn) {
    std::cout << "Saving original image. " << std::endl;
    original.save("original.pgm");    
  
    std::cout << "Saving corrupted image. " << std::endl;
    noisy.save("corrupted.pgm");
  
    std::cout << "Saving corrupted image tsv. " << std::endl;    
    noisy.save_vec("corrupted.tsv");

    // std::cout << "Saving Factors. " << std::endl;
    // {
    //   std::ofstream fout("factors.tsv");
    //   foreach(const factor_t& factor, model.factors()) {
    //     fout << factor << '\n';
    //   }
    //   fout.close();
    // }

    std::cout << "Saving Alchemy. " << std::endl;
    {
      model.save_alchemy(alchemy_fn);
    }
    std::cout << "Saving Description. " << std::endl;
    {
      std::ofstream fout("desc.txt");
      fout << *this << std::endl;
      fout.close();
    }
  
    save("problem.bin");
  } // end of save_all






#include <graphlab/macros_undef.hpp>
#endif
