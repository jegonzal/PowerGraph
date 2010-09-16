#ifndef IMAGE_DENOISE_HPP
#define IMAGE_DENOISE_HPP




#include <iostream>
#include <fstream>


#include <graphlab.hpp>

#include "data_structures.hpp"





/** Construct denoising grid model based on the image */
void construct_denoise_graph(image& img,
                             size_t num_asgs,
                             double sigma,
                             gl::graph& graph,
                             double edge_weight ,
                             const std::string& corruption ) {
  // Initialize the vertex data to somethine sensible
  vertex_data vdata;
  vdata.belief.resize(num_asgs);
  vdata.belief.uniform(-std::numeric_limits<double>::max());
  vdata.potential.resize(num_asgs);
  vdata.potential.uniform();
  vdata.potential.normalize();        
  // Add all the vertices
  double sigmaSq = sigma*sigma;
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      // initialize the potential and belief
      uint32_t pixel_id = img.vertid(i, j);
      vdata.potential.var() = vdata.belief.var() = pixel_id;
      // Set the node potential
      double obs = img.pixel(i, j);
      if(corruption == "gaussian") {
        for(size_t pred = 0; pred < num_asgs; ++pred) {
          vdata.potential.logP(pred) = 
            -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
        }
      } else if(corruption == "flip") {
        for(size_t pred = 0; pred < num_asgs; ++pred) {
          vdata.potential.logP(pred) = obs == pred? 0 : -sigma;
        }
      } else {
        std::cout << "Invalid corruption!" << std::endl;
        exit(1);
      }
      vdata.potential.normalize();

      // Set the initial assignment
      vdata.asg = vdata.potential.sample();
      
      // Store the actual data in the graph
      size_t vert_id = graph.add_vertex(vdata);

      // Color the graph with a checkerboard pattern
      graph.color(vert_id) = ((i % 2) == 0) ^ ((j % 2) == 0) ? 1 : 0;

      
      // Ensure that we are using a consistent numbering
      assert(vert_id == pixel_id);

      
    } // end of for j in cols
  } // end of for i in rows

  // Construct an edge blob
  edge_data edata;
  edata.factor_id = 0;
  edata.message.resize(num_asgs);
  edata.weight = edge_weight;
  
  // Add all the edges
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) {
        vertex_id_t target = img.vertid(i-1, j);
        edata.message.var() = target;
        graph.add_edge(vertid, target, edata);
      }

      if(i+1 < img.rows()) {
        vertex_id_t target = img.vertid(i+1, j);
        edata.message.var() = target;
        graph.add_edge(vertid, target, edata);
      }

      if(j-1 < img.cols()) {
        vertex_id_t target = img.vertid(i, j-1); 
        edata.message.var() = target;
        graph.add_edge(vertid, target, edata);
      }

      if(j+1 < img.cols()) {
        vertex_id_t target = img.vertid(i, j+1);
        edata.message.var() = target;
        graph.add_edge(vertid, target, edata);
      }
    } // end of for j in cols
  } // end of for i in rows

  // Finalize the graph
  graph.finalize();

} // End of construct graph



/** Construct ising model with no node potentials */
void construct_ising_graph(size_t rows, 
                           size_t cols,
                           size_t num_asgs,
                           gl::graph& graph) {
  // Initialize the vertex data to somethine sensible
  vertex_data vdata;
  vdata.belief.resize(num_asgs);
  vdata.belief.uniform(-std::numeric_limits<double>::max());
  vdata.potential.resize(num_asgs);
  vdata.potential.uniform();
  vdata.potential.normalize();        
  for(size_t i = 0; i < rows; ++i) {
    for(size_t j = 0; j < cols; ++j) {
      // initialize the potential and belief
      uint32_t pixel_id = image::vertid(rows, cols, i, j);
      vdata.potential.var() = vdata.belief.var() = pixel_id;
      // Set the initial assignment
      vdata.asg = vdata.potential.sample();      
      // Store the actual data in the graph
      size_t vert_id = graph.add_vertex(vdata);
      // Ensure that we are using a consistent numbering
      assert(vert_id == pixel_id);
      // Color the graph with a checkerboard pattern
      graph.color(vert_id) = ((i % 2) == 0) ^ ((j % 2) == 0) ? 1 : 0;           
    } // end of for j in cols
  } // end of for i in rows

  // Construct the edges
  edge_data edata;
  edata.factor_id = 0;
  edata.message.resize(num_asgs);
  edata.weight = 0;
  
  // Add all the edges
  for(size_t i = 0; i < rows; ++i) {
    for(size_t j = 0; j < cols; ++j) {
      size_t vertid = image::vertid(rows, cols, i, j);
      if(i-1 < rows) {
        vertex_id_t target = image::vertid(rows, cols, i-1, j);
        edata.message.var() = target;
        graph.add_edge(vertid, target, edata);
      }

      if(i+1 < rows) {
        vertex_id_t target = image::vertid(rows, cols, i+1, j);
        edata.message.var() = target;
        graph.add_edge(vertid, target, edata);
      }

      if(j-1 < cols) {
        vertex_id_t target = image::vertid(rows, cols, i, j-1); 
        edata.message.var() = target;
        graph.add_edge(vertid, target, edata);
      }

      if(j+1 < cols) {
        vertex_id_t target = image::vertid(rows, cols, i, j+1);
        edata.message.var() = target;
        graph.add_edge(vertid, target, edata);
      }
    } // end of for j in cols
  } // end of for i in rows

  // Finalize the graph
  graph.finalize();

} // End of construct graph



/**
 * This struct represents a denoise problem
 */
struct denoise_problem {
  size_t num_asgs;
  size_t rows;
  size_t cols;
  double sigma;
  double lambda;
  std::string smoothing;
  std::string drawing;
  std::string corruption;
  image original;
  image noisy;
  binary_factor edge_factor;
  gl::graph graph;


  denoise_problem() : num_asgs(0), rows(0), cols(0),
                      sigma(0), lambda(0) { }
  
  denoise_problem(size_t num_asgs, size_t rows, size_t cols,
                  double sigma, double lambda, 
                  const std::string& smoothing,
                  const std::string& drawing,
                  const std::string& corruption) :
    num_asgs(num_asgs), rows(rows), cols(cols),
    sigma(sigma), lambda(lambda),
    smoothing(smoothing),
    drawing(drawing),
    corruption(corruption),
    original(rows, cols),
    noisy(rows, cols),
    edge_factor(0, num_asgs, 0, num_asgs) {

    std::cout << "Creating a synethic image." << std::endl;
    if(drawing == "sunset") 
      original.paint_sunset(num_asgs);
    else if(drawing == "checkerboard")
      original.paint_checkerboard(num_asgs);
    else {
      std::cout << "Invalid drawing type!" << std::endl;
      exit(1);
    }
    
    std::cout << "Corrupting Image. " << std::endl;
    noisy = original;
    if(corruption == "gaussian") 
      noisy.gaussian_corrupt(sigma);
    else if(corruption == "flip")
      noisy.flip_corrupt(num_asgs, 0.5);
    else {
      std::cout << "Invalid corruption type!" << std::endl;
      exit(1);
    }

    // dummy variables 0 and 1 and num_asgs by num_asgs
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
    construct_denoise_graph(noisy, num_asgs, sigma, graph, edge_weight, corruption);
  }

  void save_all();

  
  void save(graphlab::oarchive &arc) const {
    arc << num_asgs;
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
  }
  
  void load(graphlab::iarchive &arc) {
    arc >> num_asgs;
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
  out << "Rings:          " << problem.num_asgs << std::endl
      << "Rows:           " << problem.rows << std::endl
      << "Cols:           " << problem.cols << std::endl
      << "Sigma:          " << problem.sigma << std::endl
      << "Lambda:         " << problem.lambda << std::endl
      << "Smoothing:      " << problem.smoothing << std::endl
      << "Edge Factor:  \n" << problem.edge_factor << std::endl;
  return out;
}


void denoise_problem::save_all(){ 
  std::cout << "Saving original image. " << std::endl;
  original.save("original.pgm");    
  
  std::cout << "Saving corrupted image. " << std::endl;
  noisy.save("corrupted.pgm");
  
  std::cout << "Saving corrupted image tsv. " << std::endl;    
  noisy.save_vec("corrupted.tsv");

  std::cout << "Saving Node potentials. " << std::endl;
  {
    std::ofstream fout("potentials.tsv");
    for(vertex_id_t i = 0; i < graph.num_vertices(); ++i) {
      const vertex_data& vdata = graph.vertex_data(i);
      for(size_t asg = 0; asg < vdata.potential.arity(); ++asg) {
        fout << vdata.potential.logP(asg);
        if(asg + 1 < vdata.potential.arity()) fout << '\t';
      }
      fout << '\n';
    }
    fout.close();
  }
  
  
  std::cout << "Saving Description. " << std::endl;
  {
    std::ofstream fout("desc.txt");
    fout << *this << std::endl;
    fout.close();
  }
  
  save("problem.bin");
}








#endif
