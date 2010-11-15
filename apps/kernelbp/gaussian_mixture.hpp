#ifndef GAUSSIAN_MIXTURE_HPP
#define GAUSSIAN_MIXTURE_HPP
#include <graphlab/util/random.hpp>
#include <itpp/itstat.h>
#include <itpp/itbase.h>
using namespace itpp;
template <typename D>
inline D square(const D &v) {
  return v * v;
}

/**
  Describes a Spherical gaussian
*/
template <size_t DIM>
struct SphericalGaussian{
  float center[DIM];
  float weight;
  float sigma;
};



template <>
struct SphericalGaussian<1>{
  float center;
  float weight;
  float sigma;
};

inline float gaussian_log_likelihood(float u, float sigma, float x) {
  return log(1.0/sigma) - square(x - u) / (2 * sigma * sigma);
}

inline float gaussian_likelihood(float u, float sigma, float x) {
  float inner = -square(x - u) / (2 * sigma * sigma);
  if (inner < -20) return 1E-10 / sigma;
  else return exp(-square(x - u) / (2 * sigma * sigma)) / sigma;
}

/**
  Describes a Gaussian Mixture Model over one or more dimensions
*/
template <size_t DIM>
struct GaussianMixture{
  std::vector<SphericalGaussian<DIM> > gaussians;
  GaussianMixture() { }
  GaussianMixture(const mat &center, const mat &sigma, const mat &weights) {

    // get the number of entries.
    gaussians.resize(weights.size());

    // make sure that centers is the right size
    assert((center.rows() == weights.size() && center.cols() == DIM) ||
           (center.rows() == DIM && center.cols() == weights.size()));

    for (size_t i = 0;i < gaussians.size(); ++i) {
      gaussians[i].weight = weights(i);
      gaussians[i].sigma = sigma(i);
      // do I get row or col?
      vec ctr;
      if (center.rows() == weights.size()) {
        ctr = center.get_row(i);
      }
      else if (center.cols() == weights.size()) {
        ctr = center.get_col(i);
      }

      for (size_t j = 0; j < DIM; ++j) {
        gaussians[i].center[j] = ctr[j];
      }
    }
    normalize();
  }

  template <typename ArrayType>
  float likelihood(const ArrayType &d) const {
    float ret = 0;
    for (size_t i = 0;i < gaussians.size(); ++i) {
      float curll = 0;
      for (size_t j = 0;j < DIM; ++j) {
        curll += gaussian_log_likelihood(gaussians[i].center[j],
                                         gaussians[i].sigma, d[j]);
      }
      ret += gaussians[i].weight * exp(curll);
    }
    return ret;
  }

  void normalize() {
    float t = 0 ;
    for (size_t i = 0;i < gaussians.size(); ++i) t += gaussians[i].weight;
    for (size_t i = 0;i < gaussians.size(); ++i) gaussians[i].weight /= t;
  }

  void print() const {
    for (size_t i = 0; i < gaussians.size(); ++i) {
      std::cout << gaussians[i].weight << ": [" ;
      for (size_t j = 0;j < DIM; ++j) {
        std::cout << gaussians[i].center[j];
        if (j < DIM - 1) std::cout << " ";
      }
      std::cout << "] " << gaussians[i].sigma << std::endl;
    }
  }
  
  std::vector<float> sample() const{
    float d = graphlab::random::rand01();
    float s = 0;
    // select a gaussian
    size_t sel = 0;
    for (size_t i = 0;i < gaussians.size(); ++i) {
      s += gaussians[i].weight;
      if (s >= d) {
        sel = i;
        break;
      }
    }
    std::vector<float> ret;
    ret.resize(DIM);
    for (size_t i = 0;i < DIM; ++i) {
      ret[i] = gaussians[sel].center[i] + graphlab::random::gaussian_rand() * gaussians[sel].sigma;
    }
    return ret;
  }
  
  void simplify(float threshold = 1E-8) {
    normalize();
    std::vector<SphericalGaussian<DIM> > newgaussians;
    for (size_t i = 0; i < gaussians.size(); ++i) {
      if (gaussians[i].weight > threshold) {
        newgaussians.push_back(gaussians[i]);
      }
    }
    gaussians = newgaussians;
    normalize();
  }
};



template <>
struct GaussianMixture<1>{
  std::vector<SphericalGaussian<1> > gaussians;
  GaussianMixture() { }
  GaussianMixture(const mat &center, const mat &sigma, const mat &weights) {

    // get the number of entries.
    gaussians.resize(weights.size());

    for (size_t i = 0;i < gaussians.size(); ++i) {
      gaussians[i].weight = weights(i);
      gaussians[i].sigma = sigma(i);
      gaussians[i].center = center(i);
    }
    normalize();
  }


  float likelihood(const float &d) const {
    float ret = 0;
    for (size_t i = 0;i < gaussians.size(); ++i) {
      ret += gaussians[i].weight * gaussian_likelihood(gaussians[i].center,
                                         gaussians[i].sigma, d);
    }
    return ret;
  }

  template <typename ArrayType>
  float likelihood(const std::vector<float> &d) const {
    return likelihood(d[0]);
  }

  float sample() const{
    float d = graphlab::random::rand01();
    float s = 0;
    // select a gaussian
    size_t sel = 0;
    for (size_t i = 0;i < gaussians.size(); ++i) {
      s += gaussians[i].weight;
      if (s >= d) {
        sel = i;
        break;
      }
    }

    return graphlab::random::gaussian_rand() *
           gaussians[sel].sigma + gaussians[sel].center;

  }

  void print() const {
    for (size_t i = 0; i < gaussians.size(); ++i) {
      std::cout << gaussians[i].weight << ": [";
      std::cout << gaussians[i].center;
      std::cout << "] " << gaussians[i].sigma << std::endl;
    }
  }


  void normalize() {
    float t = 0 ;
    for (size_t i = 0;i < gaussians.size(); ++i) t += gaussians[i].weight;
    for (size_t i = 0;i < gaussians.size(); ++i) gaussians[i].weight /= t;
  }

  void simplify(float threshold = 1E-8) {
    normalize();
    std::vector<SphericalGaussian<1> > newgaussians;
    for (size_t i = 0; i < gaussians.size(); ++i) {
      if (gaussians[i].weight > threshold) {
        newgaussians.push_back(gaussians[i]);
      }
    }
    gaussians = newgaussians;
    normalize();
  }
};



/**
  computes the conditional of a gaussian mixture
*/
template <size_t DIM>
struct gmm_conditional {
  GaussianMixture<DIM-1> operator()(const GaussianMixture<DIM> &g,
                                    size_t fixdim,
                                    float val) const{
    GaussianMixture<DIM-1> cond;
    for (size_t i = 0;i < g.gaussians.size(); ++i) {
      SphericalGaussian<DIM-1> newg;
      // set the new centers (its the same, just copy over,
      // ignoring the fixed dimension)
      for (size_t j = 0;j < DIM; ++j) {
        if (j != fixdim) {
          newg.center[j  - (j > fixdim)] = g.gaussians[i].center[j];
        }
      }
      // sigma is the same
      newg.sigma = g.gaussians[i].sigma;
      // weights will be different
      newg.weight = g.gaussians[i].weight *
                    gaussian_likelihood(g.gaussians[i].center[fixdim],
                                        g.gaussians[i].sigma,
                                        val);
      cond.gaussians.push_back(newg);
    }
    cond.normalize();
    return cond;
  }
};


template <>
struct gmm_conditional<2>{
  GaussianMixture<1> operator()(const GaussianMixture<2> &g,
                                    size_t fixdim,
                                    float val) const{
    GaussianMixture<1> cond;
    for (size_t i = 0;i < g.gaussians.size(); ++i) {
      SphericalGaussian<1> newg;
      // set the new centers (its the same, just copy over,
      // ignoring the fixed dimension)
      if (fixdim == 0) newg.center = g.gaussians[i].center[1];
      else if (fixdim == 1) newg.center = g.gaussians[i].center[0];
      // sigma is the same
      newg.sigma = g.gaussians[i].sigma;
      // weights will be different
      newg.weight = g.gaussians[i].weight *
                    gaussian_likelihood(g.gaussians[i].center[fixdim],
                                        g.gaussians[i].sigma,
                                        val);
      cond.gaussians.push_back(newg);
    }
    cond.normalize();
    return cond;
  }
};
#endif
