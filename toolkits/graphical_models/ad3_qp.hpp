/*  
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
 *
 * \brief This application performs MAP inference on Markov Nets 
 * provided in standard UAI file format via Dual-Decomposition. 
 *
 *
 *  \authors Dhruv Batra, André Martins, Aroma Mahendru
 */



#ifndef _AD3_QP_HPP_
#define _AD3_QP_HPP_


#include <Eigen/Eigenvalues>
#include <math.h>
#include <limits>
#include "dd_grlab.hpp"


#define NEARLY_ZERO_TOL(a,tol) (((a)<=(tol)) && ((a)>=(-(tol))))
#define NEARLY_EQ_TOL(a,b,tol) (((a)-(b))*((a)-(b))<=(tol))
#define num_max_iterations_QP_ 10
#define EXP 10



////////////////////////////////////////////////////////////////////////////////
// This class implements the Alternating Directions Dual Decompostion as 
// described in:
//
// André F. T. Martins, Mário A. T. Figueiredo, Pedro M. Q. Aguiar,
// Noah A. Smith, and Eric P. Xing.
// "Alternating Directions Dual Decomposition"
// Arxiv preprint arXiv:1212.6550, 2012.
/////////////////////////////////////////////////////////////////////////////////

struct ad3_vertex_program:public admm_vertex_program {

/**
 * \brief Maximize returns the maximum value and configuration with reference to 
 * input additional and variable log potentials. addtional log potential corresponds 
 * factor potentials and variable potential corresponds to sum of lagrange 
 * multipliers and unary potentials divided by degree of the unary vertex.
 */


void Maximize(vertex_type& vertex, vec additional_log_potentials, vec variable_log_potentials,
                Configuration &configuration,
                double *value) {
          
    vector <Configuration> states(vertex.data().nvars,-1);
    *value = -1e12;
    for (int index = 0;
         index < additional_log_potentials.size();
         ++index) {
      double score = additional_log_potentials[index];
      get_configuration_states(vertex,index, &states);
      int offset = 0;
      for (int i = 0; i < vertex.data().nvars; ++i) {
        score += variable_log_potentials[offset+states[i]];
        offset = vertex.data().cards[i];
        
      }
      
      if (configuration < 0 || score > *value) {
        configuration = index;
        *value = score;
      }
    }
    assert(configuration >= 0);
    
  }

 
 void DeleteConfiguration(Configuration &configuration) {
    configuration = -1;
  }

/**
 * \brief InvertAfterInsertion function is used to invert the matrix. Used in solveQP
 */

 bool InvertAfterInsertion(vertex_type& vertex, vector <double> & inverse_A_,
        const vector<Configuration> &active_set, const Configuration &inserted_element) {

  vector<double> inverse_A = inverse_A_;
  int size_A = active_set.size() + 1;
  vector<double> r(size_A);

  r[0] = 1.0;
  for (int i = 0; i < active_set.size(); ++i) {
    // Count how many variable values the new assignment
    // have in common with the i-th assignment.
    int num_common_values = CountCommonValues(vertex, active_set[i], inserted_element);
    r[i+1] = static_cast<double>(num_common_values);
  }

  double r0 = static_cast<double>(CountCommonValues(vertex,
      inserted_element, inserted_element));
  double s = r0;
  for (int i = 0; i < size_A; ++i) {
    if (r[i] == 0.0) continue;
    s -= r[i] * r[i] * inverse_A[i * size_A + i];
    for (int j = i+1; j < size_A; ++j) {
      if (r[j] == 0.0) continue;
      s -= 2 * r[i] * r[j] * inverse_A[i * size_A + j];
    }
  }

    if (NEARLY_ZERO_TOL(s, 1e-9)) {
         if (opts.verbose> 2) {
      cout << "Warning: updated matrix will become singular after insertion."
           << endl;
    }
    return false;
  }

  double invs = 1.0 / s;
  vector<double> d(size_A, 0.0);
  for (int i = 0; i < size_A; ++i) {
    if (r[i] == 0.0) continue;
    for (int j = 0; j < size_A; ++j) {
      d[j] += inverse_A[i * size_A + j] * r[i];
    }
  }

  int size_A_after = size_A + 1;
  inverse_A_.resize(size_A_after * size_A_after);
  for (int i = 0; i < size_A; ++i) {
    for (int j = 0; j < size_A; ++j) {
      inverse_A_[i * size_A_after + j] = inverse_A[i * size_A + j] +
          invs * d[i] * d[j];
    }
    inverse_A_[i * size_A_after + size_A] = -invs * d[i];
    inverse_A_[size_A * size_A_after + i] = -invs * d[i];
  }
  inverse_A_[size_A * size_A_after + size_A] = invs;

  return true;
}

/**
 * \brief InvertAfterRemoval function is used to invert the matrix. Used in solveQP
 */
void InvertAfterRemoval(vector <double> &inverse_A_,const vector<Configuration> &active_set,
                                       int removed_index) {
  vector<double> inverse_A = inverse_A_;
  int size_A = active_set.size() + 1;
  vector<double> r(size_A);

  ++removed_index; // Index in A has an offset of 1.
  double invs = inverse_A[removed_index * size_A + removed_index];
  assert(!NEARLY_ZERO_TOL(invs, 1e-12));
  double s = 1.0 / invs;
  vector<double> d(size_A - 1, 0.0);
  int k = 0;
  for (int i = 0; i < size_A; ++i) {
    if (i == removed_index) continue;
    d[k] = -s * inverse_A[removed_index * size_A + i];
    ++k;
  }

  int size_A_after = size_A - 1;
  inverse_A_.resize(size_A_after * size_A_after);
  k = 0;
  for (int i = 0; i < size_A; ++i) {
    if (i == removed_index) continue;
    int l = 0;
    for (int j = 0; j < size_A; ++j) {
      if (j == removed_index) continue;
      inverse_A_[k * size_A_after + l] = inverse_A[i * size_A + j] -
          invs * d[k] * d[l];
      ++l;
    }
    ++k;
  }
}

/**
 * \brief ComputeActiveSetSimilarities computes Mnz'*Mnz. Used in solveQP
 */
void ComputeActiveSetSimilarities(vertex_type& vertex,
    const vector<Configuration> &active_set,
    vector<double> *similarities) {
  int size = active_set.size();

  // Compute similarity matrix.
  similarities->resize(size * size);
  (*similarities)[0] = 0.0;
  for (int i = 0; i < active_set.size(); ++i) {
    (*similarities)[i*size + i] = static_cast<double>(
        CountCommonValues(vertex,active_set[i], active_set[i]) );
    for (int j = i+1; j < active_set.size(); ++j) {
      // Count how many variable values the i-th and j-th 
      // assignments have in common.
      int num_common_values = CountCommonValues(vertex,active_set[i], active_set[j]);
      (*similarities)[i*size + j] = num_common_values;
      (*similarities)[j*size + i] = num_common_values;
    }
  }
}

/**
 * \brief  ComputeMarginalsFromSparseDistribution computes marginalvalues for unary 
 * factor from given factor distribution.
 */
 
void ComputeMarginalsFromSparseDistribution( vertex_type& vertex, 
    const vector<Configuration> &active_set,
    const vector<double> &distribution,
    vec  &variable_posteriors,
    vec &additional_posteriors) {
    variable_posteriors.setZero();           
    additional_posteriors.setZero();  
    for (int i = 0; i < active_set.size(); ++i) {
    UpdateMarginalsFromConfiguration(vertex,active_set[i],
                                       distribution[i],
                                       variable_posteriors,
                                       additional_posteriors);
    }
  }
  
  
   // Given a configuration with a probability (weight), 
  // increment the vectors of variable and additional posteriors.
  void UpdateMarginalsFromConfiguration(vertex_type& vertex,
    const Configuration &configuration,
    double weight,
    vec &variable_posteriors,
    vec &additional_posteriors) {
    
     vector <Configuration> states(vertex.data().nvars, -1);
     get_configuration_states(vertex, configuration, &states);
     
            int offset = 0;
            
            for (int k = 0; k < vertex.data().nvars; ++k) 
            {   variable_posteriors[offset + states[k]] += weight;
                offset += vertex.data().cards[k];
            }
    additional_posteriors[configuration] += weight;
 
  }
  // Count how many common values two configurations have.
  int CountCommonValues(vertex_type& vertex,Configuration configuration1,
                        Configuration configuration2) {
    
    //assert(states1->size() == states2->size());
    int count = 0;
    vector <Configuration> states1(vertex.data().nvars, -1); 
    vector <Configuration> states2(vertex.data().nvars, -1);
    get_configuration_states(vertex, configuration1, &states1);
    get_configuration_states(vertex, configuration2, &states2);
    for(int i = 0; i< vertex.data().nvars; i++)
    {  if (states1[i] == states2[i])
      { count++;} }
    return count;
  }
  

/**
 * \brief Evaluate returns the maximum value  with reference to 
 * input additional and variable log potentials and configuration. addtional 
 * log potential corresponds factor potentials and variable potential corresponds 
 * to sum of lagrange  * multipliers and unary potentials divided by degree of 
 * the unary vertex.
 */
  
  
void Evaluate(vertex_type& vertex, vec additional_log_potentials, vec variable_log_potentials,
                const Configuration configuration,
                double *value) {
          
    vector<Configuration> states(vertex.data().nvars, -1);
    get_configuration_states(vertex, configuration, &states);
    *value = 0.0;
    int offset = 0;
    for (int i = 0;i<vertex.data().nvars; ++i) {
      *value += variable_log_potentials[offset + states[i]];
      offset = vertex.data().cards[i]; 
    }
    *value += additional_log_potentials[configuration];
  }
  
  
  
  void EigenDecompose(vector<double> *similarities,
                            vector<double> *eigenvalues) {

  int size = sqrt(similarities->size());

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  Eigen::MatrixXd sim(size, size);
  int t = 0;
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      sim(i, j) = (*similarities)[t];
      ++t;
    }
  }
  es.compute(sim);
  const Eigen::VectorXd &eigvals = es.eigenvalues(); 
  eigenvalues->resize(size);
  for (int i = 0; i < size; ++i) {
    (*eigenvalues)[i] = eigvals[i];
  }
  const Eigen::MatrixXd &eigvectors = es.eigenvectors().transpose();
  t = 0;
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      (*similarities)[t] = eigvectors(i, j);
      ++t;
    }
  }

}


// Function to solve each quadratic sub problem for dense factors. 
// It uses active set method. Caching is deactivated
// TODO: Activate caching feature

void SolveQP_dense(vertex_type& vertex,const gather_type& total,
            vec& variable_posteriors, vec& additional_posteriors) {
   vertex_data& vdata = vertex.data();                        

   vec additional_log_potentials = vdata.potentials;
   vec variable_log_potentials = total.neighbor_distribution + total.messages;      
   vector <Configuration> active_set_;
   vector<double> distribution_;
   vector<double> inverse_A_;
  // Initialize the active set.
  
   if (active_set_.size() == 0) {
    variable_posteriors.resize(variable_log_potentials.size());     
    additional_posteriors.resize(additional_log_potentials.size()); 
    distribution_.clear();
    // Initialize by solving the LP, discarding the quadratic
    // term.
    Configuration configuration = -1;
    double value;
    Maximize(vertex, additional_log_potentials, variable_log_potentials,
             configuration,
             &value);
    active_set_.push_back(configuration);
    distribution_.push_back(1.0);

    // Initialize inv(A) as [-M,1;1,0].
    inverse_A_.resize(4);
    inverse_A_[0] = static_cast<double>(
        -CountCommonValues(vertex,configuration, configuration));
    inverse_A_[1] = 1;
    inverse_A_[2] = 1;
    inverse_A_[3] = 0;
  }

  bool changed_active_set = true;
  vector<double> z;
  int num_max_iterations = num_max_iterations_QP_;
  double tau = 0;
  for (int iter = 0; iter < num_max_iterations; ++iter) {
    bool same_as_before = true;
    bool unbounded = false;
    if (changed_active_set) {
      // Recompute vector b.
      vector<double> b(active_set_.size() + 1, 0.0);
      b[0] = 1.0;
      for (int i = 0; i < active_set_.size(); ++i) {
        const Configuration &configuration = active_set_[i];
        double score;
        Evaluate(vertex, additional_log_potentials, variable_log_potentials,
                 configuration,
                 &score);
        b[i+1] = score;
      }
      // Solve the system Az = b.
      z.resize(active_set_.size());
      int size_A = active_set_.size() + 1;
      for (int i = 0; i < active_set_.size(); ++i) {
        z[i] = 0.0;
        for (int j = 0; j < size_A; ++j) {
          z[i] += inverse_A_[(i+1) * size_A + j] * b[j];
        }
      }
      tau = 0.0;
      for (int j = 0; j < size_A; ++j) {
        tau += inverse_A_[j] * b[j];
      }

      same_as_before = false;
    }

    if (same_as_before) {
      // Compute the variable marginals from the full distribution
      // stored in z.
      ComputeMarginalsFromSparseDistribution(vertex, active_set_,
                                             z,
                                             variable_posteriors,
                                             additional_posteriors);
      // Get the most violated constraint
      // (by calling the black box that computes the MAP).
      vec scores = variable_log_potentials;               
      for (int i = 0; i < scores.size(); ++i) {
        scores[i] -= variable_posteriors[i];
      }
      Configuration configuration = -1;
      double value = 0.0;
      
      Maximize(vertex,
                additional_log_potentials, scores,
               configuration,
               &value);
      double very_small_threshold = 1e-9;
      if (value <= tau + very_small_threshold) { // value <= tau.
        // We have found the solution;
        // the distribution, active set, and inv(A) are cached for the next round.
        DeleteConfiguration(configuration);
        return;
      } else {
        for (int k = 0; k < active_set_.size(); ++k) {
          // This is expensive and should just be a sanity check.
          // However, in practice, numerical issues force an already existing
          // configuration to try to be added. Therefore, we always check
          // if a configuration already exists before inserting it.
          // If it does, that means the active set method converged to a
          // solution (but numerical issues had prevented us to see it.)
          if (active_set_[k] == configuration) {                         
            if (opts.verbose > 2) {
              cout << "Warning: value - tau = "
                   << value - tau << " " << value << " " << tau
                   << endl;
            }
            // We have found the solution;
            // the distribution, active set, and inv(A)
            // are cached for the next round.
            DeleteConfiguration(configuration);

            // Just in case, clean the cache.
            // This may prevent eventual numerical problems in the future.
            for (int j = 0; j < active_set_.size(); ++j) {
              if (j == k) continue; // This configuration was deleted already.
              DeleteConfiguration(active_set_[j]);
            }
            active_set_.clear();
            inverse_A_.clear();
            distribution_.clear();

            // Return.
            return;
          }
        }
        z.push_back(0.0);
        distribution_ = z;

        // Update inv(A).
        bool singular = !InvertAfterInsertion(vertex, inverse_A_, active_set_, configuration);
        if (singular) {
          // If adding a new configuration causes the matrix to be singular,
          // don't just add it. Instead, look for a configuration in the null
          // space and remove it before inserting the new one.
          // Right now, if more than one such configuration exists, we just
          // remove the first one we find. There's a chance this could cause
          // some cyclic behaviour. If that is the case, we should randomize
          // this choice.
          // Note: This step is expensive and requires an eigendecomposition.
          // TODO: I think there is a graph interpretation for this problem.
          // Maybe some specialized graph algorithm is cheaper than doing
          // the eigendecomposition.
          vector<double> similarities(active_set_.size() * active_set_.size());
          ComputeActiveSetSimilarities(vertex, active_set_, &similarities);
          
          //cout<<"compute active similarities in solveQP .."<<endl;
          vector<double> padded_similarities((active_set_.size()+2) * 
                                             (active_set_.size()+2), 1.0);
          for (int i = 0; i < active_set_.size(); ++i) {
            for (int j = 0; j < active_set_.size(); ++j) {
              padded_similarities[(i+1)*(active_set_.size()+2) + (j+1)] =
                  similarities[i*active_set_.size() + j];
            }
          }
          padded_similarities[0] = 0.0;
          for (int i = 0; i < active_set_.size(); ++i) {
            double value = static_cast<double>(
                CountCommonValues(vertex, configuration, active_set_[i]));
            padded_similarities[(i+1)*(active_set_.size()+2) +
                                (active_set_.size()+1)] = value;
            padded_similarities[(active_set_.size()+1)*(active_set_.size()+2) +
                                (i+1)] = value;
          }
          double value = static_cast<double>(
              CountCommonValues(vertex, configuration, configuration));
          padded_similarities[(active_set_.size()+1)*(active_set_.size()+2) +
                              (active_set_.size()+1)] = value;

          vector<double> eigenvalues(active_set_.size()+2);
          EigenDecompose(&padded_similarities, &eigenvalues);
          int zero_eigenvalue = -1;
          for (int i = 0; i < active_set_.size()+2; ++i) {
            if (NEARLY_EQ_TOL(eigenvalues[i], 0.0, 1e-9)) {
              if (zero_eigenvalue >= 0) {
                // If this happens, something failed. Maybe a numerical problem
                // may cause this. In that case, just give up, clean the cache
                // and return. Hopefully the next iteration will fix it.
                cout << "Multiple zero eigenvalues: "
                     << eigenvalues[zero_eigenvalue] << " and "
                     << eigenvalues[i] << endl;
                cout << "Warning: Giving up." << endl;
                // Clean the cache.
                for (int j = 0; j < active_set_.size(); ++j) {
                  DeleteConfiguration(active_set_[j]);
                }
                active_set_.clear();
                inverse_A_.clear();
                distribution_.clear();
                return;
              }
              zero_eigenvalue = i;
            }
          }
          assert(zero_eigenvalue >= 0);
          vector<int> configurations_to_remove;
          for (int j = 1; j < active_set_.size()+1; ++j) {
            double value = padded_similarities[zero_eigenvalue*(active_set_.size()+2) + j];
            if (!NEARLY_EQ_TOL(value, 0.0, 1e-9)) {
              configurations_to_remove.push_back(j-1);
            }
          }
          if (opts.verbose > 2) {
            cout << "Pick a configuration to remove (" << configurations_to_remove.size()
                 << " out of " << active_set_.size() << ")." << endl;
          }

          assert(configurations_to_remove.size() >= 1);
          int j = configurations_to_remove[0];

          // Update inv(A).
          InvertAfterRemoval(inverse_A_, active_set_, j);

          // Remove blocking constraint from the active set.
          DeleteConfiguration(active_set_[j]); // Delete configutation.
          active_set_.erase(active_set_.begin() + j);

          singular = !InvertAfterInsertion(vertex, inverse_A_, active_set_, configuration);
          assert(!singular);
        }

        // Insert configuration to active set.
        if (opts.verbose > 2) {
          cout << "Inserted one element to the active set (iteration "
               << iter << ")." << endl;
        }
        active_set_.push_back(configuration);
        changed_active_set = true;
      }      
    } else {
      // Solution has changed from the previous iteration.
      // Look for blocking constraints.
      int blocking = -1;
      bool exist_blocking = false;
      double alpha = 1.0;
      for (int i = 0; i < active_set_.size(); ++i) {
        assert(distribution_[i] >= -1e-12);
        if (z[i] >= distribution_[i]) continue;
        if (z[i] < 0) exist_blocking = true;
        double tmp = distribution_[i] / (distribution_[i] - z[i]);
        if (blocking < 0 || tmp < alpha) {
          alpha = tmp;
          blocking = i;
        }
      }

      if (!exist_blocking) {
        // No blocking constraints.
        assert(!unbounded);
        distribution_ = z;
        alpha = 1.0;
        changed_active_set = false;
      } else {
        if (alpha > 1.0 && !unbounded) alpha = 1.0;
        // Interpolate between factor_posteriors_[i] and z.
        if (alpha == 1.0) {
          distribution_ = z;
        } else {
          for (int i = 0; i < active_set_.size(); ++i) {
            z[i] = (1 - alpha) * distribution_[i] + alpha * z[i];
            distribution_[i] = z[i];
          }
        }

        // Update inv(A).
        InvertAfterRemoval(inverse_A_, active_set_, blocking);

        // Remove blocking constraint from the active set.
        if (opts.verbose > 2) {
          cout << "Removed one element to the active set (iteration "
               << iter << ")." << endl;
        }

        DeleteConfiguration(active_set_[blocking]); // Delete configutation.
        active_set_.erase(active_set_.begin() + blocking);

        z.erase(z.begin() + blocking);
        distribution_.erase(distribution_.begin() + blocking);
        changed_active_set = true;
        for (int i = 0; i < distribution_.size(); ++i) {
          assert(distribution_[i] > -1e-16);
        }
      }
    }
  }

  // Maximum number of iterations reached.
  // Return the best existing solution by computing the variable marginals 
  // from the full distribution stored in z.
  //assert(false);
  ComputeMarginalsFromSparseDistribution(vertex, active_set_,
                                         z,
                                         variable_posteriors,
                                         additional_posteriors); 
  }
  
  
  void InsertionSort(pair<double, int> arr[], int length) {
  int i, j;
  pair<double, int> tmp;

  for (i = 1; i < length; i++) {
    j = i;
    while (j > 0 && arr[j - 1].first > arr[j].first) {
      tmp = arr[j];
      arr[j] = arr[j - 1];
      arr[j - 1] = tmp;
      j--;
    }
  }
}

  
  int project_onto_budget_constraint_cached(vec& x,
                                          int d,
                                          double budget, 
                                          vector<pair<double,int> >& y) {
  int j, k, l, level;
  double s = 0.0;
  double tau = 0.0, tightsum;
  double left, right = -std::numeric_limits<double>::infinity();

  // Load x into a reordered y (the reordering is cached).
  if (y.size() != d) {
    y.resize(d);
    for (j = 0; j < d; j++) {
      s -= x[j];
      y[j].first = -x[j];
      y[j].second = j;
    }
    sort(y.begin(), y.end());
  } else {
    for (j = 0; j < d; j++) {
      s -= x[j];
      y[j].first = -x[y[j].second];
    }
    // If reordering is cached, use a sorting algorithm 
    // which is fast when the vector is almost sorted.
    InsertionSort(&y[0], d);
  }

  tightsum = s;
  s += budget;
  
  k = l = level = 0;
  bool found = false;
  double val_a, val_b;
  while (k < d && l < d) {
    if (level != 0) {
      tau = (s - tightsum) / static_cast<double>(level);
    }
    if (k < d) val_a = y[k].first;
    val_b = 1.0 + y[l].first;
    left = right;
    if (k == d || val_b <= val_a) {
      right = val_b;
    } else {
      right = val_a;
    }
    if ((level == 0 && s == tightsum) || (level != 0 && tau <= right)) {
      // Found the right split-point!
      found = true;
      break;
    }
    if (k == d || val_b <= val_a) {
      tightsum += val_b;
      --level;
      ++l;
    } else {
      tightsum -= val_a;
      ++level;
      ++k;
    }
  }

  if (!found) {
    left = right;
    right = std::numeric_limits<double>::infinity();
  }
      
  for (j = 0; j < d; j++) {
    if (-x[j] >= right) {
      x[j] = 0.0;
    } else if (1.0 - x[j] <= left) {
      x[j] = 1.0;
    } else {
      x[j] += tau;
    }
  }

  return 0;
}

  
  int project_onto_budget_constraint(vec& x, int d, double budget) {
  int j, k, l, level;
  double s = 0.0;
  double tau = 0.0, tightsum;
  double left, right = -std::numeric_limits<double>::infinity();
  vector<double> y(d, 0.0);

  for (j = 0; j < d; j++) {
    s -= x[j];
    y[j] = -x[j];
  }
  sort(y.begin(), y.end());
  tightsum = s;
  s += budget;
  
  k = l = level = 0;
  bool found = false;
  double val_a, val_b;
  while (k < d && l < d) {
    if (level != 0) {
      tau = (s - tightsum) / static_cast<double>(level);
    }
    if (k < d) val_a = y[k];
    val_b = 1.0 + y[l];
    left = right;
    if (k == d || val_b <= val_a) {
      right = val_b;
    } else {
      right = val_a;
    }
    if ((level == 0 && s == tightsum) || (level != 0 && tau <= right)) {
      // Found the right split-point!
      found = true;
      break;
    }
    if (k == d || val_b <= val_a) {
      tightsum += val_b;
      --level;
      ++l;
    } else {
      tightsum -= val_a;
      ++level;
      ++k;
    }
  }

  if (!found) {
    left = right;
    right = std::numeric_limits<double>::infinity();
  }
      
  for (j = 0; j < d; j++) {
    if (-x[j] >= right) {
      x[j] = 0.0;
    } else if (1.0 - x[j] <= left) {
      x[j] = 1.0;
    } else {
      x[j] += tau;
    }
  }

  return 0;
}
  
  // Solve the QP subproblem for budget factor.
  // TODO Enable caching
void SolveQP_budget(vertex_type& vertex,const gather_type& total,
            vec& variable_posteriors, vec& additional_posteriors){
            
  vertex_data& vdata =  vertex.data();
  vec variable_log_potentials = total.neighbor_distribution + total.messages;                  
  vector<pair<double,int> > last_sort_;
  for (int f = 0; f < variable_log_potentials.size(); ++f) {
    variable_posteriors[f] = variable_log_potentials[f];
    if (variable_posteriors[f] < 0.0) {
      variable_posteriors[f] = 0.0;
    } else if (variable_posteriors[f] > 1.0) {
      variable_posteriors[f] = 1.0;
    }
  }

  double s = 0.0;
  for (int f = 0; f < vdata.nvars; ++f) {
    s += variable_posteriors[f];
  }

  if (s > static_cast<double>(vdata.budget)) {
    for (int f = 0; f < variable_log_potentials.size(); ++f) {
      variable_posteriors[f] = variable_log_potentials[f];
    }
    project_onto_budget_constraint(variable_posteriors, 
                                          variable_log_potentials.size(), 
                                          static_cast<double>(vdata.budget));
  }
}

// Finds best configuration of budget factors
void SolveMAP_budget(vertex_type& vertex,const gather_type& total,
            vec& variable_posteriors, vec& additional_posteriors, double& value) {
 
   vertex_data& vdata = vertex.data();
  // Create a local copy of the log potentials.
   vec log_potentials(total.messages); 
  double valaux;
 
  value = 0.0;
  
  int num_active = 0;
  double sum = 0.0;

  for (int f = 0; f < vdata.nvars; ++f) {
    valaux = log_potentials[f];
    if (valaux < 0.0) {
      variable_posteriors[f] =  0.0;
    } else {
      sum += valaux;
      variable_posteriors[f] = 1.0;
    }
    ++num_active;
  }
  if (num_active > vdata.budget) {
    vector<pair<double,int> > scores(vdata.nvars);
    for (int f = 0; f < vdata.nvars; ++f) {
      scores[f].first = -log_potentials[f];
      scores[f].second = f;
    }

    sort(scores.begin(), scores.end());
    num_active = 0;
    sum = 0.0;
    for (int k = 0; k < vdata.budget; ++k) {
      valaux = -scores[k].first;
      if (valaux < 0.0) break;
      int f = scores[k].second;
      variable_posteriors[f] = 1.0;
      sum += valaux;
      ++num_active;      
    }

    for (int k = num_active; k < vdata.nvars; ++k) {
      int f = scores[k].second;
      variable_posteriors[f] = 0.0;
    }   
  }
  
  value += sum;
  
  
}

// Finds best configuration for dense factors
void SolveMAP_dense(vertex_type& vertex,const gather_type& total,
            vec& variable_posteriors, vec& additional_posteriors, double& value){
     vertex_data& vdata = vertex.data();
     vec beliefs = vdata.potentials;         
     int num_configurations = vdata.potentials.size();
     for (int index_configuration = 0;
           index_configuration < num_configurations;
            ++index_configuration) {
        vector<int> states(vdata.nvars, -1);
        get_configuration_states(vertex, index_configuration, &states);
        int offset = 0;
        for (int k = 0; k < vdata.nvars; ++k) {
             beliefs[index_configuration] += total.messages[offset + states[k]];
             offset += vdata.cards[k];}
    } 
            
        value = beliefs.maxCoeff();
 
 }

// Finds beliefs using dense and budget factors
void compute_beliefs(vertex_type& vertex,const gather_type& total,
            vec& variable_posteriors, vec& additional_posteriors){
switch(vertex.data().factor_type){

case 0: SolveQP_dense(vertex,total, variable_posteriors, additional_posteriors);
        break;
case 1: SolveQP_budget(vertex,total, variable_posteriors, additional_posteriors);

}
};

// General solveMAP function
void SolveMAP(vertex_type& vertex,const gather_type& total,
            vec& variable_posteriors, vec& additional_posteriors, double& value){
switch(vertex.data().factor_type){

case 0: SolveMAP_dense(vertex,total, variable_posteriors, additional_posteriors, value);
        break;
case 1: SolveMAP_budget(vertex,total, variable_posteriors, additional_posteriors, value);
  }
 };

};
/* end of ad3_vertex_program*/



////////////////////////////////////////////////////////////////////////////////
// This class implements the Bethe-ADMM as  described in:
//
//   Qiang Fu, Huahua Wang and Arindam Banerjee.
// "Bethe-ADMM for Tree Decomposition based Parallel MAP Inference"
//  Conference on Uncertainty in Artificial Intelligence (UAI), 2013
//
//////////////////////////////////////////////////////////////////////////////// 

struct bethe_admm_vertex_program:public admm_vertex_program {

/* compute_grad_phi computes the gradient of bethe entropy for the factor */
  
       void compute_grad_phi(vertex_type& vertex,vec& unary_beliefs, 
                   vec& factor_beliefs, vec& unary_grad,vec& factor_grad) {

             vertex_data& vdata = vertex.data();  
            // computation for variable beliefs  
             for(int i=0; i< unary_beliefs.size(); i++){
                 unary_grad[i] *= EXP * (unary_beliefs[i]);
             }
             // computation for factor beliefs
             vector<int> states(vdata.nvars);
             for(int i=0; i< factor_beliefs.size(); i++){
                 factor_grad[i] *= factor_beliefs[i] / EXP;
                 get_configuration_states(vertex, i, &states);
                 int offset =0;
                 for(int j=0; j< vdata.nvars; j++){
                     factor_grad[i] /= unary_beliefs[offset + states[j]];
                     offset += vdata.cards[j];
                 }
             }
        }


/* run_bp computes marginal beliefs using sum-product belief propagation */
       void run_bp(vertex_type& vertex, vec& unary_pots, vec& factor_pots, 
             vec& unary_margs, vec& factor_margs) {

            vertex_data& vdata = vertex.data();        
            unary_margs.resize(unary_pots.size());
            factor_margs = factor_pots;
            vector<int> states(vdata.nvars, -1);
           // computing messages
           for(int i=0; i < vdata.nvars; i++) {
               
               vec messages = factor_pots;
               for(int j=0; j < factor_pots.size(); j++) {
                  get_configuration_states(vertex, j, &states);
                  int offset = 0;
                  for(int k = 0; k < vdata.nvars; k++) {
                     if( k != i)
                     messages[j] *= unary_pots[offset + states[k]];
                     offset += vdata.cards[k];
                  }
               }
               vector<double>  marg_messages(vdata.cards[i], 0);
               for(int j=0; j < factor_pots.size(); j++)  {
                   get_configuration_states(vertex, j, &states);
                   marg_messages[states[i]] += messages[j];
                   }
               int var_offset =0;
               for(int j=0; j < i; j++) {          
                   var_offset += vdata.cards[j];
               }
               double sum = 0;
               // computing marginal beliefs for variables
               for(int j=0; j < marg_messages.size(); j++) { 
                  unary_margs[var_offset + j] = marg_messages[j] 
                                           * unary_pots[var_offset + j];
                  sum += unary_margs[var_offset + j];
               }
               for(int j=0; j < marg_messages.size(); j++) { 
                   unary_margs[var_offset + j] /= sum ; 
               }       
          }
          // compuitng factor beliefs
          double fact_sum = 0;
          for(int i=0; i < factor_pots.size(); i++) {
              get_configuration_states(vertex, i, &states);
              int offset = 0;
              for(int j =0; j < vdata.nvars; j++) {
                 factor_margs[i] *= unary_pots[offset + states[j]];
                 offset += vdata.cards[j];
              }
              fact_sum += factor_margs[i];
          }
        
          for(int i=0; i < factor_pots.size(); i++) {
              factor_margs[i] /= fact_sum;
          }
         
       }


/* adjust_beliefs prevents overflow/ underflow of belief variable */
       void adjust_beliefs(vertex_type& vertex){
            
            vertex_data& vdata = vertex.data();
            for(int i=0; i< vdata.beliefs.size(); i++){
                if(vdata.beliefs[i] < 10e-100)
                   vdata.beliefs[i] = 10e-100;
            }  
            for(int i=0; i< vdata.factor_beliefs.size(); i++){
                if(vdata.factor_beliefs[i] < 10e-100)
                   vdata.factor_beliefs[i] = 10e-100;
            } 
       }
       
/* exponentiates potentials for bp. TODO: use faster approximation of pow */
       void exponentiate(vec& potential_vector){
       
            for(int i=0; i< potential_vector.size(); i++){
                potential_vector[i] = pow(EXP, potential_vector[i]);
            }
       }
       

/* solves QP for factor vertices using bp */
       void compute_beliefs(vertex_type& vertex,const gather_type& total,
            vec& variable_posteriors, vec& additional_posteriors){

             vertex_data& vdata = vertex.data();
             vec unary_eta, factor_eta;
             // computing eta  
             factor_eta = (vdata.potentials)/opts.alpha; 
             unary_eta = total.messages + opts.step_size * (total.neighbor_distribution - vdata.beliefs); 
             unary_eta = (unary_eta)/opts.alpha;        
             exponentiate(unary_eta);
             exponentiate(factor_eta);
             vec unary_grad, factor_grad;
             unary_grad.resize(unary_eta.size());
	     factor_grad.resize(factor_eta.size());
             compute_grad_phi(vertex, vdata.beliefs, vdata.factor_beliefs, unary_eta, factor_eta);
             //running bp on eta
             run_bp(vertex, unary_eta, factor_eta, vdata.beliefs, vdata.factor_beliefs); 
             //adjust beliefs for overflow/underflow          
             adjust_beliefs(vertex);             
        };
 /* solves MAP for factor vertices */      
        void SolveMAP(vertex_type& vertex,const gather_type& total,
            vec& variable_posteriors, vec& additional_posteriors, double& value){
             vertex_data& vdata = vertex.data();
             vec beliefs = vdata.potentials;
             vector<int> states(vdata.nvars);
             for(int i=0; i< vdata.potentials.size(); i++) {
                 get_configuration_states(vertex, i, &states);
                 int offset = 0;
                 for(int j=0; j< vdata.nvars; j++){
                     beliefs[i] += total.messages[offset + states[j]];
                     offset += vdata.cards[j];
                 }
             }
             value = beliefs.maxCoeff();  
         };

};
/* end of  bethe_admm_vertex_program */

#endif
